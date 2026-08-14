#!/usr/bin/env python3
"""Fine-tune DNABERT-2 on AMR-gene nucleotide sequences, per antibiotic.

Unlike the tabular models (which see a 0/1 gene presence matrix), this reads
the *nucleotide sequence* of each detected AMR gene (build_sequences.py) and
fine-tunes the pretrained DNABERT-2 genome language model. That lets it key on
allelic variants / point mutations within a gene that a presence flag collapses
away -- the actual reason to bring a sequence model to this problem.

Architecture: each gene sequence is encoded separately through DNABERT-2
(masked-mean pooling over tokens), and a sample is represented by the mean of
its present genes' embeddings, fed to a linear head. Encoding genes separately
(rather than concatenating them into one long string) means NO gene is lost to
truncation -- important because a sample can carry ~16 genes totalling many kb,
far past the model's context window. It also makes occlusion importance exact
and cheap: dropping a gene from the pooled mean needs only a head re-run, no
re-encoding.

NOTE this is DNABERT-2, which tokenizes with BPE (not the k-mer scheme of
DNABERT-1); DNABERT-2's whole contribution was replacing k-mers with BPE.

Produces the same six artifacts as the other models. Runs on CPU or CUDA; the
environment installs the CPU torch wheel, and force_pytorch_attention() keeps
DNABERT off its CUDA-only Triton kernel either way.
"""
import numpy as np
import torch
import torch.nn as nn
from transformers import AutoModel, AutoTokenizer

import common

MODEL_NAME = "dnabert"
BASE_MODEL = "zhihan1996/DNABERT-2-117M"

torch.manual_seed(42)
np.random.seed(42)

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def force_pytorch_attention():
    """Make DNABERT-2 use its portable PyTorch attention, not the Triton
    flash-attention kernel.

    bert_layers.py picks the flash kernel unless its module-global
    `flash_attn_qkvpacked_func` is None (the "Triton unavailable" state). That
    kernel asserts CUDA tensors, so it crashes on CPU -- and a CUDA torch build
    pulls triton in, which would otherwise select it. Nulling the global forces
    the PyTorch path, which runs on both CPU and GPU. Done after model load;
    the class reads the global at forward time. No-op if the model version
    doesn't expose it.
    """
    import sys

    for name, mod in list(sys.modules.items()):
        if name.endswith("bert_layers") and hasattr(mod, "flash_attn_qkvpacked_func"):
            mod.flash_attn_qkvpacked_func = None


class DNABertGeneClassifier(nn.Module):
    """Per-gene DNABERT-2 encoder + linear head over pooled gene embeddings."""

    def __init__(self, encoder, hidden_size):
        super().__init__()
        self.encoder = encoder
        self.head = nn.Linear(hidden_size, 2)

    def encode(self, input_ids, attention_mask):
        """Masked-mean token pooling -> one embedding per input sequence."""
        hidden = self.encoder(input_ids, attention_mask=attention_mask)[0]
        mask = attention_mask.unsqueeze(-1).float()
        return (hidden * mask).sum(1) / mask.sum(1).clamp(min=1.0)


def sample_gene_lists(df, genes):
    """Per sample, the list of (gene, sequence) for its present genes."""
    sub = df[genes].fillna("").astype(str)
    return [[(genes[j], s) for j, s in enumerate(row) if s] for row in sub.to_numpy()]


def flatten(gene_lists):
    """Flatten per-sample gene lists into parallel (text, sample_idx) arrays.
    Samples with no genes get one empty-string entry so they still pool."""
    texts, sample_idx = [], []
    for s, glist in enumerate(gene_lists):
        for _, seq in glist or [(None, "")]:
            texts.append(seq)
            sample_idx.append(s)
    return texts, sample_idx


def encode_texts(model, tokenizer, texts, max_length, batch_size, grad):
    """Encode a flat list of sequences -> [N, H] embedding tensor."""
    ctx = torch.enable_grad() if grad else torch.no_grad()
    chunks = []
    with ctx:
        for i in range(0, len(texts), batch_size):
            enc = tokenizer(texts[i:i + batch_size], padding=True, truncation=True,
                            max_length=max_length, return_tensors="pt").to(DEVICE)
            chunks.append(model.encode(enc["input_ids"], enc["attention_mask"]))
    return torch.cat(chunks, dim=0)


def pool_samples(embeds, sample_idx, n_samples, hidden_size):
    """Mean-pool gene embeddings into one vector per sample."""
    pooled = torch.zeros(n_samples, hidden_size, device=embeds.device)
    counts = torch.zeros(n_samples, device=embeds.device)
    idx = torch.tensor(sample_idx, device=embeds.device)
    pooled.index_add_(0, idx, embeds)
    counts.index_add_(0, idx, torch.ones(len(sample_idx), device=embeds.device))
    return pooled / counts.clamp(min=1.0).unsqueeze(1)


def train(model, tokenizer, gene_lists, y, hp, weights):
    opt = torch.optim.AdamW(model.parameters(), lr=hp["learning_rate"],
                            weight_decay=hp["weight_decay"])
    criterion = nn.CrossEntropyLoss(weight=weights.to(DEVICE))
    hidden = model.head.in_features
    n, bs = len(gene_lists), min(hp["batch_size"], len(gene_lists))
    yt = torch.tensor(y)
    model.train()
    for _ in range(hp["epochs"]):
        perm = torch.randperm(n).tolist()
        for i in range(0, n, bs):
            batch = perm[i:i + bs]
            texts, local_idx = flatten([gene_lists[j] for j in batch])
            embeds = encode_texts(model, tokenizer, texts, hp["max_length"],
                                  hp["encode_batch"], grad=True)
            pooled = pool_samples(embeds, local_idx, len(batch), hidden)
            opt.zero_grad()
            loss = criterion(model.head(pooled), yt[batch].to(DEVICE))
            loss.backward()
            opt.step()


def sample_embeddings(model, tokenizer, gene_lists, hp):
    """No-grad: per-gene embeddings kept per sample (for prediction + occlusion).
    flatten() emits genes in each sample's list order (or one empty entry), so
    per_sample[s][i] aligns with gene_of[s][i]."""
    texts, sample_idx = flatten(gene_lists)
    embeds = encode_texts(model, tokenizer, texts, hp["max_length"],
                          hp["encode_batch"], grad=False)
    per_sample = [[] for _ in gene_lists]
    for k, s in enumerate(sample_idx):
        per_sample[s].append(embeds[k])
    gene_of = [[g for g, _ in glist] if glist else [None] for glist in gene_lists]
    stacked = [torch.stack(v) for v in per_sample]  # every sample has >=1 entry
    return stacked, gene_of


@torch.no_grad()
def proba_from_pooled(model, pooled):
    return torch.softmax(model.head(pooled), dim=1)[:, 1].cpu().numpy()


@torch.no_grad()
def occlusion_importance(model, emb_per_sample, gene_of, genes, y_proba):
    """For each gene: among samples where present, drop it from the pooled mean
    and measure the mean fall in P(R). Only the head re-runs -- exact & cheap."""
    contrib = {g: [] for g in genes}
    for s, (embs, names) in enumerate(zip(emb_per_sample, gene_of)):
        if sum(g is not None for g in names) <= 1:
            continue  # removing the only gene leaves nothing to compare
        for k, g in enumerate(names):
            if g is None:
                continue
            keep = torch.cat([embs[:k], embs[k + 1:]], dim=0).mean(0, keepdim=True)
            contrib[g].append(y_proba[s] - proba_from_pooled(model, keep.to(DEVICE))[0])
    return [float(np.mean(contrib[g])) if contrib[g] else 0.0 for g in genes]


def parse_args():
    p = common.train_parser("Fine-tune DNABERT-2 per antibiotic")
    # Light-touch knobs. Defaults are sized for real data, where a sample can
    # carry ~20 genes of up to ~3 kb: a full-batch encode of every gene with
    # gradients would OOM, so keep samples/step small and cap per-gene length
    # (long genes truncate WITHIN the gene -- far better than dropping whole
    # genes as concatenation would).
    p.add_argument("--epochs", type=int, default=6)
    p.add_argument("--batch-size", type=int, default=2, help="samples per step")
    p.add_argument("--encode-batch", type=int, default=8, help="gene seqs per encoder forward")
    p.add_argument("--max-length", type=int, default=128,
                   help="BPE tokens/gene (~4-5 bp each); caps CPU memory on long genes")
    p.add_argument("--learning-rate", type=float, default=2e-5)
    return p.parse_args()


def main():
    args = parse_args()
    common.prepare_outputs(args)

    hp = {"epochs": args.epochs, "batch_size": args.batch_size,
          "encode_batch": args.encode_batch, "max_length": args.max_length,
          "learning_rate": args.learning_rate, "weight_decay": 1e-2,
          "base_model": BASE_MODEL}

    genes = common.read_genes(args.train_features, args.all_antibiotics)
    str_genes = {g: str for g in genes}
    train_df = common.load_labelled(args.train_features, args.antibiotic, dtype=str_genes)
    test_df = common.load_labelled(args.test_features, args.antibiotic, dtype=str_genes)
    y_train = (train_df[args.antibiotic] == "R").astype(int).to_numpy()
    y_test = (test_df[args.antibiotic] == "R").astype(int).to_numpy()

    base = {"model": MODEL_NAME, "antibiotic": args.antibiotic,
            "n_train": int(len(train_df)), "n_test": int(len(test_df)),
            "n_features": len(genes), "device": DEVICE.type}

    if not common.enough_to_fit(y_train, len(test_df), min_train=2):
        common.write_skipped(args, base, MODEL_NAME)
        return

    train_lists = sample_gene_lists(train_df, genes)
    test_lists = sample_gene_lists(test_df, genes)

    tokenizer = AutoTokenizer.from_pretrained(BASE_MODEL, trust_remote_code=True)
    encoder = AutoModel.from_pretrained(BASE_MODEL, trust_remote_code=True)
    force_pytorch_attention()
    # A sample's ~20 genes are encoded together with gradients; checkpointing
    # recomputes activations in the backward pass instead of storing them,
    # which keeps CPU memory bounded (real genes reach ~3 kb). Best-effort:
    # DNABERT-2's custom model may not support it.
    try:
        encoder.gradient_checkpointing_enable()
    except Exception:
        pass
    hidden_size = encoder.config.hidden_size
    model = DNABertGeneClassifier(encoder, hidden_size).to(DEVICE)

    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)
    train(model, tokenizer, train_lists, y_train, hp, common.class_weights(y_train))

    # Save the *full* fine-tuned model (encoder + head) so the tuned DNABERT-2
    # weights are preserved for inference. This is ~450 MB per antibiotic because
    # the base model is 117M parameters; the alternative is to save only the head
    # and lose the fine-tuned encoder, which defeats the purpose of fine-tuning.
    torch.save({"model_state_dict": model.state_dict(), "genes": genes,
                "hyperparameters": hp, "base_model": BASE_MODEL,
                "hidden_size": hidden_size}, args.model_output)
    common.write_json(args.params_output,
                      {**base, "hyperparameters": hp,
                       "train_class_balance": {"R": n_pos, "S": n_neg}})

    # Test predictions from mean-pooled gene embeddings.
    model.eval()
    emb_per_sample, gene_of = sample_embeddings(model, tokenizer, test_lists, hp)
    pooled = torch.stack([e.mean(0) for e in emb_per_sample]).to(DEVICE)
    y_proba = proba_from_pooled(model, pooled)
    y_pred = (y_proba >= 0.5).astype(int)

    metrics = common.compute_metrics(base, y_test, y_pred, y_proba)
    common.write_json(args.metrics_output, metrics)
    common.write_predictions(args.predictions_output, test_df["run"], y_test, y_pred, y_proba)

    drops = occlusion_importance(model, emb_per_sample, gene_of, genes, y_proba)
    importance = common.write_importance(args.importance_output, genes, drops)
    common.importance_plot(args.importance_plot_output, importance,
                           f"{args.antibiotic} ({MODEL_NAME}) feature importance",
                           "occlusion importance (drop in test P(R) when gene removed)")

    common.print_json(metrics)


if __name__ == "__main__":
    main()
