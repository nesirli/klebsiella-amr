"""Streamlit front end: upload an assembled genome, get a resistance profile.

Shows the rule-based call and the model calls side by side rather than picking
one, because neither wins everywhere: the rule is better for meropenem
(balanced accuracy 0.93 vs 0.86) and useless for ciprofloxacin (specificity
0.02, since K. pneumoniae intrinsically carries quinolone efflux genes). A
disagreement between them is the signal to look closer, so it is highlighted
rather than resolved.

Run with:  make app        (or: streamlit run app.py)
Login:     .streamlit/secrets.toml -> username/password, or
           $AMR_APP_USERNAME / $AMR_APP_PASSWORD
"""
import hmac
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
from streamlit.errors import StreamlitSecretNotFoundError

sys.path.insert(0, str(Path(__file__).parent))
import predict as P  # noqa: E402

# Baked into the image next to this file; overridable for local development
# against the pipeline's own results/ tree.
HERE = Path(__file__).parent
MANIFEST = os.environ.get("MANIFEST", str(HERE / "artifacts" / "manifest.json"))
SPECIES_SKETCH = os.environ.get(
    "SPECIES_SKETCH", str(HERE / "artifacts" / "kpneumoniae_ref.msh"))
AMRFINDER_DB = os.environ.get("AMRFINDER_DB", "reference/amrfinderplus")
DEMO_FASTA = os.environ.get(
    "DEMO_FASTA", str(HERE / "artifacts" / "demo_SRR30762169.fasta"))

st.set_page_config(page_title="Klebsiella AMR prediction", page_icon="🧬",
                   layout="wide")


def check_login():
    """One shared username and password. A gate for a small trusted group, not
    real authentication: no accounts, no lockout, no audit trail. Serve it over
    HTTPS, and put a proper proxy in front before exposing it more widely.

    Both are compared with compare_digest so a wrong guess takes the same time
    as a right one.
    """
    if st.session_state.get("authenticated"):
        return True

    # st.secrets.get raises StreamlitSecretNotFoundError when no secrets.toml
    # exists at all (which is the case on Railway, where the password comes
    # from the environment). Its .get default only covers missing keys, not a
    # missing file.
    def secret(key, default):
        try:
            return st.secrets.get(key, default)
        except StreamlitSecretNotFoundError:
            return default

    expected_user = secret("username", os.environ.get("AMR_APP_USERNAME", "demo"))
    expected_password = secret("password", os.environ.get("AMR_APP_PASSWORD", ""))

    st.title("🧬 Klebsiella AMR prediction")
    if not expected_password:
        st.error("No password configured. Set `AMR_APP_PASSWORD` in the "
                 "environment, or `password` in `.streamlit/secrets.toml`.")
        return False

    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    if username or password:
        ok_user = hmac.compare_digest(username, expected_user)
        ok_password = hmac.compare_digest(password, expected_password)
        if ok_user and ok_password:
            st.session_state["authenticated"] = True
            st.rerun()
        elif username and password:
            st.error("Incorrect username or password.")
    return False


@st.cache_resource
def load_manifest():
    return P.load_manifest(MANIFEST)


def call_badge(call):
    return "🔴 Resistant" if call == "R" else "🟢 Susceptible"


def main():
    st.title("🧬 Klebsiella AMR prediction")
    st.caption("Research use only. Not a diagnostic tool and not for guiding "
               "treatment.")

    if not Path(MANIFEST).exists():
        st.error(f"Missing `{MANIFEST}`. Run `make app-artifacts` first.")
        return
    manifest = load_manifest()

    with st.sidebar:
        st.subheader("Models")
        st.write(f"Trained on **{manifest['n_genes']} genes**")
        st.write(f"Manifest built {manifest['created']}")
        for name, per_abx in manifest["models"].items():
            aucs = [v["test_roc_auc"] for v in per_abx.values() if v["test_roc_auc"]]
            st.write(f"- **{name}** — mean test ROC AUC {sum(aucs)/len(aucs):.2f}")
        st.caption("Uploaded files are written to a temporary directory and "
                   "deleted when the analysis finishes.")

    uploaded = st.file_uploader(
        "Assembled genome (FASTA)", type=["fasta", "fa", "fna"],
        help="Contigs from an assembler such as SPAdes. Not raw reads.")
    use_demo = st.button(
        "Try a demo genome", disabled=uploaded is not None,
        help="Analyse the bundled demo isolate instead of uploading a file.")
    if not uploaded and not use_demo:
        st.info("Upload an assembled *Klebsiella pneumoniae* genome, or try "
                "the demo genome.")
        return
    if use_demo and not Path(DEMO_FASTA).exists():
        st.error(f"Demo genome missing from the image: `{DEMO_FASTA}`")
        return
    if use_demo:
        st.caption(f"Demo isolate: **{Path(DEMO_FASTA).name}** — a "
                   "multi-resistant test isolate from the training split's "
                   "sister dataset.")

    with tempfile.TemporaryDirectory() as tmp:
        if uploaded:
            fasta = Path(tmp) / uploaded.name
            fasta.write_bytes(uploaded.getbuffer())
        else:
            fasta = Path(DEMO_FASTA)

        # The models only know K. pneumoniae. Anything else would still get a
        # confident answer, drawn from the wrong biology -- so check the
        # organism before spending 30 seconds in AMRFinder.
        try:
            contigs, bases, gc = P.read_fasta_stats(fasta)
        except ValueError as exc:
            st.error(f"Cannot read this file: {exc}")
            st.info("Upload assembled contigs in FASTA format — not sequencing "
                    "reads (FASTQ) and not protein sequences.")
            return

        if not P.MIN_GENOME_BP <= bases <= P.MAX_GENOME_BP:
            st.error(f"This is {bases/1e6:.2f} Mb across {contigs} contigs, "
                     "which is not a bacterial genome. A *Klebsiella "
                     "pneumoniae* assembly is about 5.0–5.9 Mb.")
            return

        with st.spinner("Checking species…"):
            is_kleb, distance, species_msg = P.check_species(fasta, SPECIES_SKETCH)

        if not is_kleb:
            st.error(f"**Rejected — this does not look like *Klebsiella "
                     f"pneumoniae*.** {species_msg}")
            st.info("The models were trained only on *K. pneumoniae*. Running "
                    "them on another organism would return a confident answer "
                    "based on the wrong biology, so the upload is refused "
                    "rather than scored.")
            with st.expander("What was measured"):
                st.markdown(
                    f"- Assembly: **{bases/1e6:.2f} Mb** in **{contigs}** "
                    f"contigs, **{gc:.1%}** GC\n"
                    f"- Mash distance to a 20-genome *K. pneumoniae* reference "
                    f"sketch: **{distance:.3f}**" if distance is not None else
                    f"- Assembly: **{bases/1e6:.2f} Mb**, **{gc:.1%}** GC\n"
                    f"- Species check could not run")
            return

        st.caption(f"Species check passed — {species_msg}. "
                   f"{bases/1e6:.2f} Mb, {contigs} contigs, {gc:.1%} GC.")

        with st.spinner("Running AMRFinderPlus (about 30 seconds)…"):
            try:
                tsv = P.run_amrfinder(fasta, tmp, AMRFINDER_DB)
            except Exception as exc:  # noqa: BLE001 - surfaced to the user
                st.error(f"AMRFinderPlus failed: {exc}")
                return
            symbols, subclasses, table = P.parse_amr(tsv)

        results, matched, unknown = P.predict(manifest, symbols, subclasses)

    st.success(f"Found **{len(symbols)}** resistance genes; "
               f"**{len(matched)}** are features the models know.")
    if unknown:
        st.warning(f"{len(unknown)} gene(s) not seen in training were ignored: "
                   f"{', '.join(unknown[:8])}"
                   f"{'…' if len(unknown) > 8 else ''}. Many of these can mean "
                   "the AMRFinder database has moved on since training.")

    st.subheader("Predicted resistance")
    rows = []
    for entry in results:
        row = {"Antibiotic": entry["antibiotic"].title(),
               "AMRFinder rule": call_badge(entry["rule_call"])}
        for name, m in entry["models"].items():
            row[name] = f"{call_badge(m['call'])}  ({m['probability']:.2f})"
        calls = {entry["rule_call"], *(m["call"] for m in entry["models"].values())}
        row["Agreement"] = "✅ agree" if len(calls) == 1 else "⚠️ disagree"
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    if any(len({e["rule_call"], *(m["call"] for m in e["models"].values())}) > 1
           for e in results):
        st.info("⚠️ Where the rule and the models disagree, treat the result as "
                "uncertain and inspect the genes below.")

    for entry in results:
        if entry["low_confidence"]:
            st.warning(
                f"**{entry['antibiotic'].title()}**: unreliable for isolates "
                "collected later than the training set (test ROC AUC 0.62, "
                "versus 0.94 within the training distribution). The training "
                "years are 47% resistant and the test years 85%, so the model "
                "does not transfer across eras. Treat this row as indicative "
                "only.")

    with st.expander(f"Resistance genes found ({len(symbols)})"):
        st.dataframe(table, use_container_width=True, hide_index=True)

    with st.expander("How to read this"):
        st.markdown(
            "- **AMRFinder rule** — resistant if any gene of that drug class is "
            "present. Transparent, and the better predictor for meropenem "
            "(balanced accuracy 0.93 vs 0.86 for the models).\n"
            "- **Models** — XGBoost and LightGBM over gene presence/absence, "
            "each using a threshold fitted on the training split rather than "
            "0.5. Clearly better where the rule collapses: for ciprofloxacin "
            "the rule has specificity 0.02, because nearly every *K. "
            "pneumoniae* carries intrinsic quinolone efflux genes.\n"
            "- **Amikacin** is the most reliable (ROC AUC ~1.0); "
            "**ceftazidime** the least (see the warning above)."
        )


if __name__ == "__main__":
    if check_login():
        main()
