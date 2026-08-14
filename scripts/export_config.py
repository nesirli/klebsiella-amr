#!/usr/bin/env python3
"""Export config.yaml to a Makefile-include file."""
import argparse
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser(description="Export config.yaml to config.mk")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    c = yaml.safe_load(open(args.config))

    lines = [
        f"METADATA := {c['metadata']}",
        f"DELIMITER := {c['delimiter']}",
        f"ANTIBIOTICS := {' '.join(c['antibiotics'])}",
        f"TRAIN_CUTOFF := {c['splits']['train_cutoff']}",
        f"TEST_YEARS := {' '.join(map(str, c['splits']['test_years']))}",
        f"MAX_SAMPLES := {c.get('max_samples', -1)}",
        f"GENOME_SIZE := {c['genome']['size']}",
        f"TARGET_COVERAGE := {c['genome']['target_coverage']}",
        f"KRAKEN_DB := {c['reference']['kraken2_db']}",
        f"AMRFINDER_DB := {c['reference']['amrfinder_db']}",
        f"FASTP_THREADS := {c['resources']['fastp_threads']}",
        f"KRAKEN_THREADS := {c['resources']['kraken2_threads']}",
        f"QUAST_THREADS := {c['resources']['quast_threads']}",
        f"SPADES_THREADS := {c['resources']['spades_threads']}",
        f"SPADES_MEMORY := {c['resources']['spades_memory']}",
        f"AMRFINDER_THREADS := {c['resources']['amrfinder_threads']}",
    ]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
