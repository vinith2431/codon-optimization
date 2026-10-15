from pathlib import Path
import csv
import json
import statistics
import sys
import time

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from optimizers.genetic_algorithm import optimize_ga
from core.cai import calculate_cai
from core.codon_pair import codon_pair_score
from core.sequence import get_codons, gc_content
from core.verification import verify_translation
from constraints.homopolymer import find_homopolymers

FASTA = ROOT / "data/genes/fasta/gene08_lacz.fasta"
PAIR_FILE = ROOT / "data/codon_pair/ecoli_codon_pair_table.json"
OUT = ROOT / "results/csv"
SEEDS = [42, 43, 44, 45, 46]
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Shared settings for a fair comparison
COMMON = dict(
    organism="ecoli",
    gc_min=0.40,
    gc_max=0.60,
    target_gc=0.55,
    max_homopolymer=6,
    population_size=50,
    generations=50,
    mutation_rate=0.01,
    crossover_rate=0.80,
    return_history=True,
    device=DEVICE,
)

CONFIGS = [
    {"configuration": "A_baseline"},
    {"configuration": "B_cbp_only", "cbp": True},
    {"configuration": "C_kmer_5prime", "sequence_constraints": True},
    {"configuration": "D_combined", "cbp": True, "sequence_constraints": True},
]


def read_fasta(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing FASTA: {path}")
    seq = "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()
    if not seq or len(seq) % 3 or not set(seq) <= set("ACGT"):
        raise ValueError("FASTA must contain valid complete-codon DNA.")
    return seq


def load_pair_table():
    raw = json.loads(PAIR_FILE.read_text(encoding="utf-8"))
    table = {
        tuple(key.split("_")): float(value)
        for key, value in raw.items()
        if len(key.split("_")) == 2
    }
    if not table:
        raise ValueError(f"Invalid or empty codon-pair table: {PAIR_FILE}")
    return table


def kmer_stats(seq, k=6):
    counts = {}
    total = len(seq) - k + 1
    for i in range(max(0, total)):
        key = seq[i:i + k]
        counts[key] = counts.get(key, 0) + 1
    return (max(counts.values()) / total if total > 0 else 0.0,
            len(counts))


def run(original, pair_table, config, seed):
    kwargs = dict(COMMON, seed=seed)

    if config.get("cbp"):
        kwargs.update(
            codon_pair_table=pair_table,
            codon_pair_weight=5.0,
        )

    if config.get("sequence_constraints"):
        kwargs.update(
            kmer_size=6,
            max_kmer_frequency=0.01,
            kmer_weight=20.0,
            five_prime_window=60,
            five_prime_gc_min=0.40,
            five_prime_gc_max=0.50,
            five_prime_gc_weight=20.0,
        )

    start = time.perf_counter()
    result = optimize_ga(original, **kwargs)
    elapsed = time.perf_counter() - start
    seq = result["sequence"]

    protein_ok = verify_translation(original, seq)["valid"]
    gc = gc_content(seq)
    five_gc = gc_content(seq[:60])
    max_kmer, unique_kmers = kmer_stats(seq)
    pair_score = codon_pair_score(seq, pair_table)
    homopolymers = len(find_homopolymers(seq, 6))
    old_codons, new_codons = get_codons(original), get_codons(seq)
    edits = sum(a != b for a, b in zip(old_codons, new_codons))
    edit_rate = edits / len(old_codons)

    shared_ok = (
        protein_ok and len(seq) == len(original)
        and 0.40 <= gc <= 0.60 and homopolymers == 0
    )
    kmer_ok = max_kmer <= 0.01
    five_gc_ok = 0.40 <= five_gc <= 0.50
    active_ok = shared_ok and (
        not config.get("sequence_constraints")
        or (kmer_ok and five_gc_ok)
    )

    return {
        "configuration": config["configuration"],
        "seed": seed,
        "fitness": result["history"][-1]["best_fitness"],
        "cai": calculate_cai(seq, "ecoli"),
        "global_gc": gc,
        "five_prime_gc": five_gc,
        "codon_pair_score": pair_score,
        "max_kmer_frequency": max_kmer,
        "unique_kmers": unique_kmers,
        "edit_rate": edit_rate,
        "runtime_seconds": elapsed,
        "protein_preserved": protein_ok,
        "shared_constraints_pass": shared_ok,
        "kmer_threshold_pass": kmer_ok,
        "five_prime_gc_pass": five_gc_ok,
        "active_constraints_pass": active_ok,
        "optimized_sequence": seq,
    }


def main():
    original = read_fasta(FASTA)
    pair_table = load_pair_table()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []

    print(f"LACZ length: {len(original)} bp | device: {DEVICE.upper()}")
    print(f"Configurations: {len(CONFIGS)} | seeds each: {len(SEEDS)}")

    for config in CONFIGS:
        for seed in SEEDS:
            print(f"[{len(rows)+1}/20] {config['configuration']} seed={seed}")
            row = run(original, pair_table, config, seed)
            rows.append(row)
            print(
                f"  CAI={row['cai']:.5f} "
                f"CBP={row['codon_pair_score']:.5f} "
                f"max-kmer={row['max_kmer_frequency']:.5f} "
                f"5'-GC={row['five_prime_gc']:.3f} "
                f"checks={row['active_constraints_pass']} "
                f"time={row['runtime_seconds']:.2f}s"
            )

    detail_file = OUT / "combined_constraints.csv"
    with detail_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for config in CONFIGS:
        group = [r for r in rows if r["configuration"] == config["configuration"]]
        summary.append({
            "configuration": config["configuration"],
            "runs": len(group),
            "mean_fitness": statistics.mean(r["fitness"] for r in group),
            "std_fitness": statistics.stdev(r["fitness"] for r in group),
            "mean_cai": statistics.mean(r["cai"] for r in group),
            "mean_codon_pair_score": statistics.mean(r["codon_pair_score"] for r in group),
            "mean_max_kmer_frequency": statistics.mean(r["max_kmer_frequency"] for r in group),
            "mean_five_prime_gc": statistics.mean(r["five_prime_gc"] for r in group),
            "mean_global_gc": statistics.mean(r["global_gc"] for r in group),
            "mean_edit_rate": statistics.mean(r["edit_rate"] for r in group),
            "mean_runtime_seconds": statistics.mean(r["runtime_seconds"] for r in group),
            "protein_preservation_rate": statistics.mean(float(r["protein_preserved"]) for r in group),
            "shared_constraint_pass_rate": statistics.mean(float(r["shared_constraints_pass"]) for r in group),
            "active_constraint_pass_rate": statistics.mean(float(r["active_constraints_pass"]) for r in group),
        })

    summary_file = OUT / "combined_constraints_summary.csv"
    with summary_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    print("\nBenchmark complete.")
    print("Detailed results:", detail_file)
    print("Summary:", summary_file)
    print("Protein preservation:", sum(r["protein_preserved"] for r in rows), "/20")
    print("Shared constraints passed:", sum(r["shared_constraints_pass"] for r in rows), "/20")
    print("Active constraints passed:", sum(r["active_constraints_pass"] for r in rows), "/20")


if __name__ == "__main__":
    main()
