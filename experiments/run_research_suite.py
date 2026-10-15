"""
Integrated Research Experiment Suite.

Orchestrates all codon optimization benchmarks:
1. Fair Optimizer Comparison (Greedy, MILP, GA single-seed, GA multi-start, NSGA-II)
2. Constraint Stringency Sensitivity (Repeat threshold, 5' GC, Local GC)
3. Fitness Weight Sensitivity (CAI, Codon-Pair Bias, 5' GC, K-mer complexity)
4. Multi-Objective Pareto Benchmark (NSGA-II vs. Weighted-GA)
5. Independent Candidate Validation and Validation Failure Logging
6. Automated Research Plot Generation

Outputs:
    results/csv/optimizer_comparison.csv
    results/csv/optimizer_summary.csv
    results/csv/constraint_sensitivity.csv
    results/csv/constraint_summary.csv
    results/csv/fitness_weight_sensitivity.csv
    results/csv/fitness_weight_summary.csv
    results/csv/pareto_candidates.csv
    results/csv/pareto_summary.csv
    results/csv/validation_failures.csv
    results/plots/*.png (8 research figures)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import experiments.run_optimizer_comparison as opt_comp
import experiments.run_constraint_sensitivity as const_sens
import experiments.run_fitness_weight_sensitivity as weight_sens
import experiments.run_pareto_benchmark as pareto_bench
import experiments.generate_plots as gen_plots

OUT_DIR = ROOT / "results" / "csv"


def max_kmer_frequency(sequence: str, k: int = 6) -> float:
    if k <= 0:
        raise ValueError("k must be positive")
    if len(sequence) < k:
        return 0.0
    from collections import Counter
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    counts = Counter(kmers)
    return max(counts.values()) / len(kmers)


def five_prime_gc(sequence: str, window: int = 60) -> float:
    if window <= 0:
        raise ValueError("window must be positive")
    from core.sequence import gc_content
    return gc_content(sequence[: min(window, len(sequence))])


def longest_homopolymer(sequence: str) -> int:
    if not sequence:
        return 0
    best = run = 1
    for a, b in zip(sequence, sequence[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def read_fasta(path: Path) -> str:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing benchmark FASTA: {path}")
    seq = "".join(
        x.strip()
        for x in path.read_text(encoding="utf-8").splitlines()
        if x.strip() and not x.lstrip().startswith(">")
    ).upper().replace("U", "T")
    if not seq:
        raise ValueError(f"Empty FASTA: {path}")
    if set(seq) - set("ACGT"):
        raise ValueError(f"Non-ACGT sequence in {path}")
    if len(seq) % 3:
        raise ValueError(f"Coding length must be divisible by 3: {path}")
    return seq


def collect_validation_failures():
    """Aggregate all infeasible candidate records from detailed benchmarks into validation_failures.csv."""
    import csv

    failure_rows = []
    source_files = [
        OUT_DIR / "optimizer_comparison.csv",
        OUT_DIR / "constraint_sensitivity.csv",
        OUT_DIR / "fitness_weight_sensitivity.csv",
        OUT_DIR / "pareto_candidates.csv",
    ]

    for sf in source_files:
        if not sf.exists():
            continue
        with open(sf, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                is_feasible = row.get("feasible")
                if is_feasible in ("False", False, "0", 0):
                    failure_rows.append({
                        "source_experiment": sf.stem,
                        "gene": row.get("gene", ""),
                        "method_or_config": row.get("method", row.get("config_name", "")),
                        "seed": row.get("seed", ""),
                        "optimized_cai": row.get("optimized_cai", row.get("cai", "")),
                        "protein_preserved": row.get("protein_preserved", ""),
                        "global_gc_pass": row.get("global_gc_pass", ""),
                        "five_prime_gc_pass": row.get("five_prime_gc_pass", ""),
                        "local_gc_pass": row.get("local_gc_pass", ""),
                        "homopolymer_pass": row.get("homopolymer_pass", ""),
                        "repeat_pass": row.get("repeat_pass", ""),
                        "restriction_pass": row.get("restriction_pass", ""),
                        "failure_reasons": row.get("failure_reasons", ""),
                        "error": row.get("error", ""),
                    })

    fail_path = OUT_DIR / "validation_failures.csv"
    if failure_rows:
        fields = list(dict.fromkeys(k for r in failure_rows for k in r.keys()))
        with open(fail_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(failure_rows)
        print(f"Logged {len(failure_rows)} failure records to: {fail_path}")
    else:
        fail_path.write_text("No validation failures recorded.\n", encoding="utf-8")
        print("No validation failures recorded.")


def main():
    parser = argparse.ArgumentParser(description="Integrated Codon Optimization Research Suite")
    parser.add_argument("--smoke", action="store_true", help="Run fast verification smoke test across all modules")
    args = parser.parse_args()

    mode_str = "SMOKE TEST" if args.smoke else "FULL RESEARCH SUITE"
    start_total = time.perf_counter()

    print("*" * 80)
    print(f"STARTING CODON OPTIMIZATION RESEARCH SUITE [{mode_str}]")
    print("*" * 80)

    # 1. Fair Optimizer Comparison
    print("\n>>> STEP 1/5: Running Fair Optimizer Comparison...")
    sys.argv = ["run_optimizer_comparison.py"] + (["--smoke"] if args.smoke else [])
    opt_comp.main()

    # 2. Constraint Stringency Sensitivity
    print("\n>>> STEP 2/5: Running Constraint Sensitivity Study...")
    sys.argv = ["run_constraint_sensitivity.py"] + (["--smoke"] if args.smoke else [])
    const_sens.main()

    # 3. Fitness Weight Sensitivity
    print("\n>>> STEP 3/5: Running Fitness Weight Sensitivity Study...")
    sys.argv = ["run_fitness_weight_sensitivity.py"] + (["--smoke"] if args.smoke else [])
    weight_sens.main()

    # 4. Multi-Objective Pareto Benchmark
    print("\n>>> STEP 4/5: Running Multi-Objective Pareto Benchmark (NSGA-II)...")
    sys.argv = ["run_pareto_benchmark.py"] + (["--smoke"] if args.smoke else [])
    pareto_bench.main()

    # 5. Validation Failures & Plot Generation
    print("\n>>> STEP 5/5: Aggregating Validation Failures and Generating Figures...")
    collect_validation_failures()
    gen_plots.main()

    elapsed_total = time.perf_counter() - start_total
    print("\n" + "*" * 80)
    print(f"RESEARCH SUITE COMPLETED SUCCESSFULLY IN {elapsed_total:.2f}s")
    print(f"Outputs written to: {OUT_DIR} and {ROOT / 'results' / 'plots'}")
    print("*" * 80)


if __name__ == "__main__":
    main()
