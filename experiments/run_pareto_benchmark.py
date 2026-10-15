"""
Pareto Multi-Objective Codon Optimization Benchmark.

Compares:
1. True NSGA-II algorithm (optimizers/nsga2.py)
2. Weighted-GA Pareto front approximation (optimizers/pareto.py)

Outputs:
    results/csv/pareto_candidates.csv
    results/csv/pareto_summary.csv
"""

from pathlib import Path
import argparse
import csv
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.cai import calculate_cai
from core.sequence import clean_dna, gc_content
from core.similarity import edit_rate
from core.validator import assess_candidate, DEFAULT_CONSTRAINTS
from optimizers.nsga2 import optimize_nsga2
from optimizers.pareto import optimize_pareto_weighted_ga

DATA_DIR = ROOT / "data" / "genes" / "fasta"
OUT_DIR = ROOT / "results" / "csv"

GENES = {
    "GFP": "gene01_gfp.fasta",
    "Firefly_luciferase": "gene02_luc.fasta",
    "HBB": "gene03_hbb.fasta",
    "INS": "gene04_ins.fasta",
    "LYZ": "gene05_lyz.fasta",
    "ACTB": "gene06_actb.fasta",
    "TP53": "gene07_tp53.fasta",
    "lacZ": "gene08_lacz.fasta",
    "groL": "gene09_grol.fasta",
    "sodA": "gene10_soda.fasta",
}


def read_fasta(path: Path) -> str:
    with open(path, encoding="utf-8") as f:
        seq = "".join(line.strip() for line in f if line.strip() and not line.startswith(">"))
    return clean_dna(seq)


def write_csv(path: Path, rows: list):
    if not rows:
        print(f"No rows to write: {path.name}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r.keys()))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {path}")


def main():
    parser = argparse.ArgumentParser(description="Run Pareto benchmark comparing NSGA-II and Weighted-GA")
    parser.add_argument("--smoke", action="store_true", help="Run quick smoke test on subset of genes")
    args = parser.parse_args()

    gene_items = list(GENES.items())
    if args.smoke:
        gene_items = gene_items[:2]
        pop_size = 20
        generations = 10
    else:
        pop_size = 60
        generations = 50

    print("=" * 70)
    print("PARETO MULTI-OBJECTIVE BENCHMARK: NSGA-II vs. WEIGHTED-GA")
    print(f"Genes: {len(gene_items)}, Mode: {'SMOKE' if args.smoke else 'FULL'}")
    print(f"Population: {pop_size}, Generations: {generations}")
    print("=" * 70)

    all_candidates = []
    summaries = []

    for gene_name, fname in gene_items:
        fasta_path = DATA_DIR / fname
        if not fasta_path.exists():
            print(f"Skipping {gene_name}: file {fasta_path} not found.")
            continue

        original = read_fasta(fasta_path)
        orig_cai = calculate_cai(original, "ecoli")
        orig_gc = gc_content(original)

        print(f"\nEvaluating Gene: {gene_name} ({len(original)} bp)")
        print(f"Original sequence: CAI = {orig_cai:.4f}, GC = {orig_gc:.4f}")

        # 1. NSGA-II
        start_nsga2 = time.perf_counter()
        nsga2_res = optimize_nsga2(
            sequence=original,
            organism="ecoli",
            target_gc=0.50,
            population_size=pop_size,
            generations=generations,
            seed=42,
        )
        elapsed_nsga2 = time.perf_counter() - start_nsga2

        for cand in nsga2_res:
            row = {
                "gene": gene_name,
                "method": "NSGA-II",
                "original_cai": round(orig_cai, 6),
                "cai": cand["cai"],
                "gc_deviation": cand["gc_deviation"],
                "gc_content": cand["gc_content"],
                "edit_rate": cand["edit_rate"],
                "feasible": cand["feasible"],
                "pareto_rank": cand["pareto_rank"],
                "crowding_distance": cand["crowding_distance"],
                "runtime_seconds": round(elapsed_nsga2, 4),
            }
            all_candidates.append(row)

        feasible_nsga2 = [c for c in nsga2_res if c["feasible"]]
        best_cai_nsga2 = max([c["cai"] for c in feasible_nsga2], default=None)
        print(
            f"  [NSGA-II] Found {len(nsga2_res)} non-dominated solutions "
            f"({len(feasible_nsga2)} feasible) in {elapsed_nsga2:.2f}s. "
            f"Best feasible CAI: {best_cai_nsga2}"
        )

        # 2. Weighted-GA Pareto Approximation
        start_wga = time.perf_counter()
        wga_res = optimize_pareto_weighted_ga(
            sequence=original,
            organism="ecoli",
            target_gc=0.50,
            population_size=pop_size,
            generations=generations,
            seed=42,
        )
        elapsed_wga = time.perf_counter() - start_wga

        for cand in wga_res:
            row = {
                "gene": gene_name,
                "method": "Weighted_GA",
                "original_cai": round(orig_cai, 6),
                "cai": cand["cai"],
                "gc_deviation": cand["gc_deviation"],
                "gc_content": cand["gc_content"],
                "edit_rate": cand["edit_rate"],
                "feasible": cand["feasible"],
                "pareto_rank": cand["pareto_rank"],
                "crowding_distance": "na",
                "runtime_seconds": round(elapsed_wga, 4),
            }
            all_candidates.append(row)

        feasible_wga = [c for c in wga_res if c["feasible"]]
        best_cai_wga = max([c["cai"] for c in feasible_wga], default=None)
        print(
            f"  [Weighted-GA] Found {len(wga_res)} non-dominated solutions "
            f"({len(feasible_wga)} feasible) in {elapsed_wga:.2f}s. "
            f"Best feasible CAI: {best_cai_wga}"
        )

        summaries.append({
            "gene": gene_name,
            "nsga2_candidates": len(nsga2_res),
            "nsga2_feasible_count": len(feasible_nsga2),
            "nsga2_best_cai": best_cai_nsga2,
            "nsga2_runtime_s": round(elapsed_nsga2, 4),
            "wga_candidates": len(wga_res),
            "wga_feasible_count": len(feasible_wga),
            "wga_best_cai": best_cai_wga,
            "wga_runtime_s": round(elapsed_wga, 4),
        })

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "pareto_candidates.csv", all_candidates)
    write_csv(OUT_DIR / "pareto_summary.csv", summaries)

    print("\nPareto benchmark complete.")


if __name__ == "__main__":
    main()
