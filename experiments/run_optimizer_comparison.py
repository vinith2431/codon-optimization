"""
Comprehensive and Fair Optimizer Comparison Benchmark.

Algorithms compared under identical biological benchmarks and independent validation:
1. Greedy baseline
2. MILP (Mathematical optimization with global GC bounds)
3. MILP_constrained (MILP with global GC, 5' GC, and local GC bounds)
4. GA_single (Evolutionary single-seed baseline, seed 42)
5. GA_multistart (Independent GA runs across 5 distinct random seeds, each recorded as a separate row)
6. NSGA-II (True multi-objective Pareto optimization)

Outputs:
    results/csv/optimizer_comparison.csv  (Per-run records)
    results/csv/optimizer_summary.csv     (Aggregated metrics: feasibility rate, feasible CAI, runtime)
"""

from pathlib import Path
from typing import Optional
import argparse
import csv
import statistics
import sys
import time
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.cai import calculate_cai
from core.sequence import clean_dna, gc_content
from core.similarity import edit_rate
from core.validator import assess_candidate, DEFAULT_CONSTRAINTS
from optimizers.greedy import optimize_greedy
from optimizers.milp import optimize_milp
from optimizers.genetic_algorithm import optimize_ga
from optimizers.nsga2 import optimize_nsga2

ORGANISM = "ecoli"
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

SEEDS = [42, 123, 456, 789, 2026]


def read_fasta(path: Path) -> str:
    with open(path, encoding="utf-8") as f:
        seq = "".join(line.strip() for line in f if line.strip() and not line.startswith(">"))
    return clean_dna(seq)


def write_csv(path: Path, rows: List[Dict[str, Any]]):
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


def build_result_row(
    gene: str,
    method: str,
    original: str,
    candidate: Optional[str],
    elapsed: float,
    seed: Any = "na",
    error: str = "",
) -> Dict[str, Any]:
    orig_cai = calculate_cai(original, ORGANISM)
    orig_gc = gc_content(original)
    orig_p5 = gc_content(original[:60]) if len(original) >= 60 else orig_gc

    base_row = {
        "gene": gene,
        "method": method,
        "seed": seed,
        "original_length_bp": len(original),
        "original_cai": round(orig_cai, 6),
        "original_gc": round(orig_gc, 6),
        "original_five_prime_gc": round(orig_p5, 6),
        "runtime_seconds": round(elapsed, 4),
        "error": error,
    }

    if error or candidate is None:
        base_row.update({
            "optimized_length_bp": None,
            "optimized_cai": None,
            "delta_cai": None,
            "global_gc": None,
            "five_prime_gc": None,
            "min_local_gc": None,
            "max_local_gc": None,
            "longest_homopolymer": None,
            "repeat_count": None,
            "restriction_count": None,
            "edit_rate": None,
            "protein_preserved": False,
            "global_gc_pass": False,
            "local_gc_pass": False,
            "five_prime_gc_pass": False,
            "homopolymer_pass": False,
            "repeat_pass": False,
            "restriction_pass": False,
            "feasible": False,
            "failure_reasons": error,
        })
        return base_row

    ass = assess_candidate(original, candidate, DEFAULT_CONSTRAINTS)
    cand_cai = calculate_cai(candidate, ORGANISM)
    delta_cai = cand_cai - orig_cai

    base_row.update({
        "optimized_length_bp": len(candidate),
        "optimized_cai": round(cand_cai, 6),
        "delta_cai": round(delta_cai, 6),
        "global_gc": ass["global_gc"],
        "five_prime_gc": ass["five_prime_gc"],
        "min_local_gc": ass["min_local_gc"],
        "max_local_gc": ass["max_local_gc"],
        "longest_homopolymer": ass["longest_homopolymer"],
        "repeat_count": ass["repeat_count"],
        "restriction_count": ass["restriction_count"],
        "edit_rate": ass["edit_rate"],
        "protein_preserved": ass["protein_preserved"],
        "global_gc_pass": ass["global_gc_pass"],
        "local_gc_pass": ass["local_gc_pass"],
        "five_prime_gc_pass": ass["five_prime_gc_pass"],
        "homopolymer_pass": ass["homopolymer_pass"],
        "repeat_pass": ass["repeat_pass"],
        "restriction_pass": ass["restriction_pass"],
        "feasible": ass["feasible"],
        "failure_reasons": "; ".join(ass["failure_reasons"]),
    })
    return base_row


def generate_summary(detailed_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    methods = list(dict.fromkeys(r["method"] for r in detailed_rows))
    summary = []

    for method in methods:
        m_rows = [r for r in detailed_rows if r["method"] == method]
        total_runs = len(m_rows)
        valid_runs = [r for r in m_rows if r.get("optimized_cai") is not None]
        feasible_runs = [r for r in m_rows if r.get("feasible") is True]

        feasible_rate = (len(feasible_runs) / total_runs) if total_runs > 0 else 0.0

        feasible_cais = [r["optimized_cai"] for r in feasible_runs]
        all_cais = [r["optimized_cai"] for r in valid_runs]
        runtimes = [r["runtime_seconds"] for r in m_rows]
        edits = [r["edit_rate"] for r in valid_runs if r.get("edit_rate") is not None]

        # Top failure reasons
        all_reasons = []
        for r in m_rows:
            reasons = r.get("failure_reasons", "")
            if reasons:
                for reason in reasons.split("; "):
                    if reason.strip():
                        # Truncate detail for clean summary grouping
                        short = reason.split("(")[0].split("@")[0].strip()
                        all_reasons.append(short)

        from collections import Counter
        top_failures = [f"{k}({v})" for k, v in Counter(all_reasons).most_common(3)]

        summary.append({
            "method": method,
            "total_runs": total_runs,
            "feasible_runs": len(feasible_runs),
            "feasibility_rate_percent": round(feasible_rate * 100, 2),
            "best_feasible_cai": round(max(feasible_cais), 4) if feasible_cais else "None",
            "mean_feasible_cai": round(statistics.mean(feasible_cais), 4) if feasible_cais else "None",
            "std_feasible_cai": (
                round(statistics.stdev(feasible_cais), 4) if len(feasible_cais) > 1 else 0.0
            ),
            "mean_overall_cai": round(statistics.mean(all_cais), 4) if all_cais else "None",
            "mean_runtime_seconds": round(statistics.mean(runtimes), 4) if runtimes else 0.0,
            "mean_edit_rate_percent": round(statistics.mean(edits) * 100, 2) if edits else 0.0,
            "top_failure_reasons": "; ".join(top_failures) if top_failures else "None",
        })

    return summary


def main():
    parser = argparse.ArgumentParser(description="Run fair optimizer comparison")
    parser.add_argument("--smoke", action="store_true", help="Run fast smoke test")
    args = parser.parse_args()

    gene_items = list(GENES.items())
    if args.smoke:
        gene_items = gene_items[:2]
        ga_pop, ga_gen = 25, 10
        seeds = SEEDS[:2]
    else:
        ga_pop, ga_gen = 100, 100
        seeds = SEEDS

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    detailed_rows: List[Dict[str, Any]] = []

    print("=" * 75)
    print("CODON OPTIMIZER BENCHMARK: FAIR INDEPENDENT COMPARISON")
    print(f"Genes: {len(gene_items)}, Mode: {'SMOKE' if args.smoke else 'FULL'}")
    print(f"GA Population: {ga_pop}, Generations: {ga_gen}, Seeds: {seeds}")
    print("=" * 75)

    for gene_name, fname in gene_items:
        fasta_path = DATA_DIR / fname
        if not fasta_path.exists():
            print(f"Skipping {gene_name}: file not found.")
            continue

        original = read_fasta(fasta_path)
        print(f"\nProcessing Gene: {gene_name} ({len(original)} bp)")

        # 1. Greedy
        start = time.perf_counter()
        try:
            cand = optimize_greedy(original, ORGANISM)
            row = build_result_row(gene_name, "Greedy", original, cand, time.perf_counter() - start)
        except Exception as exc:
            row = build_result_row(gene_name, "Greedy", original, None, time.perf_counter() - start, error=str(exc))
        detailed_rows.append(row)
        print(f"  Greedy: CAI={row['optimized_cai']}, Feasible={row['feasible']}")

        # 2. MILP (standard with global GC)
        start = time.perf_counter()
        try:
            cand = optimize_milp(original, ORGANISM, gc_min=0.40, gc_max=0.60)
            row = build_result_row(gene_name, "MILP", original, cand, time.perf_counter() - start)
        except Exception as exc:
            row = build_result_row(gene_name, "MILP", original, None, time.perf_counter() - start, error=str(exc))
        detailed_rows.append(row)
        print(f"  MILP: CAI={row['optimized_cai']}, Feasible={row['feasible']}")

        # 3. MILP_constrained (with 5' GC, local GC, and restriction site cuts)
        start = time.perf_counter()
        try:
            cand = optimize_milp(
                original,
                ORGANISM,
                gc_min=0.40,
                gc_max=0.60,
                five_prime_gc_min=0.40,
                five_prime_gc_max=0.50,
                five_prime_window=60,
                local_gc_min=0.30,
                local_gc_max=0.70,
                local_gc_window=60,
                forbidden_sites=DEFAULT_CONSTRAINTS["forbidden_sites"],
                max_homopolymer=6,
            )
            row = build_result_row(gene_name, "MILP_constrained", original, cand, time.perf_counter() - start)
        except Exception as exc:
            row = build_result_row(gene_name, "MILP_constrained", original, None, time.perf_counter() - start, error=str(exc))
        detailed_rows.append(row)
        print(f"  MILP_constrained: CAI={row['optimized_cai']}, Feasible={row['feasible']}")

        # 4. GA_single (seed 42)
        start = time.perf_counter()
        try:
            cand = optimize_ga(
                dna_sequence=original,
                organism=ORGANISM,
                population_size=ga_pop,
                generations=ga_gen,
                gc_min=0.40,
                gc_max=0.60,
                five_prime_window=60,
                five_prime_gc_min=0.40,
                five_prime_gc_max=0.50,
                local_gc_window=60,
                local_gc_min=0.30,
                local_gc_max=0.70,
                repeat_length=12,
                max_homopolymer=6,
                forbidden_sites=DEFAULT_CONSTRAINTS["forbidden_sites"],
                five_prime_gc_weight=50.0,
                seed=42,
                return_history=False,
            )
            row = build_result_row(gene_name, "GA_single", original, cand, time.perf_counter() - start, seed=42)
        except Exception as exc:
            row = build_result_row(gene_name, "GA_single", original, None, time.perf_counter() - start, seed=42, error=str(exc))
        detailed_rows.append(row)
        print(f"  GA_single: CAI={row['optimized_cai']}, Feasible={row['feasible']}")

        # 5. GA_multistart (retain all seeds as separate rows!)
        for s in seeds:
            start = time.perf_counter()
            try:
                cand = optimize_ga(
                    dna_sequence=original,
                    organism=ORGANISM,
                    population_size=ga_pop,
                    generations=ga_gen,
                    gc_min=0.40,
                    gc_max=0.60,
                    five_prime_window=60,
                    five_prime_gc_min=0.40,
                    five_prime_gc_max=0.50,
                    local_gc_window=60,
                    local_gc_min=0.30,
                    local_gc_max=0.70,
                    repeat_length=12,
                    max_homopolymer=6,
                    forbidden_sites=DEFAULT_CONSTRAINTS["forbidden_sites"],
                    five_prime_gc_weight=50.0,
                    seed=s,
                    return_history=False,
                )
                row = build_result_row(gene_name, "GA_multistart", original, cand, time.perf_counter() - start, seed=s)
            except Exception as exc:
                row = build_result_row(gene_name, "GA_multistart", original, None, time.perf_counter() - start, seed=s, error=str(exc))
            detailed_rows.append(row)
            print(f"  GA_multistart (seed={s}): CAI={row['optimized_cai']}, Feasible={row['feasible']}")

        # 6. NSGA-II (select best feasible candidate)
        start = time.perf_counter()
        try:
            nsga_front = optimize_nsga2(
                sequence=original,
                organism=ORGANISM,
                target_gc=0.50,
                population_size=ga_pop,
                generations=ga_gen,
                seed=42,
            )
            feasible_cands = [c for c in nsga_front if c["feasible"]]
            if feasible_cands:
                best_cand = max(feasible_cands, key=lambda c: c["cai"])["sequence"]
            elif nsga_front:
                best_cand = nsga_front[0]["sequence"]
            else:
                best_cand = None
            row = build_result_row(gene_name, "NSGA-II", original, best_cand, time.perf_counter() - start, seed=42)
        except Exception as exc:
            row = build_result_row(gene_name, "NSGA-II", original, None, time.perf_counter() - start, seed=42, error=str(exc))
        detailed_rows.append(row)
        print(f"  NSGA-II: CAI={row['optimized_cai']}, Feasible={row['feasible']}")

    # Write per-run detailed results and aggregated summary
    summary_rows = generate_summary(detailed_rows)

    write_csv(OUT_DIR / "optimizer_comparison.csv", detailed_rows)
    write_csv(OUT_DIR / "optimizer_summary.csv", summary_rows)

    print("\n" + "=" * 75)
    print("OPTIMIZER SUMMARY BENCHMARK")
    print("=" * 75)
    print(f"{'Method':<18} | {'Runs':<5} | {'Feasible %':<10} | {'Best CAI':<9} | {'Mean CAI':<9} | {'Runtime (s)':<10}")
    print("-" * 75)
    for s in summary_rows:
        print(
            f"{s['method']:<18} | {s['total_runs']:<5} | {str(s['feasibility_rate_percent']) + '%':<10} | "
            f"{str(s['best_feasible_cai']):<9} | {str(s['mean_feasible_cai']):<9} | {s['mean_runtime_seconds']:<10.3f}"
        )
    print("=" * 75)


if __name__ == "__main__":
    main()