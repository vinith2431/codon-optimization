"""
Constraint Sensitivity Benchmark.

Systematically tests the effect of constraint stringency on codon optimization:
1. Repeated-sequence thresholds: 8 bp, 10 bp, 12 bp
2. 5' GC window ranges: strict [0.42, 0.48], standard [0.40, 0.50], relaxed [0.35, 0.55]
3. Local sliding-window GC ranges: strict [0.35, 0.65], standard [0.30, 0.70], relaxed [0.25, 0.75]

Outputs:
    results/csv/constraint_sensitivity.csv  (Per-run records)
    results/csv/constraint_summary.csv      (Aggregated statistics per configuration)
"""

from pathlib import Path
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
from optimizers.genetic_algorithm import optimize_ga

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

SEEDS = [42, 123, 456]


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


def get_constraint_configs() -> List[Dict[str, Any]]:
    """Controlled experimental constraint configurations."""
    return [
        {
            "config_name": "baseline_standard",
            "category": "baseline",
            "repeat_length": 12,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        # 1. Repeat length sensitivity
        {
            "config_name": "repeat_strict_8bp",
            "category": "repeat_threshold",
            "repeat_length": 8,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        {
            "config_name": "repeat_moderate_10bp",
            "category": "repeat_threshold",
            "repeat_length": 10,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        {
            "config_name": "repeat_relaxed_12bp",
            "category": "repeat_threshold",
            "repeat_length": 12,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        # 2. 5' GC window sensitivity
        {
            "config_name": "five_prime_strict_42_48",
            "category": "five_prime_gc",
            "repeat_length": 12,
            "five_prime_min": 0.42,
            "five_prime_max": 0.48,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        {
            "config_name": "five_prime_relaxed_35_55",
            "category": "five_prime_gc",
            "repeat_length": 12,
            "five_prime_min": 0.35,
            "five_prime_max": 0.55,
            "local_min": 0.30,
            "local_max": 0.70,
        },
        # 3. Local GC window sensitivity
        {
            "config_name": "local_gc_strict_35_65",
            "category": "local_gc",
            "repeat_length": 12,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.35,
            "local_max": 0.65,
        },
        {
            "config_name": "local_gc_relaxed_25_75",
            "category": "local_gc",
            "repeat_length": 12,
            "five_prime_min": 0.40,
            "five_prime_max": 0.50,
            "local_min": 0.25,
            "local_max": 0.75,
        },
    ]


def summarize_constraints(detailed_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    configs = list(dict.fromkeys(r["config_name"] for r in detailed_rows))
    summary = []

    for cfg in configs:
        c_rows = [r for r in detailed_rows if r["config_name"] == cfg]
        total_runs = len(c_rows)
        feasible_runs = [r for r in c_rows if r.get("feasible") is True]

        feas_rate = len(feasible_runs) / total_runs if total_runs > 0 else 0.0
        feasible_cais = [r["optimized_cai"] for r in feasible_runs if r.get("optimized_cai") is not None]
        all_cais = [r["optimized_cai"] for r in c_rows if r.get("optimized_cai") is not None]
        runtimes = [r["runtime_seconds"] for r in c_rows if r.get("runtime_seconds") is not None]
        edits = [r["edit_rate"] for r in c_rows if r.get("edit_rate") is not None]
        gc_devs = [r["gc_deviation"] for r in c_rows if r.get("gc_deviation") is not None]

        summary.append({
            "config_name": cfg,
            "category": c_rows[0].get("category", ""),
            "total_runs": total_runs,
            "feasible_runs": len(feasible_runs),
            "feasibility_rate_percent": round(feas_rate * 100, 2),
            "mean_feasible_cai": round(statistics.mean(feasible_cais), 4) if feasible_cais else "None",
            "std_feasible_cai": (
                round(statistics.stdev(feasible_cais), 4) if len(feasible_cais) > 1 else 0.0
            ),
            "mean_overall_cai": round(statistics.mean(all_cais), 4) if all_cais else "None",
            "mean_runtime_seconds": round(statistics.mean(runtimes), 4) if runtimes else 0.0,
            "mean_edit_rate_percent": round(statistics.mean(edits) * 100, 2) if edits else 0.0,
            "mean_gc_deviation": round(statistics.mean(gc_devs), 4) if gc_devs else 0.0,
        })
    return summary


def main():
    parser = argparse.ArgumentParser(description="Run constraint sensitivity experiment")
    parser.add_argument("--smoke", action="store_true", help="Run quick smoke test")
    args = parser.parse_args()

    gene_items = list(GENES.items())
    configs = get_constraint_configs()

    if args.smoke:
        gene_items = gene_items[:2]
        configs = configs[:3]  # baseline, strict 8bp repeat, relaxed 12bp repeat
        seeds = [42]
        ga_pop, ga_gen = 25, 10
    else:
        seeds = SEEDS
        ga_pop, ga_gen = 60, 50

    print("=" * 75)
    print("CONSTRAINT SENSITIVITY BENCHMARK")
    print(f"Genes: {len(gene_items)}, Configs: {len(configs)}, Seeds: {seeds}")
    print(f"GA Population: {ga_pop}, Generations: {ga_gen}")
    print("=" * 75)

    detailed_rows: List[Dict[str, Any]] = []

    for gene_name, fname in gene_items:
        fasta_path = DATA_DIR / fname
        if not fasta_path.exists():
            continue
        original = read_fasta(fasta_path)
        orig_cai = calculate_cai(original, ORGANISM)

        for cfg in configs:
            test_constraints = dict(DEFAULT_CONSTRAINTS)
            test_constraints.update({
                "repeat_length": cfg["repeat_length"],
                "five_prime_min": cfg["five_prime_min"],
                "five_prime_max": cfg["five_prime_max"],
                "local_min": cfg["local_min"],
                "local_max": cfg["local_max"],
            })

            for s in seeds:
                start = time.perf_counter()
                try:
                    cand = optimize_ga(
                        dna_sequence=original,
                        organism=ORGANISM,
                        population_size=ga_pop,
                        generations=ga_gen,
                        gc_min=test_constraints["gc_min"],
                        gc_max=test_constraints["gc_max"],
                        local_gc_window=test_constraints["local_window"],
                        local_gc_min=test_constraints["local_min"],
                        local_gc_max=test_constraints["local_max"],
                        five_prime_window=test_constraints["five_prime_window"],
                        five_prime_gc_min=test_constraints["five_prime_min"],
                        five_prime_gc_max=test_constraints["five_prime_max"],
                        repeat_length=test_constraints["repeat_length"],
                        max_homopolymer=test_constraints["max_homopolymer"],
                        forbidden_sites=test_constraints["forbidden_sites"],
                        five_prime_gc_weight=50.0,
                        seed=s,
                        return_history=False,
                    )
                    elapsed = time.perf_counter() - start
                    ass = assess_candidate(original, cand, test_constraints)
                    cand_cai = calculate_cai(cand, ORGANISM)
                    cand_gc = gc_content(cand)

                    detailed_rows.append({
                        "gene": gene_name,
                        "config_name": cfg["config_name"],
                        "category": cfg["category"],
                        "seed": s,
                        "repeat_threshold": cfg["repeat_length"],
                        "five_prime_min": cfg["five_prime_min"],
                        "five_prime_max": cfg["five_prime_max"],
                        "local_min": cfg["local_min"],
                        "local_max": cfg["local_max"],
                        "original_cai": round(orig_cai, 6),
                        "optimized_cai": round(cand_cai, 6),
                        "delta_cai": round(cand_cai - orig_cai, 6),
                        "global_gc": ass["global_gc"],
                        "gc_deviation": round(abs(cand_gc - 0.50), 6),
                        "five_prime_gc": ass["five_prime_gc"],
                        "min_local_gc": ass["min_local_gc"],
                        "max_local_gc": ass["max_local_gc"],
                        "longest_homopolymer": ass["longest_homopolymer"],
                        "repeat_count": ass["repeat_count"],
                        "restriction_count": ass["restriction_count"],
                        "edit_rate": ass["edit_rate"],
                        "protein_preserved": ass["protein_preserved"],
                        "feasible": ass["feasible"],
                        "runtime_seconds": round(elapsed, 4),
                        "failure_reasons": "; ".join(ass["failure_reasons"]),
                    })
                    print(
                        f"  [{gene_name}] {cfg['config_name']} (seed={s}): "
                        f"CAI={cand_cai:.4f}, Feasible={ass['feasible']} in {elapsed:.2f}s"
                    )
                except Exception as exc:
                    elapsed = time.perf_counter() - start
                    detailed_rows.append({
                        "gene": gene_name,
                        "config_name": cfg["config_name"],
                        "category": cfg["category"],
                        "seed": s,
                        "repeat_threshold": cfg["repeat_length"],
                        "five_prime_min": cfg["five_prime_min"],
                        "five_prime_max": cfg["five_prime_max"],
                        "local_min": cfg["local_min"],
                        "local_max": cfg["local_max"],
                        "original_cai": round(orig_cai, 6),
                        "optimized_cai": None,
                        "delta_cai": None,
                        "feasible": False,
                        "runtime_seconds": round(elapsed, 4),
                        "failure_reasons": f"Error: {exc}",
                    })

    summary_rows = summarize_constraints(detailed_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "constraint_sensitivity.csv", detailed_rows)
    write_csv(OUT_DIR / "constraint_summary.csv", summary_rows)

    print("\n" + "=" * 75)
    print("CONSTRAINT SENSITIVITY SUMMARY")
    print("=" * 75)
    print(f"{'Configuration':<28} | {'Runs':<5} | {'Feasible %':<10} | {'Mean CAI':<9} | {'Runtime (s)':<10}")
    print("-" * 75)
    for s in summary_rows:
        print(
            f"{s['config_name']:<28} | {s['total_runs']:<5} | {str(s['feasibility_rate_percent']) + '%':<10} | "
            f"{str(s['mean_feasible_cai']):<9} | {s['mean_runtime_seconds']:<10.3f}"
        )
    print("=" * 75)


if __name__ == "__main__":
    main()