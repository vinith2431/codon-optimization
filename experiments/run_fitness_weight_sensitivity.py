"""
Fitness-Weight Sensitivity Benchmark.

Systematically evaluates how the relative weights of multi-objective penalties:
- Codon Adaptation Index (CAI)
- Codon-Pair Bias (CBP)
- K-mer frequency / sequence complexity
- 5' region GC constraint
affect optimized sequence feasibility, CAI, edit rate, and composition.

Outputs:
    results/csv/fitness_weight_sensitivity.csv  (Per-run records)
    results/csv/fitness_weight_summary.csv      (Aggregated statistics per weight configuration)
"""

from pathlib import Path
import argparse
import csv
import json
import statistics
import sys
import time
from typing import Any, Dict, List, Optional

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


def load_codon_pair_table() -> Optional[dict]:
    pair_path = ROOT / "data" / "codon_pair" / "ecoli_codon_pair_table.json"
    if not pair_path.exists():
        return None
    raw = json.loads(pair_path.read_text(encoding="utf-8"))
    table = {}
    for k, v in raw.items():
        parts = k.replace("-", "_").split("_")
        if len(parts) == 2:
            c1, c2 = parts[0].upper(), parts[1].upper()
            table[(c1, c2)] = float(v)
            table[f"{c1}-{c2}"] = float(v)
            table[f"{c1}_{c2}"] = float(v)
    return table or None


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


def get_weight_configs() -> List[Dict[str, Any]]:
    return [
        {
            "config_name": "baseline",
            "cai_weight": 1.0,
            "codon_pair_weight": 0.0,
            "kmer_weight": 0.0,
            "five_prime_gc_weight": 0.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "cbp_low",
            "cai_weight": 1.0,
            "codon_pair_weight": 1.0,
            "kmer_weight": 0.0,
            "five_prime_gc_weight": 0.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "cbp_high",
            "cai_weight": 1.0,
            "codon_pair_weight": 5.0,
            "kmer_weight": 0.0,
            "five_prime_gc_weight": 0.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "five_prime_low",
            "cai_weight": 1.0,
            "codon_pair_weight": 0.0,
            "kmer_weight": 0.0,
            "five_prime_gc_weight": 10.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "five_prime_high",
            "cai_weight": 1.0,
            "codon_pair_weight": 0.0,
            "kmer_weight": 0.0,
            "five_prime_gc_weight": 50.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "kmer_active",
            "cai_weight": 1.0,
            "codon_pair_weight": 0.0,
            "kmer_weight": 20.0,
            "five_prime_gc_weight": 0.0,
            "target_gc_weight": 10.0,
            "edit_weight": 0.0,
        },
        {
            "config_name": "combined_weights",
            "cai_weight": 1.0,
            "codon_pair_weight": 3.0,
            "kmer_weight": 10.0,
            "five_prime_gc_weight": 50.0,
            "target_gc_weight": 10.0,
            "edit_weight": 1.0,
        },
    ]


def summarize_weights(detailed_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
    parser = argparse.ArgumentParser(description="Run fitness-weight sensitivity experiment")
    parser.add_argument("--smoke", action="store_true", help="Run fast smoke test")
    args = parser.parse_args()

    gene_items = list(GENES.items())
    configs = get_weight_configs()
    pair_table = load_codon_pair_table()

    if args.smoke:
        gene_items = gene_items[:2]
        configs = [c for c in configs if c["config_name"] in ("baseline", "five_prime_high", "combined_weights")]
        seeds = [42]
        ga_pop, ga_gen = 25, 10
    else:
        seeds = SEEDS
        ga_pop, ga_gen = 60, 50

    print("=" * 75)
    print("FITNESS-WEIGHT SENSITIVITY BENCHMARK")
    print(f"Genes: {len(gene_items)}, Configs: {len(configs)}, Seeds: {seeds}")
    print(f"GA Population: {ga_pop}, Generations: {ga_gen}")
    print(f"Codon Pair Table loaded: {pair_table is not None}")
    print("=" * 75)

    detailed_rows: List[Dict[str, Any]] = []

    for gene_name, fname in gene_items:
        fasta_path = DATA_DIR / fname
        if not fasta_path.exists():
            continue
        original = read_fasta(fasta_path)
        orig_cai = calculate_cai(original, ORGANISM)

        for cfg in configs:
            for s in seeds:
                start = time.perf_counter()
                try:
                    cand = optimize_ga(
                        dna_sequence=original,
                        organism=ORGANISM,
                        population_size=ga_pop,
                        generations=ga_gen,
                        cai_weight=cfg["cai_weight"],
                        codon_pair_table=pair_table,
                        codon_pair_weight=cfg["codon_pair_weight"],
                        kmer_size=6,
                        max_kmer_frequency=0.0025,
                        kmer_weight=cfg["kmer_weight"],
                        five_prime_window=60,
                        five_prime_gc_min=0.40,
                        five_prime_gc_max=0.50,
                        five_prime_gc_weight=cfg["five_prime_gc_weight"],
                        target_gc=0.50,
                        target_gc_weight=cfg["target_gc_weight"],
                        edit_weight=cfg["edit_weight"],
                        gc_min=0.40,
                        gc_max=0.60,
                        local_gc_window=60,
                        local_gc_min=0.30,
                        local_gc_max=0.70,
                        repeat_length=12,
                        max_homopolymer=6,
                        forbidden_sites=DEFAULT_CONSTRAINTS["forbidden_sites"],
                        seed=s,
                        return_history=False,
                    )
                    elapsed = time.perf_counter() - start
                    ass = assess_candidate(original, cand, DEFAULT_CONSTRAINTS)
                    cand_cai = calculate_cai(cand, ORGANISM)
                    cand_gc = gc_content(cand)

                    detailed_rows.append({
                        "gene": gene_name,
                        "config_name": cfg["config_name"],
                        "seed": s,
                        "cai_weight": cfg["cai_weight"],
                        "codon_pair_weight": cfg["codon_pair_weight"],
                        "kmer_weight": cfg["kmer_weight"],
                        "five_prime_gc_weight": cfg["five_prime_gc_weight"],
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
                        "seed": s,
                        "cai_weight": cfg["cai_weight"],
                        "codon_pair_weight": cfg["codon_pair_weight"],
                        "kmer_weight": cfg["kmer_weight"],
                        "five_prime_gc_weight": cfg["five_prime_gc_weight"],
                        "original_cai": round(orig_cai, 6),
                        "optimized_cai": None,
                        "delta_cai": None,
                        "feasible": False,
                        "runtime_seconds": round(elapsed, 4),
                        "failure_reasons": f"Error: {exc}",
                    })

    summary_rows = summarize_weights(detailed_rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "fitness_weight_sensitivity.csv", detailed_rows)
    write_csv(OUT_DIR / "fitness_weight_summary.csv", summary_rows)

    print("\n" + "=" * 75)
    print("FITNESS-WEIGHT SENSITIVITY SUMMARY")
    print("=" * 75)
    print(f"{'Configuration':<24} | {'Runs':<5} | {'Feasible %':<10} | {'Mean CAI':<9} | {'Runtime (s)':<10}")
    print("-" * 75)
    for s in summary_rows:
        print(
            f"{s['config_name']:<24} | {s['total_runs']:<5} | {str(s['feasibility_rate_percent']) + '%':<10} | "
            f"{str(s['mean_feasible_cai']):<9} | {s['mean_runtime_seconds']:<10.3f}"
        )
    print("=" * 75)


if __name__ == "__main__":
    main()
