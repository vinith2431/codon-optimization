from pathlib import Path
import argparse
import csv
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.cai import calculate_cai
from core.sequence import gc_content
from core.similarity import edit_rate
from core.verification import translate_dna, verify_translation
from core.validator import assess_candidate, DEFAULT_CONSTRAINTS
from optimizers.genetic_algorithm import optimize_ga
from optimizers.greedy import optimize_greedy
from optimizers.milp import optimize_milp

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

# Explicit experiment settings; fractions are in [0, 1].
GC_MIN, GC_MAX = 0.40, 0.60
LOCAL_WINDOW, LOCAL_MIN, LOCAL_MAX = 60, 0.30, 0.70
FIVE_PRIME_WINDOW, FIVE_PRIME_MIN, FIVE_PRIME_MAX = 60, 0.40, 0.50
MAX_HOMOPOLYMER = 6
REPEAT_LENGTH = 12

# Standard biotechnology restriction sites avoided in E. coli expression vectors
# (EcoRI, BamHI, HindIII, NdeI, XhoI)
FORBIDDEN_SITES = [
    "GAATTC",  # EcoRI
    "GGATCC",  # BamHI
    "AAGCTT",  # HindIII
    "CATATG",  # NdeI
    "CTCGAG",  # XhoI
]

CONSTRAINTS_CONFIG = {
    "gc_min": GC_MIN,
    "gc_max": GC_MAX,
    "local_window": LOCAL_WINDOW,
    "local_min": LOCAL_MIN,
    "local_max": LOCAL_MAX,
    "five_prime_window": FIVE_PRIME_WINDOW,
    "five_prime_min": FIVE_PRIME_MIN,
    "five_prime_max": FIVE_PRIME_MAX,
    "max_homopolymer": MAX_HOMOPOLYMER,
    "repeat_length": REPEAT_LENGTH,
    "forbidden_sites": FORBIDDEN_SITES,
    "check_reverse_complement": True,
}


def read_fasta(path):
    parts = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    sequence = "".join(parts).upper()
    if not sequence or set(sequence) - set("ATGC") or len(sequence) % 3:
        raise ValueError(f"Invalid coding DNA sequence: {path}")
    return sequence


def assess(original, candidate):
    """Assess candidate using centralized independent validation layer."""
    res = assess_candidate(original, candidate, CONSTRAINTS_CONFIG)
    # Format non-scalar fields for clean CSV export
    res_formatted = dict(res)
    res_formatted["global_gc_percent"] = round(res["global_gc"] * 100, 4) if res["global_gc"] is not None else None
    res_formatted["five_prime_gc_percent"] = round(res["five_prime_gc"] * 100, 4) if res["five_prime_gc"] is not None else None
    res_formatted["failure_reasons_str"] = "; ".join(res["failure_reasons"])
    res_formatted["restriction_sites_str"] = "; ".join(
        [f"{v['site']}@{v['start']}" for v in res["restriction_sites_found"]]
    )
    # Exclude complex nested lists from direct CSV columns to prevent formatting issues
    res_formatted.pop("repeat_violations", None)
    res_formatted.pop("restriction_sites_found", None)
    res_formatted.pop("failure_reasons", None)
    return res_formatted


def run_ga(original, seed, population, generations):
    result = optimize_ga(
        dna_sequence=original,
        organism="ecoli",
        population_size=population,
        generations=generations,
        mutation_rate=0.02,
        crossover_rate=0.80,
        tournament_size=3,
        elitism=2,
        cai_weight=1.0,
        gc_min=GC_MIN,
        gc_max=GC_MAX,
        local_gc_window=LOCAL_WINDOW,
        local_gc_min=LOCAL_MIN,
        local_gc_max=LOCAL_MAX,
        repeat_length=REPEAT_LENGTH,
        max_homopolymer=MAX_HOMOPOLYMER,
        forbidden_sites=FORBIDDEN_SITES,
        five_prime_window=FIVE_PRIME_WINDOW,
        five_prime_gc_min=FIVE_PRIME_MIN,
        five_prime_gc_max=FIVE_PRIME_MAX,
        gc_weight=10.0,
        local_gc_weight=5.0,
        repeat_weight=5.0,
        homopolymer_weight=5.0,
        restriction_weight=10.0,
        five_prime_gc_weight=50.0,
        edit_weight=0.0,
        seed=seed,
        return_history=True,
    )
    return result["sequence"] if isinstance(result, dict) else result


def write_csv(path, rows):
    if not rows:
        print(f"No rows to write: {path.name}")
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row.keys())))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {path}")


def result_row(gene, method, original, candidate, elapsed, seed="", error=""):
    if error or candidate is None:
        return {
            "gene": gene,
            "method": method,
            "seed": seed,
            "original_cai": round(calculate_cai(original, "ecoli"), 6),
            "optimized_cai": None,
            "cai_change": None,
            "edit_rate_percent": None,
            "runtime_seconds": round(elapsed, 4),
            "feasible": False,
            "error": error,
        }

    checks = assess(original, candidate)
    row = {
        "gene": gene,
        "method": method,
        "seed": seed,
        "original_cai": round(calculate_cai(original, "ecoli"), 6),
        "optimized_cai": round(calculate_cai(candidate, "ecoli"), 6),
        "cai_change": round(
            calculate_cai(candidate, "ecoli") - calculate_cai(original, "ecoli"), 6
        ),
        "edit_rate_percent": round(edit_rate(original, candidate) * 100, 4),
        "runtime_seconds": round(elapsed, 4),
        "error": "",
    }
    row.update(checks)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    genes = list(GENES.items())
    population, generations = (30, 10) if args.smoke else (100, 100)
    if args.smoke:
        genes = genes[:2]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    detailed, errors, summaries = [], [], []

    print("Constraint-aware multi-start benchmark")
    print(f"Genes={len(genes)}; seeds={SEEDS}")
    print(
        f"GC={GC_MIN:.2f}-{GC_MAX:.2f}; local GC={LOCAL_MIN:.2f}-{LOCAL_MAX:.2f}"
        f" / {LOCAL_WINDOW} bp; 5-prime GC={FIVE_PRIME_MIN:.2f}-{FIVE_PRIME_MAX:.2f}"
    )
    print(f"Forbidden restriction sites: {FORBIDDEN_SITES}")

    for gene, filename in genes:
        path = DATA_DIR / filename
        if not path.exists():
            errors.append({"gene": gene, "method": "input", "error": f"Missing {path}"})
            continue

        try:
            original = read_fasta(path)
            candidates = []
            print(f"\n{gene}: {len(original)} bp")

            # Independent GA runs
            for seed in SEEDS:
                try:
                    start = time.perf_counter()
                    candidate = run_ga(original, seed, population, generations)
                    elapsed = time.perf_counter() - start
                    row = result_row(
                        gene, "GA_seed", original, candidate, elapsed, seed
                    )
                    detailed.append(row)
                    candidates.append((row, candidate))
                    print(
                        f"  seed={seed}: CAI={row['optimized_cai']:.4f}, "
                        f"feasible={row['feasible']}"
                    )
                except Exception as exc:
                    elapsed = time.perf_counter() - start if 'start' in locals() else 0.0
                    err_row = result_row(gene, "GA_seed", original, None, elapsed, seed, error=str(exc))
                    detailed.append(err_row)
                    errors.append({"gene": gene, "method": f"GA_seed_{seed}", "error": str(exc)})

            feasible_candidates = [
                item for item in candidates if item[0]["feasible"]
            ]
            if feasible_candidates:
                best_row, best_sequence = max(
                    feasible_candidates, key=lambda item: item[0]["optimized_cai"]
                )
                selection_status = "best_feasible"
            elif candidates:
                best_row, best_sequence = min(
                    candidates,
                    key=lambda item: sum(
                        not item[0][key] for key in (
                            "protein_preserved", "global_gc_pass", "local_gc_pass",
                            "five_prime_gc_pass", "homopolymer_pass", "repeat_pass",
                            "restriction_pass"
                        ) if key in item[0]
                    ),
                )
                selection_status = "no_feasible_candidate_least_violations"
            else:
                best_row, best_sequence = None, None
                selection_status = "all_ga_runs_failed"

            if best_row is not None:
                selected = dict(best_row)
                selected["method"] = "GA_selected"
                selected["seed"] = best_row["seed"]
                selected["selection_status"] = selection_status
                detailed.append(selected)

            # Greedy and MILP receive exactly the same independent feasibility checks.
            for method, optimizer in (("Greedy", optimize_greedy), ("MILP", optimize_milp)):
                try:
                    start = time.perf_counter()
                    candidate = optimizer(original, "ecoli")
                    elapsed = time.perf_counter() - start
                    row = result_row(gene, method, original, candidate, elapsed)
                    detailed.append(row)
                    print(
                        f"  {method}: CAI={row['optimized_cai']:.4f}, "
                        f"feasible={row['feasible']}"
                    )
                except Exception as exc:
                    elapsed = time.perf_counter() - start if 'start' in locals() else 0.0
                    err_row = result_row(gene, method, original, None, elapsed, error=str(exc))
                    detailed.append(err_row)
                    errors.append({"gene": gene, "method": method, "error": str(exc)})

            summaries.append({
                "gene": gene,
                "ga_seeds_completed": len(candidates),
                "ga_feasible_candidates": len(feasible_candidates),
                "ga_selection_status": selection_status,
                "ga_selected_cai": best_row["optimized_cai"] if best_row else None,
                "ga_selected_feasible": best_row["feasible"] if best_row else False,
            })

        except Exception as exc:
            errors.append({"gene": gene, "method": "gene_setup", "error": str(exc)})
            print(f"ERROR {gene}: {exc}")

    write_csv(OUT_DIR / "constraint_aware_multistart.csv", detailed)
    write_csv(OUT_DIR / "constraint_aware_multistart_summary.csv", summaries)
    if errors:
        write_csv(OUT_DIR / "constraint_aware_multistart_errors.csv", errors)

    print(f"\nDetailed rows: {len(detailed)}")
    print(f"Errors: {len(errors)}")
    print("Review feasibility columns before interpreting CAI.")


if __name__ == "__main__":
    main()
