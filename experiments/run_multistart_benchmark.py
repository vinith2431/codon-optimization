from pathlib import Path
import argparse
import csv
import statistics
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from core.sequence import gc_content
from core.cai import calculate_cai
from core.similarity import edit_rate
from core.verification import translate_dna, verify_translation
from optimizers.greedy import optimize_greedy
from optimizers.milp import optimize_milp
from optimizers.genetic_algorithm import optimize_ga

DATA_DIR = PROJECT_ROOT / "data" / "genes" / "fasta"
OUTPUT_DIR = PROJECT_ROOT / "results" / "csv"

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
ORGANISM = "ecoli"


def read_fasta(path):
    parts = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    sequence = "".join(parts).upper()
    if not sequence:
        raise ValueError(f"Empty FASTA file: {path}")
    if set(sequence) - set("ATGC") or len(sequence) % 3:
        raise ValueError(f"Invalid coding DNA sequence: {path}")
    return sequence


def protein_is_preserved(original, optimized):
    if len(original) != len(optimized):
        return False
    if set(optimized) - set("ATGC"):
        return False
    if len(optimized) % 3:
        return False

    try:
        same_protein = (
            translate_dna(original) == translate_dna(optimized)
        )
        verification = verify_translation(original, optimized)
        if isinstance(verification, dict):
            verified = verification.get("valid", False)
        else:
            verified = bool(verification)
        return bool(same_protein and verified)
    except Exception:
        return False


def measure(gene, method, original, optimized, runtime, starts, seed=None):
    original_cai = calculate_cai(original, ORGANISM)
    optimized_cai = calculate_cai(optimized, ORGANISM)
    original_gc = gc_content(original) * 100
    optimized_gc = gc_content(optimized) * 100
    rate = edit_rate(original, optimized) * 100
    protein_preserved = protein_is_preserved(original, optimized)

    return {
        "gene": gene,
        "method": method,
        "seed": seed if seed is not None else "",
        "starts": starts,
        "original_cai": round(original_cai, 6),
        "optimized_cai": round(optimized_cai, 6),
        "cai_change": round(optimized_cai - original_cai, 6),
        "original_gc_percent": round(original_gc, 4),
        "optimized_gc_percent": round(optimized_gc, 4),
        "gc_change_percent": round(optimized_gc - original_gc, 4),
        "edit_rate_percent": round(rate, 4),
        "protein_preserved": protein_preserved,
        "runtime_seconds": round(runtime, 4),
        "sequence_length_bp": len(optimized),
    }


def run_ga(original, seed, population, generations):
    result = optimize_ga(
        dna_sequence=original,
        organism=ORGANISM,
        population_size=population,
        generations=generations,
        mutation_rate=0.02,
        crossover_rate=0.80,
        tournament_size=3,
        elitism=2,
        cai_weight=1.0,
        gc_weight=10.0,
        local_gc_weight=5.0,
        repeat_weight=5.0,
        homopolymer_weight=5.0,
        restriction_weight=5.0,
        edit_weight=0.0,
        seed=seed,
        return_history=True,
    )
    if isinstance(result, dict):
        return result["sequence"]
    return result


def write_csv(path, rows):
    if not rows:
        print(f"No rows to write: {path.name}")
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare single-start GA, multi-start GA, Greedy and MILP."
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Quick check on two genes with a smaller GA configuration.",
    )
    args = parser.parse_args()

    genes = list(GENES.items())
    population = 30 if args.smoke else 100
    generations = 10 if args.smoke else 100
    if args.smoke:
        genes = genes[:2]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    errors = []

    print("=" * 72)
    print("SINGLE-START VS MULTI-START CODON OPTIMIZATION")
    print(f"Genes: {len(genes)} | Multi-start seeds: {SEEDS}")
    print(f"GA population={population}, generations={generations}")
    print("=" * 72)

    for gene, filename in genes:
        path = DATA_DIR / filename
        if not path.exists():
            errors.append({
                "gene": gene, "method": "input",
                "error": f"Missing FASTA: {path}"
            })
            print(f"Skipping {gene}: missing {path}")
            continue

        try:
            original = read_fasta(path)
            original_protein = translate_dna(original)
            print(f"\n{gene}: {len(original)} bp")

            # Single-start GA
            try:
                start = time.perf_counter()
                single = run_ga(original, 42, population, generations)
                elapsed = time.perf_counter() - start
                rows.append(measure(
                    gene, "GA_single", original, single,
                    elapsed, 1, seed=42
                ))
                print(f"  GA single: CAI={calculate_cai(single, ORGANISM):.4f}, "
                      f"time={elapsed:.2f}s")
            except Exception as exc:
                errors.append({"gene": gene, "method": "GA_single", "error": str(exc)})
                print(f"  GA single failed: {exc}")

            # Multi-start GA: select the candidate with highest CAI.
            # Total runtime includes all attempts.
            try:
                start = time.perf_counter()
                candidates = []
                for seed in SEEDS:
                    candidate = run_ga(original, seed, population, generations)
                    candidates.append((
                        calculate_cai(candidate, ORGANISM),
                        seed,
                        candidate
                    ))
                    print(f"    seed {seed}: CAI={candidates[-1][0]:.4f}")

                total_elapsed = time.perf_counter() - start
                best_cai, best_seed, best_sequence = max(
                    candidates, key=lambda item: item[0]
                )
                rows.append(measure(
                    gene, "GA_multistart", original, best_sequence,
                    total_elapsed, len(SEEDS), seed=best_seed
                ))
                print(f"  GA multi: best CAI={best_cai:.4f}, "
                      f"seed={best_seed}, total={total_elapsed:.2f}s")
            except Exception as exc:
                errors.append({"gene": gene, "method": "GA_multistart", "error": str(exc)})
                print(f"  GA multi-start failed: {exc}")

            # Greedy and MILP
            for method, optimizer in (
                ("Greedy", optimize_greedy),
                ("MILP", optimize_milp),
            ):
                try:
                    start = time.perf_counter()
                    optimized = optimizer(original, ORGANISM)
                    elapsed = time.perf_counter() - start
                    rows.append(measure(
                        gene, method, original, optimized, elapsed, 1
                    ))
                    print(f"  {method}: CAI={calculate_cai(optimized, ORGANISM):.4f}, "
                          f"time={elapsed:.2f}s")
                except Exception as exc:
                    errors.append({
                        "gene": gene, "method": method, "error": str(exc)
                    })
                    print(f"  {method} failed: {exc}")

        except Exception as exc:
            errors.append({"gene": gene, "method": "gene_setup", "error": str(exc)})
            print(f"Failed to process {gene}: {exc}")

    if rows:
        detail_path = OUTPUT_DIR / "multistart_benchmark.csv"
        write_csv(detail_path, rows)

        summary = []
        groups = sorted({
            (row["gene"], row["method"]) for row in rows
        })
        for gene, method in groups:
            group = [
                row for row in rows
                if row["gene"] == gene and row["method"] == method
            ]
            cai_values = [row["optimized_cai"] for row in group]
            summary.append({
                "gene": gene,
                "method": method,
                "runs_recorded": len(group),
                "mean_optimized_cai": round(statistics.mean(cai_values), 6),
                "best_optimized_cai": round(max(cai_values), 6),
                "mean_cai_change": round(
                    statistics.mean(row["cai_change"] for row in group), 6
                ),
                "mean_gc_percent": round(
                    statistics.mean(row["optimized_gc_percent"] for row in group), 4
                ),
                "mean_edit_rate_percent": round(
                    statistics.mean(row["edit_rate_percent"] for row in group), 4
                ),
                "all_proteins_preserved": all(
                    row["protein_preserved"] for row in group
                ),
                "mean_runtime_seconds": round(
                    statistics.mean(row["runtime_seconds"] for row in group), 4
                ),
            })
        write_csv(OUTPUT_DIR / "multistart_benchmark_summary.csv", summary)

    if errors:
        write_csv(OUTPUT_DIR / "multistart_benchmark_errors.csv", errors)

    print("\n" + "=" * 72)
    print(f"Successful result rows: {len(rows)}")
    print(f"Errors: {len(errors)}")
    print("Review protein_preserved before interpreting CAI results.")
    print("Greedy/MILP and GA use different objectives; compare metrics, not raw fitness.")
    print("=" * 72)


if __name__ == "__main__":
    main()
