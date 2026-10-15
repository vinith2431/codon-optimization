"""
Pareto Multi-Objective Codon Optimization.

Provides:
- Non-dominated sorting and Pareto dominance filtering
- Weighted-GA Pareto front approximation
- High-level interface comparing NSGA-II and Weighted-GA
"""

from typing import Any, Dict, List, Optional
from core.sequence import clean_dna, gc_content
from core.similarity import edit_rate
from core.cai import calculate_cai
from core.validator import assess_candidate
from optimizers.genetic_algorithm import optimize_ga
from optimizers.nsga2 import optimize_nsga2


def dominates(
    cand_a: Dict[str, Any],
    cand_b: Dict[str, Any],
    cai_key: str = "cai",
    gc_dev_key: str = "gc_deviation",
    edit_key: str = "edit_rate",
) -> bool:
    """
    Check if candidate A dominates candidate B.
    Objectives:
    - Maximize CAI
    - Minimize GC deviation
    - Minimize edit rate

    Candidate A dominates B if A is at least as good on all objectives
    and strictly better on at least one.
    """
    a_cai = cand_a[cai_key]
    b_cai = cand_b[cai_key]
    a_gc = cand_a[gc_dev_key]
    b_gc = cand_b[gc_dev_key]
    a_edit = cand_a[edit_key]
    b_edit = cand_b[edit_key]

    not_worse = (a_cai >= b_cai) and (a_gc <= b_gc) and (a_edit <= b_edit)
    strictly_better = (a_cai > b_cai) or (a_gc < b_gc) or (a_edit < b_edit)

    return not_worse and strictly_better


def filter_non_dominated(
    candidates: List[Dict[str, Any]],
    cai_key: str = "cai",
    gc_dev_key: str = "gc_deviation",
    edit_key: str = "edit_rate",
) -> List[Dict[str, Any]]:
    """Filter a list of candidate dictionaries to only non-dominated candidates."""
    if not candidates:
        return []

    non_dominated = []
    for i, cand in enumerate(candidates):
        is_dominated = False
        for j, other in enumerate(candidates):
            if i != j and dominates(other, cand, cai_key, gc_dev_key, edit_key):
                is_dominated = True
                break
        if not is_dominated:
            non_dominated.append(cand)

    # Sort by CAI descending
    non_dominated.sort(key=lambda x: x[cai_key], reverse=True)
    return non_dominated


def optimize_pareto_weighted_ga(
    sequence: str,
    organism: str = "ecoli",
    target_gc: float = 0.50,
    constraints: Optional[Dict[str, Any]] = None,
    weights_configs: Optional[List[Dict[str, float]]] = None,
    population_size: int = 50,
    generations: int = 40,
    seed: int = 42,
) -> List[Dict[str, Any]]:
    """
    Weighted-GA Pareto approximation:
    Runs GA under several weight combinations trading off CAI, target GC, and edit rate,
    then filters the union of all results to find the non-dominated set.
    """
    clean_seq = clean_dna(sequence)

    # Default set of weight trade-offs
    if weights_configs is None:
        weights_configs = [
            # High CAI preference
            {"cai_weight": 2.0, "target_gc_weight": 5.0, "edit_weight": 0.0},
            # Balanced CAI and Target GC
            {"cai_weight": 1.0, "target_gc_weight": 15.0, "edit_weight": 0.0},
            # Low edit rate preservation
            {"cai_weight": 1.0, "target_gc_weight": 5.0, "edit_weight": 2.0},
            # Strict Target GC
            {"cai_weight": 0.5, "target_gc_weight": 25.0, "edit_weight": 0.0},
            # Balanced all-three
            {"cai_weight": 1.0, "target_gc_weight": 10.0, "edit_weight": 1.0},
        ]

    cfg = constraints or {}
    candidates = []

    for idx, w in enumerate(weights_configs):
        ga_kwargs = {
            "dna_sequence": clean_seq,
            "organism": organism,
            "target_gc": target_gc,
            "gc_min": cfg.get("gc_min", 0.40),
            "gc_max": cfg.get("gc_max", 0.60),
            "local_gc_window": cfg.get("local_window", 60),
            "local_gc_min": cfg.get("local_min", 0.30),
            "local_gc_max": cfg.get("local_max", 0.70),
            "five_prime_window": cfg.get("five_prime_window", 60),
            "five_prime_gc_min": cfg.get("five_prime_min", 0.40),
            "five_prime_gc_max": cfg.get("five_prime_max", 0.50),
            "repeat_length": cfg.get("repeat_length", 12),
            "max_homopolymer": cfg.get("max_homopolymer", 6),
            "forbidden_sites": cfg.get("forbidden_sites", []),
            "population_size": population_size,
            "generations": generations,
            "cai_weight": w.get("cai_weight", 1.0),
            "target_gc_weight": w.get("target_gc_weight", 10.0),
            "edit_weight": w.get("edit_weight", 0.0),
            "five_prime_gc_weight": 50.0,
            "seed": seed + idx,
            "return_history": False,
        }

        try:
            cand_seq = optimize_ga(**ga_kwargs)
            actual_gc = gc_content(cand_seq)
            cand_cai = calculate_cai(cand_seq, organism)
            cand_edit = edit_rate(clean_seq, cand_seq)
            cand_gc_dev = abs(actual_gc - target_gc)
            ass = assess_candidate(clean_seq, cand_seq, cfg)

            candidates.append({
                "sequence": cand_seq,
                "cai": round(cand_cai, 6),
                "gc_deviation": round(cand_gc_dev, 6),
                "gc_content": round(actual_gc, 6),
                "edit_rate": round(cand_edit, 6),
                "feasible": ass["feasible"],
                "pareto_rank": 1,
                "weight_config": w,
                "failure_reasons": ass.get("failure_reasons", []),
            })
        except Exception:
            continue

    # Filter non-dominated candidates
    return filter_non_dominated(candidates)


def optimize_pareto(
    sequence: str,
    organism: str = "ecoli",
    method: str = "nsga2",
    **kwargs,
) -> List[Dict[str, Any]]:
    """
    Main entry point for Pareto multi-objective codon optimization.

    Parameters
    ----------
    sequence : str
        Original coding DNA sequence.
    organism : str
        Target organism.
    method : str
        'nsga2' for true NSGA-II, or 'weighted_ga' for weighted-GA Pareto approximation.
    """
    if method == "nsga2":
        return optimize_nsga2(sequence=sequence, organism=organism, **kwargs)
    elif method == "weighted_ga":
        return optimize_pareto_weighted_ga(sequence=sequence, organism=organism, **kwargs)
    else:
        raise ValueError(f"Unknown Pareto method '{method}'; choose 'nsga2' or 'weighted_ga'.")
