"""
Tests for NSGA-II and Pareto Multi-Objective Optimization.
"""

import pytest
from core.verification import translate_dna, verify_translation
from optimizers.pareto import dominates, filter_non_dominated, optimize_pareto_weighted_ga
from optimizers.nsga2 import optimize_nsga2, fast_non_dominated_sort, Individual


def test_dominance_rule():
    # Candidate A: CAI=0.8, GC_dev=0.02, edit=0.1
    # Candidate B: CAI=0.7, GC_dev=0.05, edit=0.2 (A dominates B)
    cand_a = {"cai": 0.8, "gc_deviation": 0.02, "edit_rate": 0.1}
    cand_b = {"cai": 0.7, "gc_deviation": 0.05, "edit_rate": 0.2}
    assert dominates(cand_a, cand_b) is True
    assert dominates(cand_b, cand_a) is False

    # Candidate C: CAI=0.9, GC_dev=0.06, edit=0.15 (Trade-off: A has better GC & edit, C has better CAI)
    cand_c = {"cai": 0.9, "gc_deviation": 0.06, "edit_rate": 0.15}
    assert dominates(cand_a, cand_c) is False
    assert dominates(cand_c, cand_a) is False


def test_filter_non_dominated():
    candidates = [
        {"id": 1, "cai": 0.8, "gc_deviation": 0.02, "edit_rate": 0.1},
        {"id": 2, "cai": 0.7, "gc_deviation": 0.05, "edit_rate": 0.2},  # dominated by id 1
        {"id": 3, "cai": 0.9, "gc_deviation": 0.08, "edit_rate": 0.3},  # non-dominated
    ]
    front = filter_non_dominated(candidates)
    front_ids = [c["id"] for c in front]
    assert 1 in front_ids
    assert 3 in front_ids
    assert 2 not in front_ids


def test_nsga2_protein_preservation():
    seq = "ATGGCCGCGTTCAAACCCGGGTTTTAA"
    results = optimize_nsga2(
        sequence=seq,
        organism="ecoli",
        population_size=15,
        generations=5,
        seed=42,
    )
    assert len(results) > 0
    orig_protein = translate_dna(seq)
    for cand in results:
        cand_protein = translate_dna(cand["sequence"])
        assert orig_protein == cand_protein
        assert verify_translation(seq, cand["sequence"])["valid"] is True


def test_nsga2_non_dominated_front():
    seq = "ATGGCCGCGAAATTTGGT"
    results = optimize_nsga2(
        sequence=seq,
        organism="ecoli",
        population_size=20,
        generations=10,
        seed=42,
    )
    # Verify no solution in front_1 dominates another in front_1
    for i, a in enumerate(results):
        for j, b in enumerate(results):
            if i != j:
                assert dominates(a, b) is False, f"Candidate {i} dominates {j} in reported Pareto front!"


def test_nsga2_reproducibility():
    seq = "ATGGCCGCGAAATTTGGT"
    res1 = optimize_nsga2(sequence=seq, organism="ecoli", population_size=15, generations=5, seed=77)
    res2 = optimize_nsga2(sequence=seq, organism="ecoli", population_size=15, generations=5, seed=77)
    assert len(res1) == len(res2)
    for r1, r2 in zip(res1, res2):
        assert r1["sequence"] == r2["sequence"]
        assert r1["cai"] == r2["cai"]
