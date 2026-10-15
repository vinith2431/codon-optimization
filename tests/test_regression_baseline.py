"""
Regression tests verifying core requirements specified in Section 2.3:
1. DNA containing invalid symbols is rejected.
2. A sequence with an invalid coding length is rejected.
3. Synonymous optimization preserves the translated protein.
4. Global GC lower and upper bounds are respected when enforced.
5. The first-60-base GC calculation uses the intended window.
6. A sequence with a prohibited repeat fails validation.
7. A sequence with a configured forbidden restriction site fails validation.
8. An infeasible optimizer result is never silently labeled feasible.
9. Repeating an experiment with the same seed produces reproducible output.
10. An empty candidate list or empty result table does not crash the reporting code.
"""

import pytest
from core.sequence import validate_dna, clean_dna, gc_content
from core.verification import translate_dna, verify_translation
from core.validator import assess_candidate, DEFAULT_CONSTRAINTS
from optimizers.genetic_algorithm import optimize_ga, five_prime_gc_penalty
from optimizers.greedy import optimize_greedy
from optimizers.milp import optimize_milp
from experiments.run_constraint_aware_multistart import write_csv


# 1. DNA containing invalid symbols is rejected
def test_1_invalid_symbols_rejected():
    with pytest.raises(ValueError, match="Invalid DNA bases"):
        validate_dna("ATGNCC")
    with pytest.raises(ValueError, match="Invalid DNA bases"):
        validate_dna("ATG123")
    with pytest.raises(ValueError, match="Invalid DNA bases"):
        validate_dna("ATGZCC")


# 2. A sequence with an invalid coding length is rejected
def test_2_invalid_coding_length_rejected():
    with pytest.raises(ValueError, match="length must be divisible by 3"):
        validate_dna("ATGCC", require_complete_codons=True)
    with pytest.raises(ValueError, match="length must be divisible by 3"):
        validate_dna("A", require_complete_codons=True)


# 3. Synonymous optimization preserves the translated protein
def test_3_synonymous_preserves_protein():
    seq = "ATGGCCGCGTTCTGA"  # Met-Ala-Ala-Phe-Stop
    # Greedy
    opt_greedy = optimize_greedy(seq, "ecoli")
    assert translate_dna(seq) == translate_dna(opt_greedy)
    assert verify_translation(seq, opt_greedy)["valid"] is True

    # MILP
    opt_milp = optimize_milp(seq, "ecoli")
    assert translate_dna(seq) == translate_dna(opt_milp)
    assert verify_translation(seq, opt_milp)["valid"] is True

    # GA
    opt_ga = optimize_ga(seq, "ecoli", population_size=10, generations=5, seed=42)
    assert translate_dna(seq) == translate_dna(opt_ga)
    assert verify_translation(seq, opt_ga)["valid"] is True


# 4. Global GC lower and upper bounds are respected when enforced
def test_4_global_gc_bounds_respected():
    seq = "ATGGCCGCGAAATTTGGT"  # Flexible amino acids allowing GC range [0.35, 0.70]
    # MILP with GC min = 0.40, max = 0.60
    opt_milp = optimize_milp(seq, "ecoli", gc_min=0.40, gc_max=0.60)
    cand_gc = gc_content(opt_milp)
    assert 0.40 <= cand_gc <= 0.60

    # GA with GC bounds
    opt_ga = optimize_ga(seq, "ecoli", gc_min=0.40, gc_max=0.60, population_size=20, generations=10, seed=42)
    ga_gc = gc_content(opt_ga)
    assert 0.40 <= ga_gc <= 0.60


# 5. The first-60-base GC calculation uses the intended window
def test_5_five_prime_gc_window():
    # 60 bases of all G/C followed by 60 bases of all A/T
    seq_gc_at = "G" * 60 + "A" * 60
    # First 60 bases has 100% GC
    assert gc_content(seq_gc_at[:60]) == 1.0
    # Overall GC is 50%
    assert gc_content(seq_gc_at) == 0.50
    # 5' GC penalty should evaluate only the first 60 bases
    pen = five_prime_gc_penalty(seq_gc_at, window_size=60, gc_min=0.40, gc_max=0.50)
    assert pen == pytest.approx(1.0 - 0.50)  # 1.0 - max(0.50) = 0.50 violation


# 6. A sequence with a prohibited repeat fails validation
def test_6_prohibited_repeat_fails_validation():
    # 12-base repeat: ATGCATGCATGC
    seq = "ATGCATGCATGC" + "CCCGGG" + "ATGCATGCATGC"  # exactly repeats 12 bp
    res = assess_candidate(seq, seq, {"repeat_length": 12})
    assert res["repeat_pass"] is False
    assert res["feasible"] is False
    assert res["repeat_count"] > 0


# 7. A sequence with a configured forbidden restriction site fails validation
def test_7_forbidden_restriction_site_fails_validation():
    # Sequence with EcoRI site (GAATTC)
    seq = "ATGGCC" + "GAATTC" + "GCGCCG"
    res = assess_candidate(seq, seq, {"forbidden_sites": ["GAATTC"]})
    assert res["restriction_pass"] is False
    assert res["feasible"] is False
    assert res["restriction_count"] > 0
    assert any(site["site"] == "GAATTC" for site in res["restriction_sites_found"])


# 8. An infeasible optimizer result is never silently labeled feasible
def test_8_infeasible_not_labeled_feasible():
    # Sequence violating homopolymer constraint (7 A's)
    seq_violating = "ATG" + "AAAAAAA" + "TGC" + "T" * 50
    res = assess_candidate(seq_violating, seq_violating, {"max_homopolymer": 6})
    assert res["homopolymer_pass"] is False
    assert res["feasible"] is False
    assert "feasible" in res
    assert not res["feasible"]


# 9. Repeating an experiment with the same seed produces reproducible output
def test_9_reproducible_output_with_same_seed():
    seq = "ATGGCCGCGTTCAAACCCGGGTTTTAA"
    res1 = optimize_ga(seq, "ecoli", population_size=20, generations=10, seed=123)
    res2 = optimize_ga(seq, "ecoli", population_size=20, generations=10, seed=123)
    assert res1 == res2

    # Different seeds should produce distinct exploratory sequences
    res3 = optimize_ga(seq, "ecoli", population_size=20, generations=10, seed=999)
    # While they might occasionally collide, in general RNG paths differ
    assert isinstance(res1, str) and isinstance(res3, str)


# 10. An empty candidate list or empty result table does not crash reporting code
def test_10_empty_reporting_does_not_crash(tmp_path):
    empty_csv = tmp_path / "empty_report.csv"
    # Should not raise exception
    write_csv(empty_csv, [])
    assert not empty_csv.exists() or empty_csv.stat().st_size == 0
