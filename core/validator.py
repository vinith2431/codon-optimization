"""
Independent Sequence Validation Layer.

Every codon optimizer must pass through this independent validation layer.
An optimizer is never permitted to declare its candidate feasible merely
because its internal fitness function or objective score is high.
"""

from typing import Any, Dict, List, Optional
from core.sequence import clean_dna, validate_dna, gc_content
from core.similarity import edit_rate
from core.verification import translate_dna, verify_translation
from constraints.gc import gc_in_range
from constraints.local_gc import local_gc_content, local_gc_in_range
from constraints.homopolymer import find_homopolymers, homopolymer_constraint_satisfied
from constraints.repeats import find_repeats, repeat_constraint_satisfied
from constraints.restriction_sites import (
    find_restriction_sites,
    restriction_constraint_satisfied,
    COMMON_CLONING_ENZYMES,
)

# Benchmark standard constraint thresholds
DEFAULT_CONSTRAINTS = {
    "gc_min": 0.40,
    "gc_max": 0.60,
    "local_window": 60,
    "local_min": 0.30,
    "local_max": 0.70,
    "five_prime_window": 60,
    "five_prime_min": 0.40,
    "five_prime_max": 0.50,
    "max_homopolymer": 6,
    "repeat_length": 12,
    # Standard cloning enzymes avoided in recombinant E. coli expression vectors
    "forbidden_sites": [
        "GAATTC",  # EcoRI
        "GGATCC",  # BamHI
        "AAGCTT",  # HindIII
        "CATATG",  # NdeI
        "CTCGAG",  # XhoI
    ],
    "check_reverse_complement": True,
}


def longest_homopolymer_run(sequence: str) -> int:
    """Return the length of the longest identical-base homopolymer run in the sequence."""
    if not sequence:
        return 0
    longest = 1
    current = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current += 1
            if current > longest:
                longest = current
        else:
            current = 1
    return longest


def assess_candidate(
    original: str,
    candidate: str,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Perform an independent, exhaustive feasibility assessment of an optimized candidate.

    Parameters
    ----------
    original : str
        Original wild-type / baseline coding sequence.
    candidate : str
        Optimized coding sequence candidate.
    constraints : dict, optional
        Custom constraint settings. Defaults to DEFAULT_CONSTRAINTS.

    Returns
    -------
    dict
        Structured validation results with booleans, actual measurements, and failure reasons.
    """
    cfg = dict(DEFAULT_CONSTRAINTS)
    if constraints:
        cfg.update(constraints)

    original = clean_dna(original)
    candidate = clean_dna(candidate)

    failure_reasons: List[str] = []

    # 1. Sequence validity & Coding length
    valid_dna = False
    valid_length = False
    try:
        validate_dna(candidate, require_complete_codons=True)
        valid_dna = True
        valid_length = (len(candidate) % 3 == 0) and (len(candidate) == len(original))
    except Exception as exc:
        failure_reasons.append(f"Invalid DNA or coding length: {exc}")

    if len(candidate) != len(original):
        failure_reasons.append(
            f"Candidate length ({len(candidate)} bp) differs from original ({len(original)} bp)"
        )

    # 2. Protein preservation
    protein_preserved = False
    try:
        orig_protein = translate_dna(original)
        cand_protein = translate_dna(candidate)
        verif = verify_translation(original, candidate)
        verif_valid = verif.get("valid", False) if isinstance(verif, dict) else bool(verif)
        protein_preserved = (orig_protein == cand_protein) and verif_valid and valid_length
        if not protein_preserved:
            failure_reasons.append("Amino acid sequence was not preserved.")
    except Exception as exc:
        failure_reasons.append(f"Translation failed: {exc}")

    # 3. Global GC
    cand_gc = gc_content(candidate)
    gc_min = cfg.get("gc_min")
    gc_max = cfg.get("gc_max")
    global_gc_pass = gc_in_range(candidate, gc_min, gc_max)
    if not global_gc_pass:
        min_s = f"{gc_min:.2f}" if gc_min is not None else "-inf"
        max_s = f"{gc_max:.2f}" if gc_max is not None else "+inf"
        failure_reasons.append(
            f"Global GC {cand_gc:.4f} outside [{min_s}, {max_s}]"
        )

    # 4. 5' GC (First window)
    five_prime_window = cfg.get("five_prime_window", 60)
    five_min = cfg.get("five_prime_min")
    five_max = cfg.get("five_prime_max")
    five_prime_gc_val = None
    five_prime_gc_pass = True

    if len(candidate) >= five_prime_window:
        five_prime_gc_val = gc_content(candidate[:five_prime_window])
    else:
        five_prime_gc_val = cand_gc

    if five_min is not None and five_prime_gc_val < five_min:
        five_prime_gc_pass = False
    if five_max is not None and five_prime_gc_val > five_max:
        five_prime_gc_pass = False

    if not five_prime_gc_pass:
        min_s = f"{five_min:.2f}" if five_min is not None else "-inf"
        max_s = f"{five_max:.2f}" if five_max is not None else "+inf"
        failure_reasons.append(
            f"5' GC ({five_prime_gc_val:.4f} over first {min(len(candidate), five_prime_window)} bp) "
            f"outside [{min_s}, {max_s}]"
        )

    # 5. Local sliding-window GC
    local_window = cfg.get("local_window", 60)
    loc_min = cfg.get("local_min")
    loc_max = cfg.get("local_max")
    local_gc_pass = True
    min_local_gc = None
    max_local_gc = None

    if len(candidate) >= local_window:
        windows = local_gc_content(candidate, local_window)
        gc_vals = [w["gc_content"] for w in windows]
        min_local_gc = min(gc_vals) if gc_vals else cand_gc
        max_local_gc = max(gc_vals) if gc_vals else cand_gc
    else:
        min_local_gc = cand_gc
        max_local_gc = cand_gc

    if loc_min is not None and min_local_gc < loc_min:
        local_gc_pass = False
    if loc_max is not None and max_local_gc > loc_max:
        local_gc_pass = False

    if not local_gc_pass:
        min_s = f"{loc_min:.2f}" if loc_min is not None else "-inf"
        max_s = f"{loc_max:.2f}" if loc_max is not None else "+inf"
        failure_reasons.append(
            f"Local GC in [{min_local_gc:.4f}, {max_local_gc:.4f}] violates "
            f"[{min_s}, {max_s}] (window={local_window} bp)"
        )

    # 6. Homopolymers
    max_homopolymer = cfg.get("max_homopolymer")
    longest_homopolymer = longest_homopolymer_run(candidate)
    if max_homopolymer is not None:
        homopolymers = find_homopolymers(candidate, max_homopolymer)
        homopolymer_pass = len(homopolymers) == 0
        if not homopolymer_pass:
            failure_reasons.append(
                f"Longest homopolymer run is {longest_homopolymer} bp (exceeds max {max_homopolymer})"
            )
    else:
        homopolymers = []
        homopolymer_pass = True

    # 7. Repeated subsequences
    repeat_length = cfg.get("repeat_length")
    if repeat_length is not None:
        repeats = find_repeats(candidate, repeat_length)
        repeat_pass = len(repeats) == 0
        if not repeat_pass:
            failure_reasons.append(
                f"Found {len(repeats)} repeated substrings of length >= {repeat_length} bp"
            )
    else:
        repeats = []
        repeat_pass = True

    # 8. Restriction sites
    forbidden_sites = cfg.get("forbidden_sites")
    check_rc = cfg.get("check_reverse_complement", True)
    if forbidden_sites:
        restriction_violations = find_restriction_sites(
            candidate, forbidden_sites, check_reverse_complement=check_rc
        )
        restriction_pass = len(restriction_violations) == 0
        if not restriction_pass:
            found_names = [f"{v['site']}@{v['start']}" for v in restriction_violations]
            failure_reasons.append(
                f"Found forbidden restriction sites: {', '.join(found_names)}"
            )
    else:
        restriction_violations = []
        restriction_pass = True

    # 9. Conjunction of all constraints
    feasible = (
        protein_preserved
        and global_gc_pass
        and local_gc_pass
        and five_prime_gc_pass
        and homopolymer_pass
        and repeat_pass
        and restriction_pass
    )

    # Compute edit rate
    edit_rate_val = edit_rate(original, candidate)

    return {
        # Booleans
        "protein_preserved": bool(protein_preserved),
        "global_gc_pass": bool(global_gc_pass),
        "local_gc_pass": bool(local_gc_pass),
        "five_prime_gc_pass": bool(five_prime_gc_pass),
        "homopolymer_pass": bool(homopolymer_pass),
        "repeat_pass": bool(repeat_pass),
        "restriction_pass": bool(restriction_pass),
        "feasible": bool(feasible),
        # Quantitative measurements
        "global_gc": round(cand_gc, 6),
        "five_prime_gc": round(five_prime_gc_val, 6) if five_prime_gc_val is not None else None,
        "min_local_gc": round(min_local_gc, 6) if min_local_gc is not None else None,
        "max_local_gc": round(max_local_gc, 6) if max_local_gc is not None else None,
        "longest_homopolymer": int(longest_homopolymer),
        "repeat_count": len(repeats),
        "repeat_violations": repeats,
        "restriction_count": len(restriction_violations),
        "restriction_sites_found": restriction_violations,
        "edit_rate": round(edit_rate_val, 6),
        "failure_reasons": failure_reasons,
    }
