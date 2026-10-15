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
    global_gc_pass = gc_in_range(candidate, cfg["gc_min"], cfg["gc_max"])
    if not global_gc_pass:
        failure_reasons.append(
            f"Global GC {cand_gc:.4f} outside [{cfg['gc_min']:.2f}, {cfg['gc_max']:.2f}]"
        )

    # 4. 5' GC (First window)
    five_prime_window = cfg["five_prime_window"]
    five_prime_gc_val = None
    five_prime_gc_pass = False
    if len(candidate) >= five_prime_window:
        five_prime_gc_val = gc_content(candidate[:five_prime_window])
        five_prime_gc_pass = (cfg["five_prime_min"] <= five_prime_gc_val <= cfg["five_prime_max"])
        if not five_prime_gc_pass:
            failure_reasons.append(
                f"5' GC ({five_prime_gc_val:.4f} over first {five_prime_window} bp) "
                f"outside [{cfg['five_prime_min']:.2f}, {cfg['five_prime_max']:.2f}]"
            )
    else:
        # Sequence shorter than window
        five_prime_gc_val = cand_gc
        five_prime_gc_pass = (cfg["five_prime_min"] <= cand_gc <= cfg["five_prime_max"])
        if not five_prime_gc_pass:
            failure_reasons.append(
                f"5' GC ({cand_gc:.4f}) outside [{cfg['five_prime_min']:.2f}, {cfg['five_prime_max']:.2f}]"
            )

    # 5. Local sliding-window GC
    local_window = cfg["local_window"]
    local_gc_pass = False
    min_local_gc = None
    max_local_gc = None
    if len(candidate) >= local_window:
        windows = local_gc_content(candidate, local_window)
        gc_vals = [w["gc_content"] for w in windows]
        min_local_gc = min(gc_vals)
        max_local_gc = max(gc_vals)
        local_gc_pass = (min_local_gc >= cfg["local_min"]) and (max_local_gc <= cfg["local_max"])
        if not local_gc_pass:
            failure_reasons.append(
                f"Local GC in [{min_local_gc:.4f}, {max_local_gc:.4f}] violates "
                f"[{cfg['local_min']:.2f}, {cfg['local_max']:.2f}] (window={local_window} bp)"
            )
    else:
        min_local_gc = cand_gc
        max_local_gc = cand_gc
        local_gc_pass = (cfg["local_min"] <= cand_gc <= cfg["local_max"])

    # 6. Homopolymers
    homopolymers = find_homopolymers(candidate, cfg["max_homopolymer"])
    homopolymer_pass = len(homopolymers) == 0
    longest_homopolymer = longest_homopolymer_run(candidate)
    if not homopolymer_pass:
        failure_reasons.append(
            f"Longest homopolymer run is {longest_homopolymer} bp (exceeds max {cfg['max_homopolymer']})"
        )

    # 7. Repeated subsequences
    repeat_length = cfg["repeat_length"]
    repeats = find_repeats(candidate, repeat_length)
    repeat_pass = len(repeats) == 0
    if not repeat_pass:
        failure_reasons.append(
            f"Found {len(repeats)} repeated substrings of length >= {repeat_length} bp"
        )

    # 8. Restriction sites
    forbidden_sites = cfg.get("forbidden_sites", [])
    check_rc = cfg.get("check_reverse_complement", True)
    restriction_violations = find_restriction_sites(
        candidate, forbidden_sites, check_reverse_complement=check_rc
    )
    restriction_pass = len(restriction_violations) == 0
    if not restriction_pass:
        found_names = [f"{v['site']}@{v['start']}" for v in restriction_violations]
        failure_reasons.append(
            f"Found forbidden restriction sites: {', '.join(found_names)}"
        )

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
