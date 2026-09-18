"""
Global GC-content constraint utilities.
"""

from core.sequence import gc_content


def calculate_gc(sequence):
    """
    Return global GC content as a fraction.

    Example:
        0.55 means 55% GC.
    """

    return gc_content(sequence)


def gc_in_range(
    sequence,
    gc_min=None,
    gc_max=None
):
    """
    Check whether global GC content lies within
    the requested range.
    """

    gc = calculate_gc(sequence)

    if gc_min is not None and gc < gc_min:
        return False

    if gc_max is not None and gc > gc_max:
        return False

    return True


def gc_constraint_violation(
    sequence,
    gc_min=None,
    gc_max=None
):
    """
    Return the amount by which the sequence violates
    the requested GC range.

    Returns 0 when there is no violation.
    """

    gc = calculate_gc(sequence)

    violation = 0.0

    if gc_min is not None and gc < gc_min:
        violation += gc_min - gc

    if gc_max is not None and gc > gc_max:
        violation += gc - gc_max

    return violation


def gc_report(
    sequence,
    gc_min=None,
    gc_max=None
):
    """
    Generate a simple GC-content report.
    """

    gc = calculate_gc(sequence)

    return {
        "gc_content": gc,
        "gc_percentage": gc * 100,
        "minimum": gc_min,
        "maximum": gc_max,
        "valid": gc_in_range(
            sequence,
            gc_min,
            gc_max
        ),
        "violation": gc_constraint_violation(
            sequence,
            gc_min,
            gc_max
        )
    }