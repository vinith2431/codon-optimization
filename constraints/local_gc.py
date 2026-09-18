"""
Local/sliding-window GC-content constraints.
"""


def local_gc_content(sequence, window_size):
    """
    Calculate GC content for every sliding window.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    window_size : int
        Size of the sliding window.

    Returns
    -------
    list of dict
        GC information for each window.
    """

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError(
            "Sequence cannot be empty."
        )

    if window_size <= 0:
        raise ValueError(
            "Window size must be greater than zero."
        )

    if window_size > len(sequence):
        raise ValueError(
            "Window size cannot exceed sequence length."
        )

    results = []

    for start in range(
        0,
        len(sequence) - window_size + 1
    ):

        window = sequence[
            start:start + window_size
        ]

        gc_count = (
            window.count("G")
            + window.count("C")
        )

        gc = gc_count / window_size

        results.append({
            "start": start,
            "end": start + window_size,
            "sequence": window,
            "gc_content": gc,
            "gc_percentage": gc * 100
        })

    return results


def local_gc_in_range(
    sequence,
    window_size,
    gc_min=None,
    gc_max=None
):
    """
    Check whether every sliding window satisfies
    the specified GC range.
    """

    windows = local_gc_content(
        sequence,
        window_size
    )

    for window in windows:

        gc = window["gc_content"]

        if gc_min is not None and gc < gc_min:
            return False

        if gc_max is not None and gc > gc_max:
            return False

    return True


def local_gc_violations(
    sequence,
    window_size,
    gc_min=None,
    gc_max=None
):
    """
    Return all windows that violate the GC range.
    """

    windows = local_gc_content(
        sequence,
        window_size
    )

    violations = []

    for window in windows:

        gc = window["gc_content"]

        below_min = (
            gc_min is not None
            and gc < gc_min
        )

        above_max = (
            gc_max is not None
            and gc > gc_max
        )

        if below_min or above_max:

            violations.append(window)

    return violations