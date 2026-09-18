"""
Homopolymer constraint utilities.
"""


def find_homopolymers(
    sequence,
    max_length=5
):
    """
    Find homopolymer runs longer than the allowed length.

    Example:
        max_length=5

        AAAAAA -> violation
        AAAAA  -> allowed
    """

    sequence = "".join(sequence.split()).upper()

    if max_length <= 0:
        raise ValueError(
            "max_length must be greater than zero."
        )

    violations = []

    if not sequence:
        return violations

    current_base = sequence[0]
    start = 0

    for i in range(1, len(sequence) + 1):

        if (
            i == len(sequence)
            or sequence[i] != current_base
        ):

            run_length = i - start

            if run_length > max_length:

                violations.append({
                    "base": current_base,
                    "start": start,
                    "end": i,
                    "length": run_length,
                    "sequence": current_base * run_length
                })

            if i < len(sequence):

                current_base = sequence[i]
                start = i

    return violations


def has_homopolymer(
    sequence,
    max_length=5
):
    """
    Return True if an excessive homopolymer exists.
    """

    return len(
        find_homopolymers(
            sequence,
            max_length
        )
    ) > 0


def homopolymer_constraint_satisfied(
    sequence,
    max_length=5
):
    """
    Return True if the sequence satisfies the
    homopolymer constraint.
    """

    return not has_homopolymer(
        sequence,
        max_length
    )