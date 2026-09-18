"""
Repeated-sequence constraint utilities.
"""


def find_repeats(sequence, min_length=6):
    """
    Find repeated subsequences of at least min_length.

    Returns a list of dictionaries describing repeated
    sequence pairs.
    """

    sequence = "".join(sequence.split()).upper()

    if min_length <= 0:
        raise ValueError(
            "min_length must be greater than zero."
        )

    repeats = []

    seen = {}

    for i in range(len(sequence) - min_length + 1):

        fragment = sequence[
            i:i + min_length
        ]

        if fragment in seen:

            repeats.append({
                "sequence": fragment,
                "first_position": seen[fragment],
                "second_position": i
            })

        else:
            seen[fragment] = i

    return repeats


def has_repeats(sequence, min_length=6):
    """
    Return True if repeated subsequences are present.
    """

    return len(
        find_repeats(
            sequence,
            min_length
        )
    ) > 0


def repeat_constraint_satisfied(
    sequence,
    min_length=6
):
    """
    Return True when no repeated subsequences of the
    specified length are present.
    """

    return not has_repeats(
        sequence,
        min_length
    )