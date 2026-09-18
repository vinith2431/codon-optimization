"""
Restriction-site constraint utilities.
"""


def find_restriction_sites(
    sequence,
    forbidden_sites
):
    """
    Find occurrences of forbidden restriction sites.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    forbidden_sites : list[str]
        Recognition sequences that should not occur.

    Returns
    -------
    list[dict]
        Locations of forbidden sites.
    """

    sequence = "".join(sequence.split()).upper()

    if forbidden_sites is None:
        return []

    results = []

    for site in forbidden_sites:

        site = site.upper()

        if not site:
            continue

        start = 0

        while True:

            position = sequence.find(
                site,
                start
            )

            if position == -1:
                break

            results.append({
                "site": site,
                "start": position,
                "end": position + len(site)
            })

            # +1 allows overlapping matches
            start = position + 1

    return results


def contains_restriction_site(
    sequence,
    forbidden_sites
):
    """
    Return True if at least one forbidden site occurs.
    """

    return len(
        find_restriction_sites(
            sequence,
            forbidden_sites
        )
    ) > 0


def restriction_constraint_satisfied(
    sequence,
    forbidden_sites
):
    """
    Return True when no forbidden restriction site
    occurs in the sequence.
    """

    return not contains_restriction_site(
        sequence,
        forbidden_sites
    )