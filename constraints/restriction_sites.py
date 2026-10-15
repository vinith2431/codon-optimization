"""
Restriction-site constraint utilities.
"""

from core.sequence import reverse_complement

# Documented standard cloning and assembly enzymes commonly avoided in E. coli expression vectors
COMMON_CLONING_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "NcoI": "CCATGG",
    "BsaI": "GGTCTC",
}


def find_restriction_sites(
    sequence,
    forbidden_sites,
    check_reverse_complement=False
):
    """
    Find occurrences of forbidden restriction sites.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    forbidden_sites : list[str] or dict
        Recognition sequences that should not occur, or dict of enzyme names to sequences.

    check_reverse_complement : bool
        If True, also search the reverse complement for non-palindromic recognition sites.

    Returns
    -------
    list[dict]
        Locations of forbidden sites with site, start, end, and strand.
    """

    sequence = "".join(sequence.split()).upper()

    if forbidden_sites is None:
        return []

    # If passed a dict (like COMMON_CLONING_ENZYMES), use values
    if isinstance(forbidden_sites, dict):
        site_list = list(forbidden_sites.values())
    else:
        site_list = list(forbidden_sites)

    results = []

    for site in site_list:

        site = site.upper().strip()

        if not site:
            continue

        # Forward strand search (overlapping matches supported)
        start = 0
        while True:
            position = sequence.find(site, start)
            if position == -1:
                break

            results.append({
                "site": site,
                "start": position,
                "end": position + len(site),
                "strand": "+"
            })
            start = position + 1

        # Search reverse complement if requested and not self-complementary (palindromic)
        if check_reverse_complement:
            try:
                rc_site = reverse_complement(site)
            except Exception:
                rc_site = site

            if rc_site != site:
                start = 0
                while True:
                    position = sequence.find(rc_site, start)
                    if position == -1:
                        break

                    results.append({
                        "site": f"{site}(RC:{rc_site})",
                        "start": position,
                        "end": position + len(rc_site),
                        "strand": "-"
                    })
                    start = position + 1

    results.sort(key=lambda x: (x["start"], x["end"]))
    return results


def contains_restriction_site(
    sequence,
    forbidden_sites,
    check_reverse_complement=False
):
    """
    Return True if at least one forbidden site occurs.
    """

    return len(
        find_restriction_sites(
            sequence,
            forbidden_sites,
            check_reverse_complement=check_reverse_complement
        )
    ) > 0


def restriction_constraint_satisfied(
    sequence,
    forbidden_sites,
    check_reverse_complement=False
):
    """
    Return True when no forbidden restriction site
    occurs in the sequence.
    """

    return not contains_restriction_site(
        sequence,
        forbidden_sites,
        check_reverse_complement=check_reverse_complement
    )