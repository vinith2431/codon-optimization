"""
Sequence validation and basic DNA/protein sequence utilities.
"""

DNA_BASES = set("ATGC")


def clean_dna(sequence):
    """
    Remove whitespace and convert a DNA sequence to uppercase.
    """
    if sequence is None:
        raise ValueError("Sequence cannot be None.")

    sequence = "".join(sequence.split()).upper()

    return sequence


def validate_dna(sequence, require_complete_codons=True):
    """
    Validate a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    require_complete_codons : bool
        If True, sequence length must be divisible by 3.
        Set to False for general DNA operations.
    """

    sequence = clean_dna(sequence)

    if not sequence:
        raise ValueError(
            "DNA sequence cannot be empty."
        )

    invalid_bases = set(sequence) - DNA_BASES

    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases found: {sorted(invalid_bases)}"
        )

    if (
        require_complete_codons
        and len(sequence) % 3 != 0
    ):
        raise ValueError(
            "DNA sequence length must be divisible by 3."
        )

    return True


def get_codons(sequence):
    """
    Split a DNA sequence into codons.

    Codon-based operations require a complete sequence
    whose length is divisible by 3.
    """

    sequence = clean_dna(sequence)

    validate_dna(
        sequence,
        require_complete_codons=True
    )

    return [
        sequence[i:i + 3]
        for i in range(0, len(sequence), 3)
    ]


def gc_content(sequence):
    """
    Calculate GC content as a fraction between 0 and 1.
    """

    sequence = clean_dna(sequence)

    validate_dna(
        sequence,
        require_complete_codons=False
    )

    gc_count = (
        sequence.count("G")
        + sequence.count("C")
    )

    return gc_count / len(sequence)


def gc_percentage(sequence):
    """
    Calculate GC content as a percentage.
    """

    return gc_content(sequence) * 100


def dna_to_rna(sequence):
    """
    Convert DNA sequence to RNA.

    This operation does not require the sequence length
    to be divisible by 3.
    """

    sequence = clean_dna(sequence)

    validate_dna(
        sequence,
        require_complete_codons=False
    )

    return sequence.replace(
        "T",
        "U"
    )


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.

    This operation does not require a complete codon.
    """

    sequence = clean_dna(sequence)

    validate_dna(
        sequence,
        require_complete_codons=False
    )

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    return "".join(
        complement[base]
        for base in reversed(sequence)
    )


def sequence_length(sequence):
    """
    Return sequence length in nucleotides.
    """

    sequence = clean_dna(sequence)

    validate_dna(
        sequence,
        require_complete_codons=False
    )

    return len(sequence)