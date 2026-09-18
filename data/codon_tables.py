"""
Codon usage tables for the organisms used in the project.

Each amino acid maps to its synonymous codons and their usage frequencies.
These tables are used by the codon optimization algorithms.
"""

CODON_TABLES = {

    "ecoli": {
        "F": {"TTT": 0.58, "TTC": 0.42},
        "L": {"TTA": 0.14, "TTG": 0.13, "CTT": 0.12, "CTC": 0.10,
              "CTA": 0.04, "CTG": 0.47},
        "I": {"ATT": 0.49, "ATC": 0.39, "ATA": 0.11},
        "M": {"ATG": 1.00},
        "V": {"GTT": 0.28, "GTC": 0.20, "GTA": 0.17, "GTG": 0.35},
        "S": {"TCT": 0.17, "TCC": 0.15, "TCA": 0.14, "TCG": 0.14,
              "AGT": 0.16, "AGC": 0.24},
        "P": {"CCT": 0.18, "CCC": 0.13, "CCA": 0.20, "CCG": 0.50},
        "T": {"ACT": 0.19, "ACC": 0.40, "ACA": 0.17, "ACG": 0.24},
        "A": {"GCT": 0.18, "GCC": 0.26, "GCA": 0.23, "GCG": 0.33},
        "Y": {"TAT": 0.59, "TAC": 0.41},
        "H": {"CAT": 0.57, "CAC": 0.43},
        "Q": {"CAA": 0.34, "CAG": 0.66},
        "N": {"AAT": 0.47, "AAC": 0.53},
        "K": {"AAA": 0.74, "AAG": 0.26},
        "D": {"GAT": 0.63, "GAC": 0.37},
        "E": {"GAA": 0.68, "GAG": 0.32},
        "C": {"TGT": 0.46, "TGC": 0.54},
        "W": {"TGG": 1.00},
        "R": {"CGT": 0.36, "CGC": 0.37, "CGA": 0.07,
              "CGG": 0.11, "AGA": 0.07, "AGG": 0.04},
        "G": {"GGT": 0.34, "GGC": 0.37, "GGA": 0.11, "GGG": 0.17},
        "*": {"TAA": 0.61, "TAG": 0.09, "TGA": 0.30},
    },

    "bsubtilis": {
        "F": {"TTT": 0.55, "TTC": 0.45},
        "L": {"TTA": 0.14, "TTG": 0.14, "CTT": 0.13, "CTC": 0.05,
              "CTA": 0.04, "CTG": 0.50},
        "I": {"ATT": 0.35, "ATC": 0.45, "ATA": 0.20},
        "M": {"ATG": 1.00},
        "V": {"GTT": 0.30, "GTC": 0.21, "GTA": 0.15, "GTG": 0.34},
        "S": {"TCT": 0.16, "TCC": 0.15, "TCA": 0.14, "TCG": 0.10,
              "AGT": 0.15, "AGC": 0.30},
        "P": {"CCT": 0.16, "CCC": 0.12, "CCA": 0.18, "CCG": 0.54},
        "T": {"ACT": 0.17, "ACC": 0.38, "ACA": 0.17, "ACG": 0.28},
        "A": {"GCT": 0.17, "GCC": 0.27, "GCA": 0.22, "GCG": 0.34},
        "Y": {"TAT": 0.55, "TAC": 0.45},
        "H": {"CAT": 0.58, "CAC": 0.42},
        "Q": {"CAA": 0.30, "CAG": 0.70},
        "N": {"AAT": 0.45, "AAC": 0.55},
        "K": {"AAA": 0.72, "AAG": 0.28},
        "D": {"GAT": 0.62, "GAC": 0.38},
        "E": {"GAA": 0.67, "GAG": 0.33},
        "C": {"TGT": 0.45, "TGC": 0.55},
        "W": {"TGG": 1.00},
        "R": {"CGT": 0.36, "CGC": 0.37, "CGA": 0.08,
              "CGG": 0.10, "AGA": 0.05, "AGG": 0.04},
        "G": {"GGT": 0.32, "GGC": 0.38, "GGA": 0.12, "GGG": 0.18},
        "*": {"TAA": 0.63, "TAG": 0.10, "TGA": 0.27},
    },

    "paeruginosa": {
        "F": {"TTT": 0.59, "TTC": 0.41},
        "L": {"TTA": 0.10, "TTG": 0.14, "CTT": 0.10, "CTC": 0.09,
              "CTA": 0.05, "CTG": 0.52},
        "I": {"ATT": 0.47, "ATC": 0.44, "ATA": 0.09},
        "M": {"ATG": 1.00},
        "V": {"GTT": 0.22, "GTC": 0.24, "GTA": 0.12, "GTG": 0.42},
        "S": {"TCT": 0.14, "TCC": 0.17, "TCA": 0.12, "TCG": 0.15,
              "AGT": 0.14, "AGC": 0.28},
        "P": {"CCT": 0.15, "CCC": 0.15, "CCA": 0.15, "CCG": 0.55},
        "T": {"ACT": 0.16, "ACC": 0.40, "ACA": 0.13, "ACG": 0.31},
        "A": {"GCT": 0.15, "GCC": 0.29, "GCA": 0.18, "GCG": 0.38},
        "Y": {"TAT": 0.58, "TAC": 0.42},
        "H": {"CAT": 0.55, "CAC": 0.45},
        "Q": {"CAA": 0.32, "CAG": 0.68},
        "N": {"AAT": 0.47, "AAC": 0.53},
        "K": {"AAA": 0.70, "AAG": 0.30},
        "D": {"GAT": 0.60, "GAC": 0.40},
        "E": {"GAA": 0.65, "GAG": 0.35},
        "C": {"TGT": 0.47, "TGC": 0.53},
        "W": {"TGG": 1.00},
        "R": {"CGT": 0.34, "CGC": 0.39, "CGA": 0.07,
              "CGG": 0.12, "AGA": 0.05, "AGG": 0.03},
        "G": {"GGT": 0.30, "GGC": 0.40, "GGA": 0.11, "GGG": 0.19},
        "*": {"TAA": 0.65, "TAG": 0.08, "TGA": 0.27},
    },

    "mtuberculosis": {
        "F": {"TTT": 0.62, "TTC": 0.38},
        "L": {"TTA": 0.18, "TTG": 0.14, "CTT": 0.10, "CTC": 0.06,
              "CTA": 0.06, "CTG": 0.46},
        "I": {"ATT": 0.45, "ATC": 0.43, "ATA": 0.12},
        "M": {"ATG": 1.00},
        "V": {"GTT": 0.28, "GTC": 0.23, "GTA": 0.16, "GTG": 0.33},
        "S": {"TCT": 0.20, "TCC": 0.16, "TCA": 0.14, "TCG": 0.10,
              "AGT": 0.15, "AGC": 0.25},
        "P": {"CCT": 0.20, "CCC": 0.14, "CCA": 0.22, "CCG": 0.44},
        "T": {"ACT": 0.20, "ACC": 0.37, "ACA": 0.18, "ACG": 0.25},
        "A": {"GCT": 0.20, "GCC": 0.27, "GCA": 0.23, "GCG": 0.30},
        "Y": {"TAT": 0.60, "TAC": 0.40},
        "H": {"CAT": 0.58, "CAC": 0.42},
        "Q": {"CAA": 0.35, "CAG": 0.65},
        "N": {"AAT": 0.49, "AAC": 0.51},
        "K": {"AAA": 0.76, "AAG": 0.24},
        "D": {"GAT": 0.64, "GAC": 0.36},
        "E": {"GAA": 0.69, "GAG": 0.31},
        "C": {"TGT": 0.50, "TGC": 0.50},
        "W": {"TGG": 1.00},
        "R": {"CGT": 0.32, "CGC": 0.35, "CGA": 0.10,
              "CGG": 0.12, "AGA": 0.07, "AGG": 0.04},
        "G": {"GGT": 0.34, "GGC": 0.35, "GGA": 0.13, "GGG": 0.18},
        "*": {"TAA": 0.60, "TAG": 0.10, "TGA": 0.30},
    },
}


def get_codon_table(organism):
    """
    Return the codon usage table for the requested organism.
    """
    organism = organism.lower()

    if organism not in CODON_TABLES:
        raise ValueError(
            f"Unknown organism '{organism}'. "
            f"Available organisms: {list(CODON_TABLES.keys())}"
        )

    return CODON_TABLES[organism]


def get_codons(amino_acid, organism):
    """
    Return all synonymous codons for an amino acid.
    """
    table = get_codon_table(organism)
    amino_acid = amino_acid.upper()

    if amino_acid not in table:
        raise ValueError(f"Unknown amino acid: {amino_acid}")

    return table[amino_acid]


def get_preferred_codon(amino_acid, organism):
    """
    Return the highest-frequency codon for an amino acid.
    """
    codons = get_codons(amino_acid, organism)

    return max(codons, key=codons.get)