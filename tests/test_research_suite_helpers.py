import pytest
from experiments.run_research_suite import max_kmer_frequency, five_prime_gc, longest_homopolymer, read_fasta

def test_overlapping_kmer_frequency():
    assert max_kmer_frequency('AAAAAA', 3) == 1.0
    assert max_kmer_frequency('ACGTACGT', 3) == pytest.approx(2 / 6)

def test_kmer_size_validation():
    with pytest.raises(ValueError): max_kmer_frequency('ACGT', 0)

def test_five_prime_gc_prefix():
    assert five_prime_gc('GGGGAAAATTTT', 4) == 1.0
    assert five_prime_gc('GGGGAAAATTTT', 8) == pytest.approx(0.5)

def test_longest_homopolymer():
    assert longest_homopolymer('') == 0
    assert longest_homopolymer('ACCCCGG') == 4
    assert longest_homopolymer('ACGT') == 1

def test_fasta_requires_coding_length(tmp_path):
    p=tmp_path/'bad.fasta'; p.write_text('>x\nATGAAAT\n', encoding='utf-8')
    with pytest.raises(ValueError, match='divisible by 3'): read_fasta(p)

def test_fasta_parsing(tmp_path):
    p=tmp_path/'ok.fasta'; p.write_text('>x\nATG\n AAA\n', encoding='utf-8')
    assert read_fasta(p) == 'ATGAAA'
