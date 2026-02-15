# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    seq1, _ = read_fasta("../data/test_seq1.fa")
    seq2, _ = read_fasta("../data/test_seq2.fa")

    nw = NeedlemanWunsch("../substitution_matrices/BLOSUM62.mat", gap_open=-4, gap_extend=-1)
    nw.align(seq1, seq2)

    n, m = len(seq1), len(seq2)
    M = nw._align_matrix

    # basic shape + corner
    assert M.shape == (n + 1, m + 1)
    assert M[0, 0] == pytest.approx(0.0)

    # check borders
    assert M[1, 0] == pytest.approx(-4.0)
    assert M[0, 1] == pytest.approx(-4.0)

    if n >= 2:
        assert M[2, 0] == pytest.approx(-8.0)
    if m >= 2:
        assert M[0, 2] == pytest.approx(-8.0)

def test_nw_backtrace():

    seq3, _ = read_fasta("../data/test_seq3.fa")
    seq4, _ = read_fasta("../data/test_seq4.fa")

    nw = NeedlemanWunsch("../substitution_matrices/BLOSUM62.mat", gap_open=-4, gap_extend=-1)

    score, a_aln, b_aln = nw.align(seq3, seq4)

    assert score == pytest.approx(18.0)
    assert a_aln == "MAVHQLIRRP"
    assert b_aln == "M---QLIRHP"