# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Create aligner
    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-4, gap_extend=-1)

    species = [
        ("Gallus gallus", gg_seq),
        ("Mus musculus", mm_seq),
        ("Balaeniceps rex", br_seq),
        ("Tursiops truncatus", tt_seq),
    ]

    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    results = []
    for name, seq in species:
        score, _, _ = nw.align(hs_seq, seq)
        results.append((score, name))

    results.sort(key=lambda x: x[0], reverse=True)

    print("Species ordered by similarity to human BRD2:")
    for score, name in results:
        print(f"{name}\t{score}")

    print()

    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("Alignment scores vs human BRD2:")
    for score, name in results:
        print(f"Human vs {name}: {score}")

    pass
    

if __name__ == "__main__":
    main()
