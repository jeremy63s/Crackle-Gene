from Bio.Blast import NCBIWWW, NCBIXML

def extract_protein_seq(info_matrix, a, b, c):
    """
    Given info_matrix and slice indices a (row), b (start col), c (end col),
    return the protein_seq from that ORF.
    """
    # slice out the triplet letters
    AA_triplocate = info_matrix[a, b:c]
    # take every 3rd position, starting at 0
    AA_seq = AA_triplocate[0::3]
    # join into a string
    return ''.join(AA_seq)

def process_result_matrix(info_matrix, result_matrix, do_blast=True):
        """
        For each row [orf_code, left, right] in result_matrix:
          - map orf_code → a row index
          - extract protein_seq
          - optionally BLAST and print top hits
        """
        a_map = {
            1: 14,
            2: 15,
            3: 16,
            4: 20,
            5: 21,
            6: 22
        }
        for idx, (orf_code, left, right) in enumerate(result_matrix):
            try:
                a_row = a_map[int(orf_code)]
            except KeyError:
                print(f"Row {idx}: unknown ORF code {orf_code}, skipping.")
                continue

            # extract
            protein_seq = extract_protein_seq(info_matrix, a_row, int(left), int(right))
            print(f"\nRow {idx} (ORF {orf_code}): cols {left}–{right}")
            print("Protein sequence:", protein_seq)

            if do_blast and protein_seq:
                print("Submitting to BLAST…")
                record = blast_protein_seq(protein_seq)
                print("Top 5 BLAST hits for this query:")
                for aln in record.alignments[:5]:
                    hsp = aln.hsps[0]
                    print(f" {aln.accession}\t{aln.hit_def}\tE={hsp.expect}")
            elif do_blast:
                print("Empty sequence—skipping BLAST.")
    # Save the matrix to an absolute path accessible to your system.

def extract_protein_seq(info_matrix, a, b, c):
    """
    Given info_matrix and slice indices a (row), b (start col), c (end col),
    return the protein_seq from that ORF.
    """
    # slice out the triplet letters
    AA_triplocate = info_matrix[a, b:c]
    # take every 3rd position, starting at 0
    AA_seq = AA_triplocate[0::3]
    # join into a string
    return ''.join(AA_seq)


def blast_protein_seq(seq, email=None):
        """
        Run NCBIWWW.qblast on a single protein sequence (BLASTP against nr),
        returning the parsed record.
        """
        fasta = f">query\n{seq}"
        handle = NCBIWWW.qblast("blastp", "nr", fasta, format_type="XML")
        return NCBIXML.read(handle)


def process_result_matrix(info_matrix, result_matrix, do_blast=True):
    """
    For each row [orf_code, left, right] in result_matrix:
      - map orf_code → a row index
      - extract protein_seq
      - optionally BLAST and print top hits
    """
    a_map = {
        1: 14,
        2: 15,
        3: 16,
        4: 20,
        5: 21,
        6: 22
    }
    for idx, (orf_code, left, right) in enumerate(result_matrix):
        try:
            a_row = a_map[int(orf_code)]
        except KeyError:
            print(f"Row {idx}: unknown ORF code {orf_code}, skipping.")
            continue

        # extract
        protein_seq = extract_protein_seq(info_matrix, a_row, int(left), int(right))
        print(f"\nRow {idx} (ORF {orf_code}): cols {left}–{right}")
        print("Protein sequence:", protein_seq)

        if do_blast and protein_seq:
            print("Submitting to BLAST…")
            record = blast_protein_seq(protein_seq)
            print("Top 5 BLAST hits for this query:")
            for aln in record.alignments[:5]:
                hsp = aln.hsps[0]
                print(f" {aln.accession}\t{aln.hit_def}\tE={hsp.expect}")
        elif do_blast:
            print("Empty sequence—skipping BLAST.")
# Save the matrix to an absolute path accessible to your system.