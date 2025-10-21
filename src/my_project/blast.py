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



    # Save the matrix to an absolute path accessible to your system.
# blast.py
from typing import Callable, Optional, Dict

def process_result_matrix(
    info_matrix,
    result_matrix,
    do_blast: bool = True,
    on_log: Optional[Callable[[str], None]] = None,
) -> Dict[int, str]:
    """
    For each [orf_code, left, right] row in result_matrix:
      - map orf_code → a protein row index
      - extract a protein subsequence
      - if do_blast: run BLAST, return top name per row
    Returns: {row_index: "Top hit name (E=...)"}
    """
    def log(msg: str) -> None:
        if on_log:
            on_log(msg)   # goes to Streamlit UI via your queue
        else:
            print(msg)    # fallback to terminal

    a_map = {1: 14, 2: 15, 3: 16, 4: 20, 5: 21, 6: 22}
    blast_names: Dict[int, str] = {}

    # NOTE: numpy ints → coerce to int to avoid weird keys
    for idx, (orf_code, left, right) in enumerate(result_matrix):
        try:
            a_row = a_map[int(orf_code)]
        except KeyError:
            log(f"Row {idx}: unknown ORF code {orf_code}, skipping.")
            continue

        protein_seq = extract_protein_seq(info_matrix, a_row, int(left), int(right))
        log(f"\nRow {idx} (ORF {int(orf_code)}): cols {int(left)}–{int(right)}")
        log(f"Protein sequence: {protein_seq[:60]}{'...' if len(protein_seq) > 60 else ''}")

        if do_blast and protein_seq:
            log("Submitting to BLAST…")
            record = blast_protein_seq(protein_seq)

            if getattr(record, "alignments", None):
                top = record.alignments[0]
                hsp = top.hsps[0]
                name = getattr(top, "hit_def", "") or getattr(top, "title", "") or str(top)
                blast_names[idx] = f"{name} (E={getattr(hsp, 'expect', 'NA')})"

                log("Top 5 BLAST hits for this query:")
                for aln in record.alignments[:5]:
                    h = aln.hsps[0]
                    n2 = getattr(aln, "hit_def", "") or getattr(aln, "title", "") or aln.accession
                    log(f"  {aln.accession}\t{n2}\tE={getattr(h, 'expect', 'NA')}")
            else:
                log("No BLAST hits found.")
        elif do_blast:
            log("Empty sequence—skipping BLAST.")

    return blast_names




def blast_protein_seq(seq, email=None):
        """
        Run NCBIWWW.qblast on a single protein sequence (BLASTP against nr),
        returning the parsed record.
        """
        fasta = f">query\n{seq}"
        handle = NCBIWWW.qblast("blastp", "nr", fasta, format_type="XML")
        return NCBIXML.read(handle)

