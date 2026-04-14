from Bio.Blast import NCBIWWW, NCBIXML

def extract_protein_seq(info_matrix, a, b, c, tracker_row=None, frame_num=None):
    AA_triplocate = info_matrix[a, b:c]
    AA_seq = AA_triplocate[0::3]

    if tracker_row is not None and frame_num is not None:
        from my_project.microservices import frame_active
        trk_slice = info_matrix[tracker_row, b:c:3]
        result = []
        for aa, trk_val in zip(AA_seq, trk_slice):
            if not frame_active(int(trk_val), frame_num):
                break
            if aa not in ('-', 'X', '*'):
                result.append(aa)
        return ''.join(result)

    return ''.join(ch for ch in AA_seq if ch not in ('-', 'X', '*'))



    # Save the matrix to an absolute path accessible to your system.
# blast.py
from typing import Callable, Optional, Dict

def process_result_matrix(
    info_matrix,
    result_matrix,
    do_blast: bool = True,
    on_log: Optional[Callable[[str], None]] = None,
) -> Dict[int, str]:
    def log(msg: str) -> None:
        if on_log:
            on_log(msg)
        else:
            print(msg)

    a_map     = {1: 14, 2: 15, 3: 16, 4: 20, 5: 21, 6: 22}
    evo_a_map = {1: 17, 2: 18, 3: 19, 4: 23, 5: 24, 6: 25}
    trk_map     = {1: (8, 1),  2: (8, 2),  3: (8, 3),  4: (11, 1), 5: (11, 2), 6: (11, 3)}
    evo_trk_map = {1: (9, 1),  2: (9, 2),  3: (9, 3),  4: (13, 1), 5: (13, 2), 6: (13, 3)}

    blast_names: Dict[int, str] = {}

    for idx, (orf_code, left, right) in enumerate(result_matrix):
        code = int(orf_code)
        try:
            a_row     = a_map[code]
            evo_a_row = evo_a_map[code]
        except KeyError:
            log(f"Row {idx}: unknown ORF code {orf_code}, skipping.")
            continue

        trk_row,     frame_num  = trk_map.get(code,     (8,  1))
        evo_trk_row, evo_frame  = evo_trk_map.get(code, (9,  1))

        ref_seq = extract_protein_seq(
            info_matrix, a_row, int(left), int(right),
            tracker_row=trk_row, frame_num=frame_num
        )
        evo_seq = extract_protein_seq(
            info_matrix, evo_a_row, int(left), int(right),
            tracker_row=evo_trk_row, frame_num=evo_frame
        )

        log(f"\nRow {idx} (ORF {code}): cols {int(left)}–{int(right)}")
        log(f"Ref protein:  {ref_seq[:60]}{'...' if len(ref_seq) > 60 else ''}")
        log(f"Evo protein:  {evo_seq[:60]}{'...' if len(evo_seq) > 60 else ''}")

        if do_blast:
            results = []
            for label, seq in [("REF", ref_seq), ("EVO", evo_seq)]:
                if not seq:
                    log(f"{label}: empty sequence — skipping BLAST.")
                    continue
                log(f"Submitting {label} to BLAST…")
                record = blast_protein_seq(seq)
                if getattr(record, "alignments", None):
                    top = record.alignments[0]
                    hsp = top.hsps[0]
                    name = getattr(top, "hit_def", "") or getattr(top, "title", "") or str(top)
                    hit = f"{label}: {name} (E={getattr(hsp, 'expect', 'NA')})"
                    results.append(hit)
                    log(f"Top 5 BLAST hits for {label}:")
                    for aln in record.alignments[:5]:
                        h = aln.hsps[0]
                        n2 = getattr(aln, "hit_def", "") or getattr(aln, "title", "") or aln.accession
                        log(f"  {aln.accession}\t{n2}\tE={getattr(h, 'expect', 'NA')}")
                else:
                    log(f"{label}: no BLAST hits found.")

            if results:
                blast_names[idx] = "\n".join(results)

    return blast_names




def blast_protein_seq(seq, email=None):
        """
        Run NCBIWWW.qblast on a single protein sequence (BLASTP against nr),
        returning the parsed record.
        """
        fasta = f">query\n{seq}"
        handle = NCBIWWW.qblast("blastp", "nr", fasta, format_type="XML")
        return NCBIXML.read(handle)

