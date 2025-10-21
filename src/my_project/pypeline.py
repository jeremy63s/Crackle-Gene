# src/my_project/pipeline.py
from __future__ import annotations
from typing import Callable, Dict, Any, List, Tuple
import numpy as np

# Always use package-qualified imports to avoid duplicate modules
from my_project.sequence_utils import (
    divide_sequence_with_offset, reverse_complement,
    translate_nuc_row_active, create_annotation_array,
    CODON_TABLE, array_to_str, parse_orf_frame, process_orf_array
)
from my_project.alignment import diagonal_check, compare_sections, get_alignment
from my_project.matrix_utils import build_chunk_matrix, build_final_matrix
from my_project.io_helpers import extract_extended_context
from my_project.blast import process_result_matrix
import orfipy_core
from Bio.Seq import Seq
try:
    from tqdm import tqdm
except Exception:  # tqdm missing is fine
    tqdm = None

class TqdmToCallback:
    """
    Wraps a tqdm progress bar but also calls a callback with (phase, n, total)
    every time the bar updates. If tqdm isn't available, still calls the callback.
    """
    def __init__(self, phase: str, total: int, cb: ProgressCB | None, **tqdm_kwargs):
        self.phase = phase
        self.cb = cb
        self.total = int(total) if total else 1
        self.n = 0
        self._tqdm = None
        if tqdm is not None:
            self._tqdm = tqdm(total=self.total, desc=tqdm_kwargs.get("desc", phase), leave=False)

    def update(self, n: int = 1):
        self.n += int(n)
        if self._tqdm is not None:
            self._tqdm.update(n)
        if self.cb:
            self.cb(self.phase, self.n, self.total)

    def close(self):
        if self._tqdm is not None:
            self._tqdm.close()
        # ensure final tick is shown
        if self.cb and self.n < self.total:
            self.cb(self.phase, self.total, self.total)

ProgressCB = Callable[[str, int, int], None]  # phase, i, total
LogCB      = Callable[[str], None]

def run_pipeline(

    rotated_ref: str,
    rotated_evolved: str,
    n_sections: int,
    threshold: float,
    do_full: bool,
    on_progress: ProgressCB | None = None,
    on_log: LogCB | None = None,
    do_blast: bool = False,
) -> Dict[str, Any]:
    """Pure pipeline. No Streamlit. Returns data only."""
    def p(phase: str, i: int, total: int) -> None:
        if on_progress: on_progress(phase, i, total)
    def log(msg: str) -> None:
        if on_log: on_log(msg)

    # 1) Divide sequences
    ref_sections = divide_sequence_with_offset(rotated_ref, n_sections)
    evo_sections = divide_sequence_with_offset(rotated_evolved, n_sections)

    # 2) Diagonal pre-check
    # 2) Diagonal pre-check
    pct_met, max_consec_fail, raw_scores, diag_scores = diagonal_check(
        evo_sections,
        ref_sections,
        threshold,
        progress_cb=(lambda i, total: on_progress and on_progress("diagonal", i, total))  # <-- NEW
    )


    log(f"Diagonal: {pct_met:.1f}% met; longest fail run {max_consec_fail}")

    # 3) Matches
    matches: List[Dict[str, Any]] = []

    if do_full:
        total_pairs = max(1, len(evo_sections) * len(ref_sections))

        # Mirror terminal tqdm into Streamlit: pass a lambda that relays each pair
        matches = compare_sections(
            evo_sections, ref_sections, threshold,
            rotated_evolved, rotated_ref,
            progress_cb=lambda k, total: on_progress and on_progress("compare", k, total)
        )
        # If compare_sections has inner loops, best is to pass a step callback down into it.
        # For now, mirror coarse progress so UI moves (start+finish).
        total_pairs = max(1, len(evo_sections) * len(ref_sections))
        # a tqdm bar that also ticks the Streamlit bar:
        bar = TqdmToCallback("compare", total_pairs, on_progress, desc="compare")
        # compare_sections computes internally; we can't tick per pair without modifying it:
        #p("compare", 0, total_pairs)  # 0%
        #matches = compare_sections(
        #    evo_sections, ref_sections, threshold,
        #    rotated_evolved, rotated_ref
        #)
        # jump to 100%
        #bar.update(total_pairs)  # mirrors terminal & UI to full
        #bar.close()
    else:
        total = len(evo_sections)
        bar = TqdmToCallback("compare", max(1, total), on_progress, desc="compare")
        for i, norm_score in enumerate(diag_scores, start=1):
            bar.update(1)  # <- advances terminal tqdm AND calls on_progress("compare", i, total)

            if norm_score < threshold:
                continue

            raw = raw_scores[i - 1]
            evo_sec, evo_off = evo_sections[i - 1]
            ref_sec, ref_off = ref_sections[i - 1]
            (aligned_evo, aligned_ref, midline,
             align_score, bq, eq, br, er) = get_alignment(evo_sec, ref_sec)
            ext_evo = extract_extended_context(rotated_evolved, evo_off, bq, eq)
            ext_ref = extract_extended_context(rotated_ref,   ref_off, br, er)
            matches.append({
                "evo_index": i,
                "ref_index": i,
                "raw_alignment_score": raw,
                "normalized_score": norm_score,
                "alignment_score": align_score,
                "aligned_evo": aligned_evo,
                "aligned_ref": aligned_ref,
                "midline": midline,
                "begin_query": bq, "end_query": eq,
                "begin_ref": br,   "end_ref": er,
                "evo_offset": evo_off, "ref_offset": ref_off,
                "extended_evo": ext_evo, "extended_ref": ext_ref,
                "is_match": 1, "match_index": len(matches),
            })
        bar.close()

    # 4) Build matrices
    chunk_matrices: List[np.ndarray] = []
    matrices_total = max(1, len(matches))
    bar_m = TqdmToCallback("matrices", matrices_total, on_progress, desc="matrices")
    for _j, m in enumerate(matches, start=1):
        bar_m.update(1)  # <- syncs terminal + UI
        chunk_matrices.append(build_chunk_matrix(m))
    bar_m.close()

    final_matrix = np.array(build_final_matrix(chunk_matrices))
    match_map = {(m["evo_index"], m["ref_index"]): m for m in matches}

    # 5) Augment with AA / reverse complements
    data = final_matrix
    ref_seq_list     = list(data[1, :])
    evolved_seq_list = list(data[3, :])

    aa_ref             = translate_nuc_row_active(ref_seq_list, CODON_TABLE)
    aa_evolved         = translate_nuc_row_active(evolved_seq_list, CODON_TABLE)
    ref_revcomp        = reverse_complement(ref_seq_list)
    aa_ref_revcomp     = translate_nuc_row_active(ref_revcomp, CODON_TABLE)
    evolved_revcomp    = reverse_complement(evolved_seq_list)
    aa_evolved_revcomp = translate_nuc_row_active(evolved_revcomp, CODON_TABLE)

    data2 = np.vstack([data, aa_ref, aa_evolved, ref_revcomp, aa_ref_revcomp, evolved_revcomp, aa_evolved_revcomp])

    # 6) Normalize dashes to N on specified rows
    N_matrix = data2.copy()
    for ridx in (1, 3, 10, 12):
        if ridx < N_matrix.shape[0]:
            N_matrix[ridx, :] = np.where(N_matrix[ridx, :] == '-', 'N', N_matrix[ridx, :])

    # 7) Protein annotation rows (12 more)
    seq_indices = {"forward_ref": 1, "forward_evolved": 3, "revcomp_ref": 10, "revcomp_evolved": 12}
    sequences = {name: array_to_str(N_matrix[idx, :]) for name, idx in seq_indices.items()}

    results_ann = {}
    for name, seq in sequences.items():
        seq_len = len(seq)
        frame_annotations = {0: create_annotation_array(seq_len), 1: create_annotation_array(seq_len), 2: create_annotation_array(seq_len)}
        for start, stop, strand_val, description in orfipy_core.orfs(seq, minlen=3, maxlen=1_000_000, strand='f'):
            frame_val = parse_orf_frame(description)
            if frame_val is None: continue
            frame_idx = frame_val - 1
            protein_seq = str(Seq(seq[start:stop]).translate(to_stop=True))
            for i, aa in enumerate(protein_seq):
                codon_start = start + i * 3
                for pos in range(codon_start, min(codon_start + 3, seq_len)):
                    frame_annotations[frame_idx][pos] = aa
        results_ann[name] = frame_annotations

    protein_matrix_rows: List[np.ndarray] = []
    for seq_name in ["forward_ref", "forward_evolved", "revcomp_ref", "revcomp_evolved"]:
        fa = results_ann[seq_name]
        for frame in range(3):
            protein_matrix_rows.append(fa[frame])
    protein_matrix = np.array(protein_matrix_rows)

    info_matrix = np.vstack((N_matrix, protein_matrix))

    # 8) ORF diffs → result_matrix
    a_map = {1: 14, 2: 15, 3: 16, 4: 20, 5: 21, 6: 22}
    # --- ORF diffs → result_matrix ---
    orf_mapping = {
        "ORF1_forward": {"code": 1, "rows": (14, 17)},
        "ORF2_forward": {"code": 2, "rows": (15, 18)},
        "ORF3_forward": {"code": 3, "rows": (16, 19)},
        "ORF1_reverse": {"code": 4, "rows": (20, 23)},
        "ORF2_reverse": {"code": 5, "rows": (21, 24)},
        "ORF3_reverse": {"code": 6, "rows": (22, 25)},
    }

    all_results: List[Tuple[int, int, int]] = []
    for _, props in orf_mapping.items():
        code = props["code"]
        r1, r2 = props["rows"]
        # where rows differ ⇒ ORF segment
        arr = np.where(info_matrix[r1, :] != info_matrix[r2, :])[0].astype(int)
        rows = process_orf_array(info_matrix, arr, r1, r2, code)
        all_results.extend(rows)

    result_matrix = np.array(all_results) if all_results else np.empty((0, 3), dtype=int)

    # --- BLAST (names back to UI) ---
    blast_names = {}
    try:
        blast_names = process_result_matrix(
            info_matrix,
            result_matrix,
            do_blast=bool(do_blast),  # <-- forward the checkbox flag
            on_log=on_log,  # <-- forward logger so messages appear in Streamlit
        ) or {}
    except Exception as e:
        if on_log:
            on_log(f"BLAST step skipped/failed: {e}")

    # --- return everything the UI needs ---
    return {
        "ref_sections": ref_sections,
        "evo_sections": evo_sections,
        "diag_pct_met": pct_met,
        "diag_max_consec_fail": max_consec_fail,
        "matches": matches,
        "chunk_matrices": chunk_matrices,
        "final_matrix": data2,
        "match_map": match_map,
        "info_matrix": info_matrix,
        "result_matrix": result_matrix,
        "blast_names": blast_names,  # <-- Tab 2 reads names from here
    }

