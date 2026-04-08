# src/my_project/pipeline.py
from __future__ import annotations
from typing import Callable, Dict, Any, List, Tuple, Optional
import numpy as np
import pandas as pd

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
#prodigal
from my_project import microservices as ms_pipeline


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


# ---------------------------------------------------------------------------
# Isolation-Forest anomaly detection helpers
# ---------------------------------------------------------------------------

def _gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    total = len(seq) - seq.count('N')
    return gc / total if total > 0 else 0.0


def _count_orfs(orf_rows) -> int:
    total = 0
    for row in orf_rows:
        total += sum(1 for aa in row if aa not in ('', '-', None, '*'))
    return total


def _similarity(ref_seq: str, evo_seq: str) -> float:
    matches = sum(1 for a, b in zip(ref_seq, evo_seq) if a == b and a != '-')
    total   = sum(1 for a, b in zip(ref_seq, evo_seq) if a != '-' and b != '-')
    return matches / total if total > 0 else 0.0


def _effective_length(seq: str) -> int:
    return sum(1 for c in seq if c != '-')


def run_isolation_forest(
    info_matrix: np.ndarray,
    n_sections: int,
    contamination: str | float = "auto",
    n_estimators: int = 200,
    random_state: int = 42,
) -> Dict[str, Any]:
    """
    Compute per-chunk features from *info_matrix* and fit an IsolationForest.

    Returns a dict with:
        feature_df   – DataFrame of raw features + anomaly_score + label
        anomaly_scores – np.ndarray in [0,1], 1 = most anomalous
        anomaly_labels – np.ndarray of {-1, 1}  (-1 = anomaly)
    """
    from sklearn.ensemble import IsolationForest
    from sklearn.preprocessing import StandardScaler

    chunk_boundaries = np.array_split(np.arange(info_matrix.shape[1]), n_sections)

    feature_df = pd.DataFrame(
        index=range(n_sections),
        columns=['similarity_score', 'length_ratio', 'gc_difference', 'orf_density_difference'],
        dtype=float,
    )

    for n, cols in enumerate(chunk_boundaries):
        chunk = info_matrix[:, cols]
        chunk_len = len(cols)

        ref_seq = "".join(str(v) for v in chunk[1].tolist())
        evo_seq = "".join(str(v) for v in chunk[3].tolist())

        ref_len = _effective_length(ref_seq)
        evo_len = _effective_length(evo_seq)

        feature_df.loc[n, 'similarity_score']      = _similarity(ref_seq, evo_seq)
        feature_df.loc[n, 'length_ratio']           = ref_len / (evo_len + 1e-9)
        feature_df.loc[n, 'gc_difference']          = _gc_content(ref_seq) - _gc_content(evo_seq)

        ref_orf_density = _count_orfs([chunk[14].tolist(), chunk[15].tolist(), chunk[16].tolist()]) / chunk_len
        evo_orf_density = _count_orfs([chunk[17].tolist(), chunk[18].tolist(), chunk[19].tolist()]) / chunk_len
        feature_df.loc[n, 'orf_density_difference'] = ref_orf_density - evo_orf_density

    X = feature_df.values.astype(float)
    # Replace any NaN/inf that could crash sklearn
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = IsolationForest(
        n_estimators=n_estimators,
        contamination=contamination,
        random_state=random_state,
        n_jobs=-1,
    )
    labels = clf.fit_predict(X_scaled)
    scores = clf.score_samples(X_scaled)

    scores_norm   = (scores - scores.min()) / (scores.max() - scores.min() + 1e-12)
    anomaly_score = 1.0 - scores_norm

    feature_df['anomaly_score'] = anomaly_score
    feature_df['label']         = labels

    return {
        "feature_df":      feature_df,
        "anomaly_scores":  anomaly_score,
        "anomaly_labels":  labels,
        "chunk_boundaries": [c.tolist() for c in chunk_boundaries],
    }


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
    final_matrix = np.array(build_final_matrix(chunk_matrices))

    # recompute frameshift counter as true cumulative across entire sequence
    fs_cumulative = []
    current = 0
    for i in range(final_matrix.shape[1]):
        ref_base = str(final_matrix[1, i])
        evo_base = str(final_matrix[3, i])
        ref_is_gap = ref_base in ('-', 'N') or ref_base.upper() not in {'A','C','G','T','N'}
        evo_is_gap = evo_base in ('-', 'N') or evo_base.upper() not in {'A','C','G','T','N'}
        if ref_is_gap and not evo_is_gap:
            current -= 1
        elif evo_is_gap and not ref_is_gap:
            current += 1
        fs_cumulative.append(current)

    final_matrix[5, :] = fs_cumulative
    match_map = {(m["evo_index"], m["ref_index"]): m for m in matches}

    # 5) Augment with AA / reverse complements
    data = final_matrix
    ref_seq_list     = list(data[1, :])
    evolved_seq_list = list(data[3, :])

    ref_seq_list     = list(data[1, :])
    evolved_seq_list = list(data[3, :])

    ref_revcomp     = reverse_complement(ref_seq_list)
    evolved_revcomp = reverse_complement(evolved_seq_list)

    seq_len = data.shape[1]
    placeholder = [0] * seq_len  # rows 8,9,11,13 will be filled after prodigal

    data2 = np.vstack([
        data,
        placeholder,      # row 8  — ref active ORF tracker (filled below)
        placeholder,      # row 9  — evo active ORF tracker (filled below)
        ref_revcomp,      # row 10 — revcomp ref nucleotides
        placeholder,      # row 11 — revcomp ref active ORF tracker (filled below)
        evolved_revcomp,  # row 12 — revcomp evo nucleotides
        placeholder,      # row 13 — revcomp evo active ORF tracker (filled below)
    ])

    # 6) Normalize dashes to N on specified rows
    N_matrix = data2.copy()
    for ridx in (1, 3, 10, 12):
        if ridx < N_matrix.shape[0]:
            N_matrix[ridx, :] = np.where(N_matrix[ridx, :] == '-', 'N', N_matrix[ridx, :])
    log("Running Prodigal ORF prediction…")

    tracker_seq_map = {
        "forward_ref":     (1,  8),
        "forward_evolved": (3,  9),
        "revcomp_ref":     (10, 11),
        "revcomp_evolved": (12, 13),
    }

    seq_names = list(tracker_seq_map.keys())

    for seq_idx, (seq_name, (nuc_row_idx, tracker_row_idx)) in enumerate(tracker_seq_map.items()):
        p("prodigal", seq_idx, 4)  # emit progress before starting each sequence
        log(f"Prodigal: processing {seq_name}…")

        nuc_row = data2[nuc_row_idx, :]

        clean_seq = ""
        pos_map   = []
        for i, ch in enumerate(nuc_row):
            if ch not in ('-', 'N'):
                clean_seq += ch
                pos_map.append(i)

        try:
            orfs = ms_pipeline.run_prodigal(clean_seq, seq_name=seq_name)
            log(f"Prodigal: {len(orfs)} ORFs found in {seq_name}.")
        except Exception as e:
            log(f"Prodigal failed for {seq_name}: {e} — skipping tracker.")
            orfs = []

        tracker = ms_pipeline.build_active_orf_tracker(
            seq     = clean_seq,
            pos_map = pos_map,
            seq_len = seq_len,
            orfs    = orfs,
        )

        data2[tracker_row_idx, :] = tracker

    p("prodigal", 4, 4)  # complete
    log("Prodigal ORF prediction complete.")

    for tracker_row in (8, 9, 11, 13):
        N_matrix[tracker_row, :] = data2[tracker_row, :]
    # 7) Protein annotation rows (12 more)

    seq_map_translation = [
        (1,  "forward_ref"),
        (3,  "forward_evolved"),
        (10, "revcomp_ref"),
        (12, "revcomp_evolved"),
    ]

    protein_matrix_rows: List[np.ndarray] = []
    for nuc_row_idx, seq_name in seq_map_translation:
        nuc_row = data2[nuc_row_idx, :]

        clean_seq = ""
        pos_map = []
        for i, ch in enumerate(nuc_row):
            if ch not in ('-', 'N'):
                clean_seq += ch
                pos_map.append(i)

        for frame_offset in range(3):
            aa_array = ['X'] * seq_len
            for codon_idx in range(frame_offset, len(clean_seq) - 2, 3):
                codon = clean_seq[codon_idx:codon_idx + 3]
                aa = CODON_TABLE.get(codon.upper(), None)
                if aa is None:
                    translated = 'X'
                elif aa == '*':
                    translated = '*'
                else:
                    translated = aa
                for offset in range(3):
                    clean_pos = codon_idx + offset
                    if clean_pos < len(pos_map):
                        aa_array[pos_map[clean_pos]] = translated
            protein_matrix_rows.append(aa_array)

    protein_matrix = np.array(protein_matrix_rows)

    info_matrix = np.vstack((N_matrix, protein_matrix))
# --- BLAST (names back to UI) ---
    blast_names = {}
    try:
        blast_names = process_result_matrix(
            info_matrix,
            np.empty((0, 3), dtype=int),
            do_blast=bool(do_blast),
            on_log=on_log,
        ) or {}
    except Exception as e:
        if on_log:
            on_log(f"BLAST step skipped/failed: {e}")

        ...
    # 8) ORF diffs → result_matrix
    all_mutated_orfs = []
    try:
        log("Finding mutated ORFs…")
        forward_mutated = ms_pipeline.find_mutated_orfs(info_matrix, forward=True)
        revcomp_mutated = ms_pipeline.find_mutated_orfs(info_matrix, forward=False)
        all_mutated_orfs = forward_mutated + revcomp_mutated
        log(f"Found {len(all_mutated_orfs)} mutated ORF pairs "
            f"({sum(1 for r in all_mutated_orfs if r['pair_tier']==1)} exact, "
            f"{sum(1 for r in all_mutated_orfs if r['pair_tier']==2)} overlapping, "
            f"{sum(1 for r in all_mutated_orfs if r['pair_tier']==3)} orphaned).")
    except Exception as e:
        import traceback
        log(f"find_mutated_orfs failed: {e}\n{traceback.format_exc()}")
# test
    all_results = []
    for entry in all_mutated_orfs:
        is_revcomp = not entry.get("forward", True)
        ref_orf = entry.get("ref_orf")
        evo_orf = entry.get("evo_orf")
        tier = entry["pair_tier"]

        # code = ref frame (1-3 forward, 4-6 revcomp), -1 if no ref orf
        if ref_orf is not None:
            frame, ref_start, ref_end = ref_orf
            code = frame + (3 if is_revcomp else 0)
        else:
            frame, evo_start_tmp, _ = evo_orf
            code = frame + (3 if is_revcomp else 0)

        ref_start = ref_orf[1] if ref_orf is not None else -1
        ref_end   = ref_orf[2] if ref_orf is not None else -1
        evo_start = evo_orf[1] if evo_orf is not None else -1
        evo_end   = evo_orf[2] if evo_orf is not None else -1

        all_results.append([code, ref_start, ref_end, evo_start, evo_end, tier])

    result_matrix = np.array(all_results, dtype=int) if all_results else np.empty((0, 6), dtype=int)

    # 9) Isolation-Forest anomaly detection
    isolation_forest_results: Optional[Dict[str, Any]] = None
    isolation_forest_error: Optional[str] = None
    try:
        log("Running Isolation Forest anomaly detection…")
        isolation_forest_results = run_isolation_forest(info_matrix, n_sections)
        n_anomalies = int((isolation_forest_results["anomaly_labels"] == -1).sum())
        log(f"Isolation Forest complete: {n_anomalies}/{n_sections} anomalous chunks detected.")
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        log(f"find_mutated_orfs failed: {e}\n{tb}")
        print(f"CRITICAL ERROR IN find_mutated_orfs:\n{tb}", flush=True)

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
        "blast_names": blast_names,
        "isolation_forest": isolation_forest_results,
        "isolation_forest_error": isolation_forest_error,
        "mutated_orfs": all_mutated_orfs,
    }

