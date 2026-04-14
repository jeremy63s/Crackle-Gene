# microservices.py
from __future__ import annotations
import re
from typing import Dict, Tuple
# for prodigal
import subprocess
import tempfile
import os
import numpy as np

# IUPAC DNA (includes RNA U) + gap '-'
_AMBIG = set("ACGTURYSWKMBDHVN-")

def clean_seq(s: str) -> str:
    """Uppercase, strip non-letters/dashes, convert U->T."""
    s = re.sub(r"[^A-Za-z-]", "", (s or "")).upper()
    return s.replace("U", "T")

def validate_seq(seq: str) -> Tuple[bool, str]:
    """Basic validation for DNA strings with IUPAC ambiguity codes."""
    if not seq:
        return False, "Sequence is empty."
    bad = {c for c in seq if c not in _AMBIG}
    if bad:
        return False, f"Invalid characters: {''.join(sorted(bad))}"
    return True, ""

def parse_fasta(text: str) -> Dict[str, str]:
    """
    Parse FASTA or raw text into {header: cleaned_seq}.
    If no FASTA header ('>') is present, returns a single '(pasted sequence)' record.
    """
    text = (text or "").strip()
    if not text:
        return {}
    lines = text.splitlines()
    if not any(line.startswith(">") for line in lines):
        return {"(pasted sequence)": clean_seq(text)}

    records: Dict[str, str] = {}
    header, buf = None, []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = clean_seq("".join(buf))
            header = (line[1:].strip() or "unnamed")
            buf = []
        else:
            buf.append(line)
    if header is not None:
        records[header] = clean_seq("".join(buf))
    return records

def gc_content(seq: str) -> float:
    """GC percent over ACGT-only positions."""
    if not seq:
        return 0.0
    core = [b for b in seq if b in "ACGT"]
    return (sum(b in "GC" for b in core) / max(1, len(core))) * 100.0

def preview(seq: str, head: int = 60) -> str:
    """Compact preview string of a long sequence."""
    if len(seq) <= head:
        return seq
    return f"{seq[:head]}"

# microservices.py

import re

def get_sequence(seq_text: str) -> str:
    """
    Streamlit-friendly: accept either raw sequence text or FASTA text
    (already provided by the UI), return an uppercased, cleaned sequence.

    - No terminal prompts.
    - If multiple FASTA records are present, the first one is used.
    - If validate_seq() exists in this module, it's applied.
    """
    if not isinstance(seq_text, str):
        raise TypeError("get_sequence expects a string (raw or FASTA text).")

    s = seq_text.strip()
    if not s:
        return ""

    # If you already have parse_fasta() in this module (you do, per UI), use it:
    try:
        recs = parse_fasta(s)  # expects {header: sequence}
    except NameError:
        recs = {}

    if recs:
        # pick the first record
        _, seq = next(iter(recs.items()))
    else:
        # treat as raw sequence: strip headers/whitespace and non-letters
        # (allow '-' for gaps if needed)
        seq = re.sub(r"^>.*$", "", s, flags=re.MULTILINE)  # drop any FASTA headers if present
        seq = re.sub(r"\s+", "", seq)                      # remove whitespace
        seq = re.sub(r"[^A-Za-z\-*]", "", seq)             # keep letters, gaps, optional '*'

    # If you have a validator, use it
    try:
        ok, msg = validate_seq(seq)
        if not ok:
            raise ValueError(msg)
    except NameError:
        pass  # no validator available; skip

    return seq.upper()


def rotate_sequence(seq: str, start_subseq: str, *, allow_missing: bool = True) -> str:
    """
    Rotate a circular sequence so that it begins with start_subseq.
    - If start_subseq is not found and allow_missing=True, return the original seq (UI-friendly).
    - If allow_missing=False, raise ValueError (CLI/strict behavior).

    The sequence is treated as circular by doubling it.
    """
    if not seq or not start_subseq:
        return seq

    seq_len = len(seq)
    sub_len = len(start_subseq)
    if sub_len == 0 or sub_len > 2 * seq_len:
        return seq

    extended = seq + seq  # simulate circularity

    for i in range(seq_len):
        if extended[i:i + sub_len] == start_subseq:
            rotated = extended[i:i + seq_len]
            # Ensure exact prefix (handles repeats)
            if not rotated.startswith(start_subseq):
                idx = rotated.find(start_subseq)
                if idx != -1:
                    rotated = rotated[idx:] + rotated[:idx]
            return rotated

    if allow_missing:
        return seq
    raise ValueError("Starting subsequence not found in the sequence (even when treated as circular).")

# Suggest a reasonable number of chunks for the size of their genome

def default_n_sections(seq_len: int) -> int:
    if seq_len < 100:
        return 1
    elif seq_len < 100_000:
        return max(1, seq_len // 100)
    elif seq_len < 1_000_000:
        return 1000
    else:
        return max(1, seq_len // 1000)

# make a df with only anomolous chunks

def anomalous_rows(feature_df) -> "pd.DataFrame":
    return feature_df[feature_df["label"] != 1]


# ms for tab 2 of MAINAPP
def highlighted_block(a: str, b: str, a_label="REF", b_label="EVO", start_pos: int = 0) -> str:
    import html
    la = "".join(
        f"<span class='{'mm' if (i < len(b) and ch != b[i]) else 'ok'}'>{html.escape(ch)}</span>"
        for i, ch in enumerate(a)
    )
    lb = "".join(
        f"<span class='{'mm' if (i < len(a) and ch != a[i]) else 'ok'}'>{html.escape(ch)}</span>"
        for i, ch in enumerate(b)
    )
    lab_a = f"{a_label} {start_pos:>7d}: "
    lab_b = f"{b_label} {start_pos:>7d}: "
    return (
        "<div class='seqwrap'>"
        f"<div class='rowline'><span class='lab'>{lab_a}</span><span class='diffblock'>{la}</span></div>"
        f"<div class='rowline'><span class='lab'>{lab_b}</span><span class='diffblock'>{lb}</span></div>"
        "</div>"
    )

def diff_indices(code: int, start: int, end: int, info_matrix) -> list:
    r_row, e_row = (1, 3) if code in (1, 2, 3) else (10, 12)
    ref = info_matrix[r_row, start:end]
    evo = info_matrix[e_row, start:end]
    return [start + i for i, (a, b) in enumerate(zip(ref, evo)) if a != b]

def mut_type_from_frameshift(i: int, fs_row, revcomp: bool = False) -> str:
    try:
        seq_len = len(fs_row)
        if revcomp:
            i = seq_len - 1 - i
        curr = int(fs_row[i])
        prev = int(fs_row[i - 1]) if i > 0 else curr
        delta = curr - prev
    except Exception:
        return "point"
    if delta == 0:   return "point"
    if delta == -1:  return "insertion"
    if delta == 1:   return "deletion"
    return f"delta={delta}"

    # for active ORF tracker and prodigal implementation
def frame_active(value: int, frame: int) -> bool:
    """Check if a given frame (1,2,3) is active in the encoded value."""
    return str(frame) in str(value)




def run_prodigal(seq: str, seq_name: str = "seq") -> list[dict]:
    """
    Run Prodigal on a clean (no gaps) sequence string.
    Returns list of dicts with keys: start, stop, strand, frame
    Only returns forward strand ORFs — caller handles strand separately.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">{seq_name}\n{seq}\n")
        fasta_path = f.name

    gff_path = fasta_path + ".gff"

    try:
        subprocess.run(
            ["prodigal", "-i", fasta_path, "-f", "gff", "-o", gff_path,
             "-p", "meta", "-q"],
            check=True,
            capture_output=True,
        )
        orfs = []
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                start  = int(parts[3]) - 1  # GFF is 1-based → 0-based
                stop   = int(parts[4])
                strand = parts[6]
                frame  = (start % 3) + 1    # 1, 2, or 3
                orfs.append({
                    "start":  start,
                    "stop":   stop,
                    "strand": strand,
                    "frame":  frame,
                })
        return orfs
    finally:
        os.unlink(fasta_path)
        if os.path.exists(gff_path):
            os.unlink(gff_path)


def build_active_orf_tracker(
    seq: str,
    pos_map: list[int],
    seq_len: int,
    orfs: list[dict],
) -> list[int]:
    """
    Build an active ORF tracker row of length seq_len (original coordinate space).

    Encoding at each position:
        0         = no active ORF
        1/2/3     = one active ORF in that reading frame
        12/23/13  = two simultaneously active ORFs
        112/123   = nested ORFs (digits concatenated, sorted)

    A stop codon on frame X wipes ALL active ORFs on frame X.
    Values are stored at original (gap-inclusive) coordinates via pos_map.
    Only forward strand ORFs should be passed — caller handles strand.
    """
    # work in clean-seq coordinate space first
    clean_len = len(pos_map)
    # per-frame open ORF counter in clean-seq space
    frame_counts = {1: 0, 2: 0, 3: 0}

    # build a list of events sorted by position
    # event: (clean_pos, type, frame)
    # type: 'start' or 'stop'
    events: list[tuple[int, str, int]] = []
    for orf in orfs:
        if orf["strand"] != "+":
            continue
        frame = orf["frame"]
        events.append((orf["start"], "start", frame))
        events.append((orf["stop"], "stop",  frame))  

    events.sort(key=lambda e: e[0])

    # build event map: clean_pos → list of events at that position
    from collections import defaultdict
    event_map: dict[int, list[tuple[str, int]]] = defaultdict(list)
    for pos, etype, frame in events:
        event_map[pos].append((etype, frame))

    # walk clean sequence, track active ORFs, write to original coords
    tracker = [0] * seq_len

    for clean_idx in range(clean_len):
        # apply any events at this position
        if clean_idx in event_map:
            for etype, frame in event_map[clean_idx]:
                if etype == "start":
                    frame_counts[frame] += 1
                elif etype == "stop":
                    frame_counts[frame] = 0  # wipe ALL on this frame

        # encode current state
        digits = []
        for frame in (1, 2, 3):
            digits.extend([frame] * frame_counts[frame])
        value = int("".join(str(d) for d in digits)) if digits else 0

        # write to original coordinate
        original_pos = pos_map[clean_idx]
        tracker[original_pos] = value

    # gap positions stay 0 (already initialized)
# mark gap positions with -1 so boundary detection isn't fooled by gaps

    tracker_arr = np.array(tracker)
    covered = np.zeros(seq_len, dtype=bool)
    covered[pos_map] = True
    tracker_arr[~covered] = -1
    return tracker_arr.tolist()

from collections import defaultdict

def _parse_active_value(value: int) -> list[int]:
    """
    Parse an active ORF tracker integer into a list of active frame numbers.
    Example: 113 -> [1, 1, 3], 23 -> [2, 3], 0 -> []
    """
    if value <= 0:
        return []
    return [int(d) for d in str(value)]


def _find_orf_boundaries(
    tracker_row: np.ndarray,
    mutation_idx: int,
    frame: int,
    seq_len: int,
) -> tuple[int, int] | None:
    """
    Given a mutation index and a specific frame number, search backwards and
    forwards in tracker_row to find where that frame's ORF begins and ends.
    Skips -1 (gap sentinel) values.
    Returns (start, end) in global info_matrix coordinates, or None if not found.
    """
    # search backwards for start
    start = None
    for i in range(mutation_idx, -1, -1):
        val = int(tracker_row[i])
        if val == -1:
            continue  # skip gap sentinels
        frames_here = _parse_active_value(val)
        if frame not in frames_here:
            start = i + 1
            break
    if start is None:
        start = 0  # hit the beginning of the sequence

    # search forwards for end
    end = None
    for i in range(mutation_idx, seq_len):
        val = int(tracker_row[i])
        if val == -1:
            continue  # skip gap sentinels
        frames_here = _parse_active_value(val)
        if frame not in frames_here:
            end = i
            break
    if end is None:
        end = seq_len  # hit the end of the sequence

    if start >= end:
        return None
    return (start, end)


def _extract_affected_orfs(
    tracker_row: np.ndarray,
    mutation_idx: int,
    seq_len: int,
) -> list[tuple[int, int, int]]:
    val = int(tracker_row[mutation_idx])
    if val <= 0:  # 0 = no ORF, -1 = gap — skip entirely
        return []
    # only search for frames present at THIS mutation index
    # do not drift into neighboring ORFs
    frames_at_mutation = _parse_active_value(val)
    seen = set()
    results = []
    for frame in frames_at_mutation:
        boundaries = _find_orf_boundaries(tracker_row, mutation_idx, frame, seq_len)
        if boundaries is None:
            continue
        # verify the mutation index is actually within the boundaries
        if not (boundaries[0] <= mutation_idx <= boundaries[1]):
            continue  # boundary drifted — skip
        key = (frame, boundaries[0], boundaries[1])
        if key not in seen:
            seen.add(key)
            results.append(key)
    return results

def diff_indices_span(code: int, start: int, end: int, info_matrix) -> list:
    r_row, e_row = (1, 3) if code in (1, 2, 3) else (10, 12)
    if start < 0 or end < 0:
        return []
    ref = info_matrix[r_row, start:end]
    evo = info_matrix[e_row, start:end]
    return [start + i for i, (a, b) in enumerate(zip(ref, evo))
            if a != b and a not in ('-', 'N') and b not in ('-', 'N')]

#new pairing system helpers
def num_digits_vec(arr):
    """Vectorized digit counter for ORF tracker arrays. -1 stays -1."""
    result = np.zeros(len(arr), dtype=int)
    mask = arr > 0
    if mask.any():
        result[mask] = np.array([len(str(int(v))) for v in arr[mask]])
    result[arr == -1] = -1
    return result

def get_transition_sequence(window):
    stripped = [(i, v) for i, v in enumerate(window) if v != -1]
    transitions = []
    for pos in range(1, len(stripped)):
        _, curr = stripped[pos]
        _, prev = stripped[pos - 1]
        if curr != prev:
            transitions.append((pos, num_digits(prev), num_digits(curr)))
    return transitions

def truncate_at_neg1(window, start_from=2):
    for i in range(start_from, len(window)):
        if window[i] == -1:
            return window[:i]
    return window

def num_digits(val):
    if val == 0:
        return 1
    return len(str(int(abs(val))))

def should_remove_pair(win_a, win_b):
    wa, wb = list(win_a), list(win_b)
    a1 = wa[1] if len(wa) > 1 else 0
    b1 = wb[1] if len(wb) > 1 else 0
    if a1 > 0 or b1 > 0:
        return False
    wa = truncate_at_neg1(wa, start_from=2)
    wb = truncate_at_neg1(wb, start_from=2)
    trans_a = get_transition_sequence(wa)
    trans_b = get_transition_sequence(wb)
    return trans_a == trans_b

# ── AA row lookup ────────────────────────────────────────────────────────────
_AA_ROW_MAP = {
    (8,  1): 14, (8,  2): 15, (8,  3): 16,
    (9,  1): 17, (9,  2): 18, (9,  3): 19,
    (11, 1): 20, (11, 2): 21, (11, 3): 22,
    (13, 1): 23, (13, 2): 24, (13, 3): 25,
}

_REF_TRACKER_ROWS  = {8, 11}
_EVO_TRACKER_ROWS  = {9, 13}
_FWD_TRACKER_ROWS  = {8, 9}
_REV_TRACKER_ROWS  = {11, 13}

def aa_row_for(tracker_row: int, frame: int) -> int:
    return _AA_ROW_MAP[(tracker_row, frame)]

def is_ref_row(tracker_row: int) -> bool:
    return tracker_row in _REF_TRACKER_ROWS

def is_fwd_row(tracker_row: int) -> bool:
    return tracker_row in _FWD_TRACKER_ROWS

def score_orf_pairing(
    info_matrix: np.ndarray,
    ref_row: int, ref_frame: int, ref_start: int, ref_end: int,
    evo_row: int, evo_frame: int, evo_start: int, evo_end: int,
) -> dict:
    """
    Score a candidate ref/evo ORF pairing by counting identical AAs
    in the overlapping coordinate window.
    Returns dict with overlap_start, overlap_end, match_count, overlap_len, score.
    """
    overlap_start = max(ref_start, evo_start)
    overlap_end   = min(ref_end,   evo_end)

    if overlap_start >= overlap_end:
        return {
            "overlap_start": overlap_start,
            "overlap_end":   overlap_end,
            "match_count":   0,
            "overlap_len":   0,
            "score":         0.0,
        }

    ref_aa_row = aa_row_for(ref_row, ref_frame)
    evo_aa_row = aa_row_for(evo_row, evo_frame)

    _skip = {'X', '*', '-', None, 'N', ''}

    match_count = 0
    overlap_len = 0
    for i in range(overlap_start, overlap_end):
        r = info_matrix[ref_aa_row, i]
        e = info_matrix[evo_aa_row, i]
        if r in _skip or e in _skip:
            continue
        overlap_len += 1
        if r == e:
            match_count += 1

    score = match_count / overlap_len if overlap_len > 0 else 0.0
    return {
        "overlap_start": overlap_start,
        "overlap_end":   overlap_end,
        "match_count":   match_count,
        "overlap_len":   overlap_len,
        "score":         score,
    }


def compute_pairing_scores(
    info_matrix: np.ndarray,
    mut_idx: int,
    orfs_at_index: list[dict],
) -> list[dict]:
    """
    For a mutation index with >2 ORFs, compute scores for all valid
    ref/evo pairings within the same direction.
    Returns list of pairing dicts sorted by score descending.
    """
    ref_orfs = [o for o in orfs_at_index if is_ref_row(o["row"])]
    evo_orfs = [o for o in orfs_at_index if not is_ref_row(o["row"])]

    pairings = []
    for r in ref_orfs:
        for e in evo_orfs:
            # must be same direction
            if is_fwd_row(r["row"]) != is_fwd_row(e["row"]):
                continue
            s = score_orf_pairing(
                info_matrix,
                r["row"], r["orf"], r["start"], r["end"],
                e["row"], e["orf"], e["start"], e["end"],
            )
            pairings.append({
                "mut_idx":    mut_idx,
                "ref_row":    r["row"],
                "ref_frame":  r["orf"],
                "ref_start":  r["start"],
                "ref_end":    r["end"],
                "evo_row":    e["row"],
                "evo_frame":  e["orf"],
                "evo_start":  e["start"],
                "evo_end":    e["end"],
                **s,
            })

    pairings.sort(key=lambda x: x["score"], reverse=True)
    return pairings


def build_seq_diffs(info_matrix: np.ndarray) -> list:
    """Build and filter seq_diffs records (forward + revcomp)."""
    from collections import defaultdict

    def _build_records(ref_nuc_row, evo_nuc_row, trk_rows):
        mutation_sites = np.where(
            info_matrix[ref_nuc_row, :].astype(object) !=
            info_matrix[evo_nuc_row, :].astype(object)
        )[0]
        records = []
        for idx in mutation_sites:
            for row in trk_rows:
                end = min(int(idx) + 21, info_matrix.shape[1])
                window = info_matrix[row, max(0, int(idx)-1):end]
                records.append((idx, row, window))
        return records

    records_f = _build_records(1,  3,  [8, 9])
    records_r = _build_records(10, 12, [11, 13])
    combined  = records_f + records_r

    paired = defaultdict(list)
    for idx, row, window in combined:
        paired[int(idx)].append((idx, row, window))

    sorted_indices = sorted(paired.keys())

    def split_consecutive(idxs):
        if not idxs:
            return []
        runs, current = [], [idxs[0]]
        for i in range(1, len(idxs)):
            if idxs[i] - idxs[i-1] == 1:
                current.append(idxs[i])
            else:
                runs.append(current)
                current = [idxs[i]]
        runs.append(current)
        return runs

    def pair_w1_sig(pair):
        sig = {}
        for idx, row, window in pair:
            sig[row] = int(window[1]) if len(window) > 1 else 0
        return tuple(sorted(sig.items()))

    # removal filter
    kept_after_removal = set()
    for idx in sorted_indices:
        pair = paired[idx]
        if len(pair) != 2:
            kept_after_removal.add(idx)
            continue
        (_, _, win_a), (_, _, win_b) = pair
        if not should_remove_pair(list(win_a), list(win_b)):
            kept_after_removal.add(idx)

    sorted_kept = sorted(kept_after_removal)
    if not sorted_kept:
        return []

    # consecutive mutation filter
    runs = split_consecutive(sorted_kept)
    final_kept = set()
    for run in runs:
        if len(run) == 1:
            final_kept.add(run[0])
            continue
        sub_run_start = 0
        current_sig = pair_w1_sig(paired[run[0]])
        for i in range(1, len(run)):
            sig = pair_w1_sig(paired[run[i]])
            if sig != current_sig:
                final_kept.add(run[sub_run_start])
                final_kept.add(run[i - 1])
                sub_run_start = i
                current_sig = sig
        final_kept.add(run[sub_run_start])
        final_kept.add(run[-1])

    result = []
    for idx in sorted(final_kept):
        result.extend(paired[idx])
    return result


def build_new_result_matrix(
    seq_diffs: list,
    info_matrix: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """
    Build new_result_matrix with columns [mut_idx, row, orf, start, end].
    Also returns pairing_scores dict keyed by mut_idx for ambiguous sites.
    """
    seq_len = info_matrix.shape[1]
    rows = []

    # group by mut_idx to detect ambiguous sites
    from collections import defaultdict
    by_idx = defaultdict(list)

    for idx, trk_row, window in seq_diffs:
        mut_idx = int(idx)
        tracker = info_matrix[trk_row, :].astype(int)
        search_idx = mut_idx
        while search_idx < seq_len and tracker[search_idx] == -1:
            search_idx += 1
        bounds = _extract_affected_orfs(tracker, search_idx, seq_len)
        if not bounds:
            rows.append([mut_idx, trk_row, -1, -1, -1])
        else:
            for frame, start, end in bounds:
                rows.append([mut_idx, trk_row, frame, start, end])
                by_idx[mut_idx].append({
                    "row":   trk_row,
                    "orf":   frame,
                    "start": start,
                    "end":   end,
                })

    # compute pairing scores for ambiguous indices
    pairing_scores = {}
    for mut_idx, orfs in by_idx.items():
        ref_orfs = [o for o in orfs if is_ref_row(o["row"])]
        evo_orfs = [o for o in orfs if not is_ref_row(o["row"])]
        # ambiguous = more than one valid pairing exists
        valid_pairs = [
            (r, e) for r in ref_orfs for e in evo_orfs
            if is_fwd_row(r["row"]) == is_fwd_row(e["row"])
        ]
        if len(valid_pairs) > 1:
            pairing_scores[mut_idx] = compute_pairing_scores(
                info_matrix, mut_idx, orfs
            )

    mat = np.array(rows, dtype=object) if rows else np.empty((0, 5), dtype=object)
    return mat, pairing_scores