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
        events.append((orf["stop"] - 3, "stop",  frame))  # stop codon starts 3 before end

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

def _frameshift_at(info_matrix: np.ndarray, pos: int, forward: bool) -> int:
    """
    Returns (ref_nt_count - evo_nt_count) % 3 for non-spacer nucleotides
    from position 0 up to but not including pos.
    For forward sequences uses rows 1 and 3, for revcomp uses rows 10 and 12.
    """
    if forward:
        ref_row, evo_row = 1, 3
    else:
        ref_row, evo_row = 10, 12
    ref_count = sum(1 for v in info_matrix[ref_row, :pos] if v not in ('-', 'N', 'X'))
    evo_count = sum(1 for v in info_matrix[evo_row, :pos] if v not in ('-', 'N', 'X'))
    return (ref_count - evo_count) % 3

def find_mutated_orfs(
    info_matrix: np.ndarray,
    forward: bool = True,
) -> list[dict]:
    if forward:
        ref_nuc_row, evo_nuc_row = 1, 3
        ref_trk_row, evo_trk_row = 8, 9
    else:
        ref_nuc_row, evo_nuc_row = 10, 12
        ref_trk_row, evo_trk_row = 11, 13

    seq_len = info_matrix.shape[1]
    ref_nuc = info_matrix[ref_nuc_row]
    evo_nuc = info_matrix[evo_nuc_row]
    ref_trk = info_matrix[ref_trk_row].astype(int)
    evo_trk = info_matrix[evo_trk_row].astype(int)

    mismatch_indices = [
        i for i in range(seq_len)
        if ref_nuc[i] != evo_nuc[i]
    ]

    ref_orfs_at_mutations: dict[tuple[int,int,int], list[int]] = defaultdict(list)
    evo_orfs_at_mutations: dict[tuple[int,int,int], list[int]] = defaultdict(list)

    for idx in mismatch_indices:
        for orf in _extract_affected_orfs(ref_trk, idx, seq_len):
            ref_orfs_at_mutations[orf].append(idx)
        for orf in _extract_affected_orfs(evo_trk, idx, seq_len):
            evo_orfs_at_mutations[orf].append(idx)

    ref_orf_set = set(ref_orfs_at_mutations.keys())
    evo_orf_set = set(evo_orfs_at_mutations.keys())

    results = []
    paired_evo = set()

    # tier 1: same absolute start AND end coordinates regardless of frame
    for ref_orf in sorted(ref_orf_set):
        ref_frame, ref_start, ref_end = ref_orf
        for evo_orf in sorted(evo_orf_set):
            evo_frame, evo_start, evo_end = evo_orf
            if ref_start == evo_start and ref_end == evo_end:
                results.append({
                    "ref_orf": ref_orf,
                    "evo_orf": evo_orf,
                    "pair_tier": 1,
                    "forward": forward,
                    "mutation_indices": sorted(set(
                        ref_orfs_at_mutations[ref_orf] +
                        evo_orfs_at_mutations[evo_orf]
                    )),
                })
                paired_evo.add(evo_orf)
                break

    paired_ref = {r["ref_orf"] for r in results}
    unpaired_ref = ref_orf_set - paired_ref
    unpaired_evo = evo_orf_set - paired_evo

    # tier 2: start OR stop aligns, and relative frameshift is consistent
    for ref_orf in sorted(unpaired_ref):
        ref_frame, ref_start, ref_end = ref_orf
        best_evo = None
        best_score = -1

        for evo_orf in sorted(unpaired_evo):
            evo_frame, evo_start, evo_end = evo_orf

            start_matches = ref_start == evo_start
            end_matches = ref_end == evo_end

            if not start_matches and not end_matches:
                continue

            # check relative frameshift at the matching codon boundary
            check_pos = ref_start if start_matches else ref_end
            fs = _frameshift_at(info_matrix, check_pos, forward)
            expected_evo_frame = ((ref_frame - 1 - fs) % 3) + 1
            if evo_frame != expected_evo_frame:
                continue

            # prefer end match over start match as tiebreaker
            score = 2 if end_matches else 1
            if score > best_score:
                best_score = score
                best_evo = evo_orf

        if best_evo is not None:
            results.append({
                "ref_orf": ref_orf,
                "evo_orf": best_evo,
                "pair_tier": 2,
                "forward": forward,
                "mutation_indices": sorted(set(
                    ref_orfs_at_mutations[ref_orf] +
                    evo_orfs_at_mutations[best_evo]
                )),
            })
            paired_evo.add(best_evo)
            unpaired_evo.discard(best_evo)

    paired_ref = {r["ref_orf"] for r in results}
    unpaired_ref = ref_orf_set - paired_ref

    # tier 3: one ORF starts inside the other's span
    for ref_orf in sorted(unpaired_ref):
        ref_frame, ref_start, ref_end = ref_orf
        best_evo = None
        best_overlap = 0

        for evo_orf in sorted(unpaired_evo):
            evo_frame, evo_start, evo_end = evo_orf
            # check if either starts inside the other
            ref_starts_inside_evo = evo_start <= ref_start <= evo_end
            evo_starts_inside_ref = ref_start <= evo_start <= ref_end
            if not ref_starts_inside_evo and not evo_starts_inside_ref:
                continue
            overlap = min(ref_end, evo_end) - max(ref_start, evo_start)
            if overlap > best_overlap:
                best_overlap = overlap
                best_evo = evo_orf

        if best_evo is not None:
            results.append({
                "ref_orf": ref_orf,
                "evo_orf": best_evo,
                "pair_tier": 3,
                "forward": forward,
                "mutation_indices": sorted(set(
                    ref_orfs_at_mutations[ref_orf] +
                    evo_orfs_at_mutations[best_evo]
                )),
            })
            paired_evo.add(best_evo)
            unpaired_evo.discard(best_evo)

    paired_ref = {r["ref_orf"] for r in results}
    unpaired_ref = ref_orf_set - paired_ref

    # remaining orphans
    for ref_orf in sorted(unpaired_ref):
        results.append({
            "ref_orf": ref_orf,
            "evo_orf": None,
            "pair_tier": 3,
            "forward": forward,
            "mutation_indices": ref_orfs_at_mutations[ref_orf],
        })
    for evo_orf in sorted(unpaired_evo):
        results.append({
            "ref_orf": None,
            "evo_orf": evo_orf,
            "pair_tier": 3,
            "forward": forward,
            "mutation_indices": evo_orfs_at_mutations[evo_orf],
        })

    return results

def diff_indices_span(code: int, start: int, end: int, info_matrix) -> list:
    r_row, e_row = (1, 3) if code in (1, 2, 3) else (10, 12)
    if start < 0 or end < 0:
        return []
    ref = info_matrix[r_row, start:end]
    evo = info_matrix[e_row, start:end]
    return [start + i for i, (a, b) in enumerate(zip(ref, evo))
            if a != b and a not in ('-', 'N') and b not in ('-', 'N')]