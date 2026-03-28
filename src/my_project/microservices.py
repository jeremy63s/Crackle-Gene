# microservices.py
from __future__ import annotations
import re
from typing import Dict, Tuple

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

def mut_type_from_frameshift(i: int, fs_row) -> str:
    try:
        curr = int(fs_row[i])
        prev = int(fs_row[i - 1]) if i > 0 else curr
        delta = curr - prev
    except Exception:
        return "point"
    if delta == 0:   return "point"
    if delta == -1:  return "insertion"
    if delta == 1:   return "deletion"
    return f"delta={delta}"