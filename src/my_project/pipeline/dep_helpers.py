import streamlit as st
import re
import io

def _clean_seq(seq: str) -> str:
    """Uppercase, remove non-letters/dashes, convert U->T."""
    s = re.sub(r"[^A-Za-z-]", "", seq).upper()
    return s.replace("U", "T")

def _validate_seq(seq: str) -> tuple[bool, str]:
    bad = set(ch for ch in seq if ch not in AMBIGUOUS_DNA)
    if bad:
        return False, f"Found invalid characters: {''.join(sorted(bad))}"
    if len(seq) == 0:
        return False, "Sequence is empty after cleaning."
    return True, ""

def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    bases = [b for b in seq if b in "ACGT"]
    return (sum(b in "GC" for b in bases) / max(1, len(bases))) * 100.0

def _preview(seq: str, head=60, tail=60) -> str:
    if len(seq) <= head + tail:
        return seq
    return f"{seq[:head]} … {seq[-tail:]}"

def _parse_fasta(text: str) -> dict:
    """
    Return dict {header: sequence}. If no '>' is found, treat as single unnamed sequence.
    """
    text = text.strip()
    if not text:
        return {}

    if not any(line.startswith(">") for line in text.splitlines()):
        cleaned = _clean_seq(text)
        return {"(pasted sequence)": cleaned}

    records = {}
    header = None
    seq_lines = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = _clean_seq("".join(seq_lines))
            header = line[1:].strip() or "unnamed"
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        records[header] = _clean_seq("".join(seq_lines))
    return records

@st.cache_data(show_spinner=False)
def _read_uploaded_text(file) -> str:
    return file.getvalue().decode("utf-8", errors="replace") if file else ""

def _picker(label: str):
    """
    UI block that returns (chosen_header, chosen_seq) for one input (ref/evo),
    reading from either an upload or a paste.
    """
    st.subheader(label)
    t_upload, t_paste = st.tabs(["Upload FASTA", "Paste FASTA or raw sequence"])

    chosen_header = None
    chosen_seq = ""

def _picker(label: str):
    """
    UI block that returns (chosen_header, chosen_seq) for one input (ref/evo),
    reading from either an upload or a paste.
    """
    st.subheader(label)
    t_upload, t_paste = st.tabs(["Upload FASTA", "Paste FASTA or raw sequence"])

    chosen_header = None
    chosen_seq = ""
