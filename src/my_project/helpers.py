import streamlit as st
import re
import io

AMBIGUOUS_DNA = set("ACGTURYSWKMBDHVN-")

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
#    t_upload, t_paste = st.tabs(["Upload FASTA", "Paste FASTA or raw sequence"])

    chosen_header = None
    chosen_seq = ""

    with t_upload:
        file = st.file_uploader(f"{label} file", type=["fa","fasta","fna","ffn","faa","txt"], key=f"{label}_file")
        if file:
            text = _read_uploaded_text(file)
            recs = _parse_fasta(text)
            if not recs:
                st.error("Could not read any sequence from this file.")
            else:
                headers = list(recs.keys())
                if len(headers) > 1:
                    h = st.selectbox("Choose record", headers, key=f"{label}_sel")
                else:
                    h = headers[0]
                seq = recs[h]
                ok, msg = _validate_seq(seq)
                if not ok:
                    st.error(msg)
                else:
                    chosen_header, chosen_seq = h, seq
                    col1, col2, col3 = st.columns([1,1,2])
                    col1.metric("Length", f"{len(seq):,}")
                    col2.metric("GC %", f"{_gc_content(seq):.1f}")
                    with col3:
                        st.caption("Preview")
                        st.code(_preview(seq), language="text")

    with t_paste:
        pasted = st.text_area(
            f"{label} (paste FASTA or raw)",
            height=180,
            placeholder=">my_sequence\nACGTACGT...\n(or paste raw sequence without header)",
            key=f"{label}_paste"
        )
        if pasted.strip():
            recs = _parse_fasta(pasted)
            if not recs:
                st.error("Nothing parsed from pasted text.")
            else:
                headers = list(recs.keys())
                if len(headers) > 1:
                    h2 = st.selectbox("Choose record", headers, key=f"{label}_paste_sel")
                else:
                    h2 = headers[0]
                seq2 = recs[h2]
                ok2, msg2 = _validate_seq(seq2)
                if not ok2:
                    st.error(msg2)
                else:
                    chosen_header, chosen_seq = h2, seq2
                    col1, col2, col3 = st.columns([1,1,2])
                    col1.metric("Length", f"{len(seq2):,}")
                    col2.metric("GC %", f"{_gc_content(seq2):.1f}")
                    with col3:
                       st.caption("Preview")
                       st.code(_preview(seq2), language="text")

    return chosen_header, chosen_seq

ref_h, ref_seq = _picker("Reference sequence")
#st.divider()
evo_h, evo_seq = _picker("Evolved sequence")
#st.divider()

# --- Starting subsequence picker (same upload/paste UI) ---
start_h, start_seq = _picker("Starting subsequence")
st.divider()

# ---------- Validation & Save ----------
def _presence(name, seq):
    return f"✅ {name}: {len(seq):,} bp loaded" if seq else f"⚠️ {name} not provided"

st.markdown(
    f"{_presence('Reference', ref_seq)}  \n"
    f"{_presence('Evolved', evo_seq)}  \n"
    f"{_presence('Starting subsequence', start_seq)}"
)

problems = []
if ref_seq and start_seq and start_seq not in ref_seq:
    problems.append("Starting subsequence not found in **Reference**.")
if evo_seq and start_seq and start_seq not in evo_seq:
    problems.append("Starting subsequence not found in **Evolved**.")
for p in problems:
    st.warning(p)

# Save to session (note: no widget uses key 'start_subseq', so this is safe)
#if st.button("Save to session", type="primary", width="stretch"):
#    if not (ref_seq and evo_seq):
#        st.error("Please provide both Reference and Evolved sequences.")
#    else:
#        st.session_state["ref_header"] = ref_h or "reference"
#        st.session_state["ref_seq"] = ref_seq
#        st.session_state["evo_header"] = evo_h or "evolved"
#        st.session_state["evo_seq"] = evo_seq
#        st.session_state["start_header"] = start_h or "starting_subsequence"
#        st.session_state["start_subseq"] = start_seq
#        st.success("Saved to session state.")

def _presence(name, seq):
    return f"✅ {name}: {len(seq):,} bp loaded" if seq else f"⚠️ {name} not provided"