# UI.py
import io
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
import re
import io
from helpers import _picker
from helpers import _clean_seq
from helpers import _gc_content
from helpers import _preview
from helpers import _presence


st.set_page_config(page_title="Upload Sequences", page_icon="🧬", layout="wide")
st.title("🧬 Upload / Paste Sequences")
st.set_page_config(page_title="Genome Comperator", layout="wide")

# ---------------- Sidebar ----------------
#with st.sidebar:
#    st.header("Controls")
#    name = st.text_input("Your name", value="Guest")
#    slider_n = st.slider("Rows for demo DataFrame", 10, 1000, 200, step=10)
#    noise = st.slider("Noise level", 0.0, 1.0, 0.2, step=0.05)
#    st.caption("Change values and watch the app rerun.")

# ---------------- Session State ----------------
#if "clicks" not in st.session_state:
#    st.session_state.clicks = 0

#def bump():
#    st.session_state.clicks += 1


st.title("Welcome to Genome Comperator :)")
#st.write(f"Hello, **{name}**.")

#st.button("Increment counter", on_click=bump)
#st.write(f"Counter value: `{st.session_state.clicks}`")

#st.divider()

# ---------------- Tabs ----------------
tab1, tab2, tab3, tab4 = st.tabs(["Upload sequences", "Data & Charts", "Files", "State & Caching"])

# ---- Tab 1: Widgets + Form ----
with tab1:
    ref_h, ref_seq = _picker("Reference sequence")
    st.divider()
    evo_h, evo_seq = _picker("Evolved sequence")
    st.divider()



    start_h, start_seq = _picker("Starting subsequence")
    st.divider()


    # ---------- Validation & Save ----------



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
    if st.button("Save to session", type="primary", width="stretch"):
        if not (ref_seq and evo_seq):
            st.error("Please provide both Reference and Evolved sequences.")
        else:
            st.session_state["ref_header"] = ref_h or "reference"
            st.session_state["ref_seq"] = ref_seq
            st.session_state["evo_header"] = evo_h or "evolved"
            st.session_state["evo_seq"] = evo_seq
            st.session_state["start_header"] = start_h or "starting_subsequence"
            st.session_state["start_subseq"] = start_seq
            st.success("Saved to session state.")

    # Optional: show a compact summary
    with st.expander("Session summary"):
        for k in ("ref_header", "evo_header", "start_subseq"):
            st.write(f"- {k}: {st.session_state.get(k, '(none)')}")
        for name, seq in (("ref_seq", st.session_state.get("ref_seq")), ("evo_seq", st.session_state.get("evo_seq"))):
            if seq:
                st.write(f"- {name}: {len(seq):,} bp, GC {_gc_content(seq):.1f}%")
                st.code(_preview(seq), language="text")
