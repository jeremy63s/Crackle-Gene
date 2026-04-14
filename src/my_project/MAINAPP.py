from __future__ import annotations
import os
import html
import json
import re
import threading, queue, time

os.environ.setdefault("MPLBACKEND", "Agg")  # force headless before pyplot

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import parasail
import string
import streamlit as st
from typing import Tuple, List, Optional
from my_project.pypeline import run_pipeline
import streamlit.components.v1 as components
import microservices as ms  # pure helpers (no Streamlit inside)
from my_project.blast import process_result_matrix


INFO_MATRIX_CHEATSHEET = [
   (0,  "Alignment data (chunkwise)"),
   (1,  "Reference nucleotide sequence"),
   (2,  "Reference nucleotide index"),
   (3,  "Evolved nucleotide sequence"),
   (4,  "Evolved nucleotide index"),
   (5,  "Frameshift counter (del:+1, ins:-1)"),
   (6,  "Mutation type (in/del:0, point:-1)"),
   (7,  "Global index"),
   (8,  "Reference active ORF tracker"),
   (9,  "Evolved active ORF tracker"),
   (10, "RevComp of reference seq"),
   (11, "RevComp reference active ORF tracker"),
   (12, "RevComp of evolved seq"),
   (13, "RevComp evolved active ORF tracker"),
   (14, "Reference ORF frame 1"),
   (15, "Reference ORF frame 2"),
   (16, "Reference ORF frame 3"),
   (17, "Evolved ORF frame 1"),
   (18, "Evolved ORF frame 2"),
   (19, "Evolved ORF frame 3"),
   (20, "RevComp reference ORF frame 1"),
   (21, "RevComp reference ORF frame 2"),
   (22, "RevComp reference ORF frame 3"),
   (23, "RevComp evolved ORF frame 1"),
   (24, "RevComp evolved ORF frame 2"),
   (25, "RevComp evolved ORF frame 3"),
]


#st.set_page_config(page_title="Crackle Gene", page_icon="🧬", layout="wide")
st.title("Welcome to Crackle Gene 🧬 :)")


@st.cache_data(show_spinner=False)
def _decode_uploaded(uploaded) -> str:
   return uploaded.getvalue().decode("utf-8", errors="replace") if uploaded else ""


def sequence_picker(label: str, key_prefix: str, session_key: str):
   """
   Render one input with two tabs (Upload FASTA / Paste RAW or FASTA).
   Saves the final, cleaned sequence to st.session_state[session_key].
   key_prefix must be unique per picker: e.g. 'ref', 'evo', 'start'.
   Returns (header, sequence).
   """
   st.subheader(label)
   t_upload, t_paste = st.tabs(["Upload FASTA", "Paste FASTA or raw sequence"])


   chosen_header, chosen_seq = None, ""


   # --- Upload tab ---
   with t_upload:
       f = st.file_uploader(
           f"{label} file",
           type=["fa", "fasta", "fna", "ffn", "faa", "txt"],
           key=f"{key_prefix}_file",
       )
       if f:
           text = _decode_uploaded(f)
           recs = ms.parse_fasta(text)
           if not recs:
               st.error("Could not read any sequence from this file.")
           else:
               headers = list(recs.keys())
               h = st.selectbox("Choose record", headers, key=f"{key_prefix}_sel") if len(headers) > 1 else headers[0]
               seq = recs[h]
               ok, msg = ms.validate_seq(seq)
               if not ok:
                   st.error(msg)
               else:
                   chosen_header, chosen_seq = h, seq
                   st.session_state[session_key] = seq  # auto-save
                   c1, c2, c3 = st.columns([1, 1, 2])
                   c1.metric("Length", f"{len(seq):,}")
                   c2.metric("GC %", f"{ms.gc_content(seq):.1f}")
                   with c3:
                       st.caption("Preview")
                       st.code(ms.preview(seq), language="text")


   # --- Paste tab ---
   with t_paste:
       pasted = st.text_area(
           f"{label} (paste FASTA or raw)",
           height=160,
           placeholder=">my_sequence\nACGT...\n(or paste raw sequence without header)",
           key=f"{key_prefix}_paste",
       )
       if pasted.strip():
           recs = ms.parse_fasta(pasted)
           if not recs:
               st.error("Nothing parsed from pasted text.")
           else:
               headers = list(recs.keys())
               h2 = st.selectbox("Choose record", headers, key=f"{key_prefix}_paste_sel") if len(headers) > 1 else headers[0]
               seq2 = recs[h2]
               ok2, msg2 = ms.validate_seq(seq2)
               if not ok2:
                   st.error(msg2)
               else:
                   chosen_header, chosen_seq = h2, seq2
                   st.session_state[session_key] = seq2  # auto-save
                   c1, c2, c3 = st.columns([1, 1, 2])
                   c1.metric("Length", f"{len(seq2):,}")
                   c2.metric("GC %", f"{ms.gc_content(seq2):.1f}")
                   with c3:
                       st.caption("Preview")
                       st.code(ms.preview(seq2), language="text")


   st.session_state.setdefault(session_key, "")
   return chosen_header, st.session_state[session_key]


# ---------- Tabs ----------
page = st.sidebar.radio(
    "Navigation",
    ["Upload sequences", "Submit to blast", "View", "BLAST results", "Debug"],
    key="nav_page",
)


if page == "Upload sequences":
    st.header("Upload / Paste Sequences")
    st.caption(f"BLAST flag = {bool(st.session_state.get('do_blast', False))}")
    st.caption(f"run_pipeline from {run_pipeline.__module__}")

    ref_h, ref_seq = sequence_picker("Reference sequence", key_prefix="ref", session_key="ref_seq")
    st.divider()
    evo_h, evolved_seq = sequence_picker("Query sequence", key_prefix="evo", session_key="evolved_seq")
    st.divider()
    start_h, start_subseq = sequence_picker("Starting subsequence", key_prefix="start", session_key="start_subseq")

    st.divider()
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        _ref_seq = st.session_state.get("ref_seq") or ""
        _default_n = ms.default_n_sections(len(_ref_seq)) if _ref_seq else 1
        n_sections = st.number_input(
            "Number of chunks (n)",
            min_value=1,
            value=_default_n,
            step=1,
            help=f"Auto-set from reference length ({len(_ref_seq):,} bp). Adjust manually if needed.",
        )
    with c2:
        use_isolation_forest = st.toggle(
            "Use Isolation Forest (auto)",
            value=True,
            help="Default: Isolation Forest detects mutated chunks automatically from multi-feature analysis. "
                 "Turn off to use the manual normalized score threshold instead.",
        )
        threshold = st.number_input(
            "Normalized score threshold",
            min_value=0.0, max_value=1.0, value=0.995, step=0.001,
            disabled=use_isolation_forest,
            help="Only used when Isolation Forest is disabled.",
        )
    with c3:
        do_full = st.checkbox("Run full n×n comparison", value=False, help="If off, only diagonal is evaluated.")
        do_blast = st.checkbox("Run NCBI BLAST annotation", value=False,
                               help="If on, pipeline will annotate protein names.")
        st.session_state["do_blast"] = bool(do_blast)
    st.session_state["use_isolation_forest"] = bool(use_isolation_forest)

    run_alignment_button = st.button("Run Alignment", key="run_alignment_btn", type="primary")
    if run_alignment_button:
        ref    = st.session_state.get("ref_seq", "") or ""
        evo    = st.session_state.get("evolved_seq", "") or ""
        start  = st.session_state.get("start_subseq", "") or ""
        missing = [n for n, v in [("Reference", ref), ("Evolved", evo), ("Starting subsequence", start)] if not v]
        if missing:
            st.error(f"Please provide: {', '.join(missing)}.")
            st.stop()

        ref_seq      = ms.get_sequence(ref)
        evolved_seq  = ms.get_sequence(evo)
        start_subseq = ms.get_sequence(start)
        if start_subseq not in ref_seq:
            st.warning("Start subsequence not in **Reference**; using unrotated ref.")
        if start_subseq not in evolved_seq:
            st.warning("Start subsequence not in **Evolved**; using unrotated evolved.")
        rotated_ref     = ms.rotate_sequence(ref_seq, start_subseq)
        rotated_evolved = ms.rotate_sequence(evolved_seq, start_subseq)

        prog_diag     = st.empty()
        prog_compare  = st.empty()
        prog_matrices = st.empty()
        prog_prodigal = st.empty()

        q: "queue.Queue[tuple]" = queue.Queue()
        _result = {}

        def _on_progress(phase: str, i: int, total: int):
            q.put(("progress", phase, int(i), int(total)))

        def _on_log(msg: str):
            q.put(("log", str(msg)))

        def _worker():
            try:
                res = run_pipeline(
                    rotated_ref=rotated_ref,
                    rotated_evolved=rotated_evolved,
                    n_sections=int(n_sections),
                    threshold=float(threshold),
                    do_full=bool(do_full),
                    on_progress=_on_progress,
                    on_log=_on_log,
                    do_blast=bool(st.session_state.get("do_blast", False)),
                )
                _result["val"] = res
                q.put(("done", None))
            except Exception as e:
                import traceback
                tb = traceback.format_exc()
                q.put(("error", e))

        t = threading.Thread(target=_worker, daemon=True)
        t.start()

        log_area = st.empty()
        while True:
            try:
                item = q.get(timeout=0.1)
            except queue.Empty:
                if not t.is_alive():
                    break
                time.sleep(0.01)
                continue

            kind = item[0]
            if kind == "progress":
                _, phase, i, total = item
                total = max(int(total), 1)
                frac  = min(max(int(i), 0) / total, 1.0)
                if phase == "diagonal":
                    prog_compare.empty(); prog_matrices.empty(); prog_prodigal.empty()
                    prog_diag.progress(frac, text=f"Diagonal check: {i}/{total}")
                elif phase == "compare":
                    prog_diag.empty(); prog_matrices.empty(); prog_prodigal.empty()
                    prog_compare.progress(frac, text=f"Comparing chunks: {i}/{total}")
                elif phase == "matrices":
                    prog_diag.empty(); prog_compare.empty(); prog_prodigal.empty()
                    prog_matrices.progress(frac, text=f"Building matrices: {i}/{total}")
                elif phase == "prodigal":
                    prog_diag.empty(); prog_compare.empty(); prog_matrices.empty()
                    prog_prodigal.progress(frac, text=f"Prodigal ORF prediction: {i}/4 sequences")
            elif kind == "log":
                _, msg = item
                with log_area.container():
                    st.write(msg)
            elif kind == "error":
                _, err = item
                prog_diag.empty(); prog_compare.empty()
                prog_matrices.empty(); prog_prodigal.empty()
                st.exception(err)
                st.stop()
            elif kind == "done":
                break
            time.sleep(0.01)

        prog_diag.empty(); prog_compare.empty()
        prog_matrices.empty(); prog_prodigal.empty()

        results = _result.get("val")
        if results is None:
            st.error("Pipeline did not return results.")
            st.stop()

        st.session_state["results"] = results
        st.success("Alignment completed.")

    # ---------- RENDER SECTION ----------
    results = st.session_state.get("results")
    if results:
        with st.expander("Sections (debug)"):
            st.write("Reference sections:",
                     [(i + 1, off, len(sec)) for i, (sec, off) in enumerate(results["ref_sections"])])
            st.write("Evolved sections:",
                     [(i + 1, off, len(sec)) for i, (sec, off) in enumerate(results["evo_sections"])])

        _use_if = st.session_state.get("use_isolation_forest", True)
        if _use_if and results.get("isolation_forest") is not None:
            _if     = results["isolation_forest"]
            _scores = _if["anomaly_scores"]
            _labels = _if["anomaly_labels"]
            _n_anom = int((_labels == -1).sum())

            st.subheader("Isolation Forest — Anomaly Scores per Chunk")
            st.caption(f"{_n_anom} of {len(_labels)} chunks flagged as anomalous (score closer to 1.0 = more anomalous).")

            _fig_if, _ax_if = plt.subplots(figsize=(12, 3))
            _colors = ["#e63946" if lbl == -1 else "#457b9d" for lbl in _labels]
            _ax_if.bar(range(len(_scores)), _scores, color=_colors, width=1.0)
            _ax_if.set_xlabel("Chunk index")
            _ax_if.set_ylabel("Anomaly score")
            _ax_if.set_title("Isolation Forest anomaly score (red = anomalous chunk)")
            _ax_if.set_xlim(0, len(_scores))
            _ax_if.set_ylim(0, 1)
            st.pyplot(_fig_if)
            plt.close(_fig_if)

            with st.expander("Anomaly feature table"):
                st.dataframe(ms.anomalous_rows(_if["feature_df"]), use_container_width=True)
        elif _use_if:
            st.warning("Isolation Forest results unavailable. Re-run the alignment.")
            if results.get("isolation_forest_error"):
                st.error(results["isolation_forest_error"])

        st.write("final_matrix shape:", results["final_matrix"].shape)
        st.write("info_matrix shape:",  results["info_matrix"].shape)

        if results["result_matrix"].size:
            st.write("result_matrix shape:", results["result_matrix"].shape)
        else:
            st.warning("No ORF boundaries computed.")
if page == "Submit to blast":
    st.header("Results")

    res = st.session_state.get("results")
    if not res:
        st.info("Run the alignment first to see results.")
        st.stop()

    rm             = res.get("result_matrix")
    info_matrix    = res.get("info_matrix")
    pairing_scores = res.get("pairing_scores", {})

    if rm is None or (hasattr(rm, '__len__') and len(rm) == 0):
        st.warning("No ORF boundaries computed yet.")
        st.stop()

    # ── build index-grouped structure ───────────────────────────────────────
    from collections import defaultdict
    sites = defaultdict(list)
    for r in rm:
        mut_idx, trk_row, orf, start, end = int(r[0]), int(r[1]), int(r[2]), int(r[3]), int(r[4])
        sites[mut_idx].append({
            "row": trk_row, "orf": orf, "start": start, "end": end
        })

    sorted_sites = sorted(sites.keys())

    def _direction_label(trk_row):
        return "fwd" if ms.is_fwd_row(trk_row) else "rev"

    def _role_label(trk_row):
        return "ref" if ms.is_ref_row(trk_row) else "evo"

    def _orf_label(o):
        return f"{_role_label(o['row'])} F{o['orf']} [{o['start']},{o['end']}] ({_direction_label(o['row'])})"

    # ── dropdown ────────────────────────────────────────────────────────────
    site_labels = []
    for idx in sorted_sites:
        orfs = sites[idx]
        orf_summary = ", ".join(_orf_label(o) for o in orfs)
        site_labels.append(f"Index {idx} — {orf_summary}")

    selected_labels = st.multiselect(
        "Select mutation sites",
        options=site_labels,
        default=[],
        key="mut_site_multiselect",
    )
    selected_indices = [sorted_sites[site_labels.index(l)] for l in selected_labels]

    st.markdown("---")

    # ── BLAST submission ─────────────────────────────────────────────────────
    st.session_state.setdefault("blast_logs", [])
    log_placeholder = None

    blast_names = res.get("blast_names") or {}

    st.markdown("""<style>
    .seqwrap { width:100%; max-width:100%; overflow-x:auto; overflow-y:hidden;
               border:1px solid rgba(0,0,0,0.08); padding:12px 14px; border-radius:10px; background:#fff; }
    .rowline { white-space:nowrap; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
               "Liberation Mono", monospace; font-size:16px; line-height:1.7; }
    .lab { color:#666; margin-right:12px; display:inline-block; min-width:160px; }
    .diffblock { display:inline; }
    .ok { color:#111; }
    .mm { background:#ffd24d; color:#000 !important; border-radius:2px;
          box-shadow: inset 0 0 0 1px rgba(0,0,0,.04); }
    .blastbox { margin-top:10px; padding:12px 14px; border-radius:10px;
                border:1px solid rgba(12,50,135,0.2); background:#f4f7ff; }
    .blasttitle { font-weight:600; color:#123d7a; margin-bottom:4px; }
    .blastbody { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
                 "Liberation Mono", monospace; font-size:14px; line-height:1.5;
                 white-space:pre-wrap; word-break:break-word; color:#0b1533; }
    </style>""", unsafe_allow_html=True)

    log_container = st.container()
    with log_container:
        show_logs = st.checkbox("Show BLAST log output", key="blast_show_logs")
        if show_logs:
            log_placeholder = st.empty()
            existing_log = "\n".join(st.session_state.get("blast_logs", [])) or "(No BLAST activity yet.)"
            log_placeholder.code(existing_log, language="text")

    # ── per-site display ─────────────────────────────────────────────────────
    for mut_idx in selected_indices:
        orfs = sites[mut_idx]
        st.subheader(f"Mutation index {mut_idx}")

        ref_orfs = [o for o in orfs if ms.is_ref_row(o["row"]) and o["orf"] != -1]
        evo_orfs = [o for o in orfs if not ms.is_ref_row(o["row"]) and o["orf"] != -1]

        # ── pairing selector ────────────────────────────────────────────────
        scores = pairing_scores.get(mut_idx)
        if scores:
            st.caption("Multiple pairings possible — ranked by AA overlap:")
            score_rows = []
            for s in scores:
                score_rows.append({
                    "ref":     f"row {s['ref_row']} F{s['ref_frame']} [{s['ref_start']},{s['ref_end']}]",
                    "evo":     f"row {s['evo_row']} F{s['evo_frame']} [{s['evo_start']},{s['evo_end']}]",
                    "matches": s["match_count"],
                    "overlap": s["overlap_len"],
                    "score":   f"{s['score']:.1%}",
                })
            st.dataframe(pd.DataFrame(score_rows), use_container_width=True, hide_index=True)

            # build selector options from scores
            pairing_options = [
                f"ref row {s['ref_row']} F{s['ref_frame']} ↔ evo row {s['evo_row']} F{s['evo_frame']}"
                for s in scores
            ]
            chosen_pairing_label = st.selectbox(
                "Choose pairing to view",
                options=pairing_options,
                key=f"pairing_sel_{mut_idx}",
            )
            chosen_score = scores[pairing_options.index(chosen_pairing_label)]
            ref_orf = {"row": chosen_score["ref_row"], "orf": chosen_score["ref_frame"],
                       "start": chosen_score["ref_start"], "end": chosen_score["ref_end"]}
            evo_orf = {"row": chosen_score["evo_row"], "orf": chosen_score["evo_frame"],
                       "start": chosen_score["evo_start"], "end": chosen_score["evo_end"]}
        elif ref_orfs and evo_orfs:
            # check same direction
            fwd_ref = [o for o in ref_orfs if ms.is_fwd_row(o["row"])]
            fwd_evo = [o for o in evo_orfs if ms.is_fwd_row(o["row"])]
            rev_ref = [o for o in ref_orfs if not ms.is_fwd_row(o["row"])]
            rev_evo = [o for o in evo_orfs if not ms.is_fwd_row(o["row"])]
            # prefer fwd pair, fall back to rev
            if fwd_ref and fwd_evo:
                ref_orf, evo_orf = fwd_ref[0], fwd_evo[0]
            elif rev_ref and rev_evo:
                ref_orf, evo_orf = rev_ref[0], rev_evo[0]
            else:
                ref_orf = ref_orfs[0] if ref_orfs else None
                evo_orf = evo_orfs[0] if evo_orfs else None
        else:
            ref_orf = ref_orfs[0] if ref_orfs else None
            evo_orf = evo_orfs[0] if evo_orfs else None

        # ── sequence display ─────────────────────────────────────────────────
        view_mode = st.radio(
            "View mode",
            ["DNA sequence", "Amino acid sequence"],
            horizontal=True,
            key=f"view_mode_{mut_idx}",
        )

        def _get_dna_slice(orf):
            if orf is None:
                return ""
            nuc_row = {8: 1, 9: 3, 11: 10, 13: 12}[orf["row"]]
            return "".join(info_matrix[nuc_row, orf["start"]:orf["end"]].tolist())

        def _get_aa_slice(orf):
            if orf is None:
                return ""
            aa_row = ms.aa_row_for(orf["row"], orf["orf"])
            frame_num = orf["orf"]
            trk_row = orf["row"]
            result = []
            for i, v in enumerate(info_matrix[aa_row, orf["start"]:orf["end"]:3].tolist()):
                pos = orf["start"] + i * 3
                if pos < info_matrix.shape[1] and ms.frame_active(
                    int(info_matrix[trk_row, pos]), frame_num
                ):
                    result.append(str(v))
            return "".join(result) or "(no ORF in this frame)"

        if view_mode == "DNA sequence":
            ref_slice = _get_dna_slice(ref_orf) if ref_orf else ""
            evo_slice = _get_dna_slice(evo_orf) if evo_orf else ""
            ref_label, evo_label = "REF (DNA)", "EVO (DNA)"
        else:
            ref_slice = _get_aa_slice(ref_orf) if ref_orf else ""
            evo_slice = _get_aa_slice(evo_orf) if evo_orf else ""
            ref_label, evo_label = "REF (AA)", "EVO (AA)"

        # align by coordinate — pad the one that starts later
        if ref_orf and evo_orf:
            ref_start = ref_orf["start"]
            evo_start = evo_orf["start"]
            if view_mode == "DNA sequence":
                if evo_start < ref_start:
                    ref_slice = "-" * (ref_start - evo_start) + ref_slice
                elif ref_start < evo_start:
                    evo_slice = "-" * (evo_start - ref_start) + evo_slice
            else:
                if evo_start < ref_start:
                    pad = (ref_start - evo_start) // 3
                    ref_slice = "-" * pad + ref_slice
                elif ref_start < evo_start:
                    pad = (evo_start - ref_start) // 3
                    evo_slice = "-" * pad + evo_slice
            display_start = min(ref_start, evo_start)
        else:
            display_start = (ref_orf or evo_orf)["start"] if (ref_orf or evo_orf) else 0

        if not ref_slice:
            ref_slice = "(no ORF detected)"
        if not evo_slice:
            evo_slice = "(no ORF detected)"

        st.markdown(
            ms.highlighted_block(ref_slice, evo_slice,
                                  a_label=ref_label, b_label=evo_label,
                                  start_pos=display_start),
            unsafe_allow_html=True,
        )

        # ── BLAST per ORF ────────────────────────────────────────────────────
        st.markdown("**Submit to BLAST:**")
        blast_candidates = []
        if ref_orf and ref_orf["orf"] != -1:
            blast_candidates.append(("ref", ref_orf))
        if evo_orf and evo_orf["orf"] != -1:
            blast_candidates.append(("evo", evo_orf))
        # also list any unpaired ORFs
        paired_rows = set()
        if ref_orf:
            paired_rows.add((ref_orf["row"], ref_orf["orf"]))
        if evo_orf:
            paired_rows.add((evo_orf["row"], evo_orf["orf"]))
        for o in orfs:
            if (o["row"], o["orf"]) not in paired_rows and o["orf"] != -1:
                role = "ref" if ms.is_ref_row(o["row"]) else "evo"
                blast_candidates.append((role, o))

        for role, orf in blast_candidates:
            blast_key = f"blast_{mut_idx}_{orf['row']}_{orf['orf']}"
            col1, col2 = st.columns([3, 1])
            col1.caption(_orf_label(orf))
            if col2.button("Submit", key=f"btn_{blast_key}"):
    # build orf_code: frame 1-3 for forward, 4-6 for revcomp
                is_rev = not ms.is_fwd_row(orf["row"])
                orf_code = orf["orf"] + (3 if is_rev else 0)
                
                single_matrix = np.array([[
                    orf_code, orf["start"], orf["end"]
                ]], dtype=int)
                
                def _log_blast(msg):
                    st.session_state["blast_logs"] = \
                        st.session_state.get("blast_logs", []) + [msg]
                    if log_placeholder:
                        log_placeholder.code(
                            "\n".join(st.session_state["blast_logs"]),
                            language="text"
                        )
                with st.spinner(f"BLASTing {_orf_label(orf)}…"):
                    try:
                        blast_result = process_result_matrix(
                            info_matrix, single_matrix,
                            do_blast=True, on_log=_log_blast
                        )
                        hit = blast_result.get(0, "No BLAST hits found.")
                        key_str = f"{mut_idx}_{orf['row']}_{orf['orf']}"
                        existing = dict(res.get("blast_names") or {})
                        existing[key_str] = hit
                        res["blast_names"] = existing
                        st.session_state["results"]["blast_names"] = existing
                        st.success(f"BLAST complete: {hit[:80]}…" if len(hit) > 80 else f"BLAST complete: {hit}")
                    except Exception as exc:
                        st.error(f"BLAST failed: {exc}")
            # show existing result if present
            key_str = f"{mut_idx}_{orf['row']}_{orf['orf']}"
            if key_str in blast_names:
                safe = html.escape(blast_names[key_str]).replace("\n", "<br>")
                st.markdown(
                    f"<div class='blastbox'><div class='blasttitle'>{_orf_label(orf)}</div>"
                    f"<div class='blastbody'>{safe}</div></div>",
                    unsafe_allow_html=True,
                )

        st.markdown("---")
if page == "BLAST results":
   st.header("BLAST Results")
   res = st.session_state.get("results")
   blast_names = {}
   if res:
       blast_names = res.get("blast_names") or {}


   if not blast_names:
       st.info("No BLAST submissions yet.")
   else:
       st.caption("NCBI BLAST annotations by mutation number.")
       for idx in sorted(blast_names):
           name = blast_names.get(idx) or "No BLAST hits found."
           safe_name = html.escape(name).replace("\n", "<br>")
           st.markdown(
               f"<div class='blastbox'><div class='blasttitle'>Mutation {idx + 1}</div>"
               f"<div class='blastbody'>{safe_name}</div></div>",
               unsafe_allow_html=True,
           )


if page == "View":
    st.header("Selected mutations summary")
    res = st.session_state.get("results")
    if not res:
        st.info("Run the alignment in Tab 1 to populate this view.")
        st.stop()

    rm             = res.get("result_matrix")
    info_matrix    = res.get("info_matrix")
    pairing_scores = res.get("pairing_scores", {})

    if rm is None or (hasattr(rm, '__len__') and len(rm) == 0) or info_matrix is None:
        st.info("No mutation data available yet.")
        st.stop()

    from collections import defaultdict
    sites = defaultdict(list)
    for r in rm:
        mut_idx, trk_row, orf, start, end = int(r[0]), int(r[1]), int(r[2]), int(r[3]), int(r[4])
        sites[mut_idx].append({"row": trk_row, "orf": orf, "start": start, "end": end})

    sorted_sites = sorted(sites.keys())

    def _direction_label(trk_row):
        return "fwd" if ms.is_fwd_row(trk_row) else "rev"

    def _role_label(trk_row):
        return "ref" if ms.is_ref_row(trk_row) else "evo"

    def _orf_label(o):
        return f"{_role_label(o['row'])} F{o['orf']} [{o['start']},{o['end']}] ({_direction_label(o['row'])})"

    site_labels = []
    for idx in sorted_sites:
        orfs = sites[idx]
        orf_summary = ", ".join(_orf_label(o) for o in orfs)
        site_labels.append(f"Index {idx} — {orf_summary}")

    selected_labels = st.multiselect(
        "Select mutation sites",
        options=site_labels,
        default=[],
        key="view_site_multiselect",
    )
    selected_indices = [sorted_sites[site_labels.index(l)] for l in selected_labels]

    if not selected_indices:
        st.info("Select at least one mutation site above to view sequences.")
        st.stop()

    # ── controls ─────────────────────────────────────────────────────────────
    controls_left, controls_mid, controls_right = st.columns([0.33, 0.27, 0.4])
    with controls_left:
        show_blast = st.checkbox("Show BLAST results", key="view_show_blast")
    with controls_mid:
        pdf_mode = st.checkbox("PDF friendly layout", key="view_pdf_mode")
    mode_options = ["DNA sequence", "Amino acid sequence"]
    with controls_right:
        mode = st.radio(
            "Sequence display",
            mode_options,
            horizontal=True,
            key="view_mode",
        )

    blast_names = res.get("blast_names") or {}
    wrap_width = 70 if pdf_mode else None

    st.markdown(
        """
        <style>
        .pdf-block { margin-bottom:16px; padding:14px 16px; border-radius:12px;
                     border:1px solid rgba(0,0,0,0.12); background:#fff; }
        .pdf-title { font-weight:600; font-size:18px; color:#111; margin-bottom:10px; }
        .seqwrap { width:100%; max-width:100%; overflow-x:auto; overflow-y:hidden;
                   border:1px solid rgba(0,0,0,0.08); padding:12px 14px; border-radius:10px;
                   background:#fff; display:block; }
        .rowline { white-space:nowrap;
                   font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono",monospace;
                   font-size:14px; line-height:1.7; color:#111; display:block; width:max-content; }
        .lab { color:#666; margin-right:12px; display:inline-block; min-width:160px; }
        .diffblock { display:inline; color:#111; }
        .ok { color:#111; }
        .mm { background:#ffd24d; color:#000 !important; border-radius:2px; }
        .blastbox { margin-top:10px; padding:12px 14px; border-radius:10px;
                    border:1px solid rgba(12,50,135,0.2); background:#f4f7ff; }
        .blasttitle { font-weight:600; color:#123d7a; margin-bottom:4px; }
        .blastbody { font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono",monospace;
                     font-size:14px; line-height:1.5; white-space:pre-wrap;
                     word-break:break-word; color:#0b1533; }
        .seqwrap.pdf .pdf-pair { margin-bottom:14px; }
        .seqwrap.pdf .pdf-pair:last-child { margin-bottom:0; }
        @media print {
            .pdf-block { page-break-before:always; margin:0 0 8px 0; padding:10px 12px; }
            .pdf-block:first-child { page-break-before:auto; }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    def _get_dna_slice(orf):
        if orf is None:
            return ""
        nuc_row = {8: 1, 9: 3, 11: 10, 13: 12}[orf["row"]]
        return "".join(info_matrix[nuc_row, orf["start"]:orf["end"]].tolist())

    def _get_aa_slice(orf):
        if orf is None:
            return ""
        aa_row = ms.aa_row_for(orf["row"], orf["orf"])
        frame_num = orf["orf"]
        trk_row = orf["row"]
        result = []
        for i, v in enumerate(info_matrix[aa_row, orf["start"]:orf["end"]:3].tolist()):
            pos = orf["start"] + i * 3
            if pos < info_matrix.shape[1] and ms.frame_active(
                int(info_matrix[trk_row, pos]), frame_num
            ):
                result.append(str(v))
        return "".join(result) or "(no ORF in this frame)"

    def _chunk_size(length: int) -> int:
        return wrap_width if wrap_width else max(length, 1)

    def _build_highlight_html(ref_seq: str, evo_seq: str, ref_lbl: str, evo_lbl: str,
                               start_pos: int, pdf: bool, chunk: int) -> str:
        if not ref_seq and not evo_seq:
            return "<div class='seqwrap'>No sequence data available.</div>"
        container_class = "seqwrap pdf" if pdf else "seqwrap"
        seq_len = max(len(ref_seq), len(evo_seq), 1)
        blocks = []
        for offset in range(0, seq_len, chunk):
            chunk_ref = ref_seq[offset:offset + chunk]
            chunk_evo = evo_seq[offset:offset + chunk]
            base_index = start_pos + offset
            def _spanify(a: str, b: str) -> str:
                cls = "mm" if b and a != b else "ok"
                return f"<span class='{cls}'>{html.escape(a)}</span>"
            la = "".join(_spanify(a, chunk_evo[i] if i < len(chunk_evo) else "") for i, a in enumerate(chunk_ref))
            lb = "".join(_spanify(b, chunk_ref[i] if i < len(chunk_ref) else "") for i, b in enumerate(chunk_evo))
            lab_a = f"{ref_lbl} {base_index:>7d}: "
            lab_b = f"{evo_lbl} {base_index:>7d}: "
            pair_start = "<div class='pdf-pair'>" if pdf else ""
            pair_end   = "</div>" if pdf else ""
            blocks.append(
                f"{pair_start}"
                f"<div class='rowline'><span class='lab'>{html.escape(lab_a)}</span>"
                f"<span class='diffblock'>{la}</span></div>"
                f"<div class='rowline'><span class='lab'>{html.escape(lab_b)}</span>"
                f"<span class='diffblock'>{lb}</span></div>"
                f"{pair_end}"
            )
        return f"<div class='{container_class}'>{''.join(blocks)}</div>"

    def _build_copy_block(header: str, ref_lbl: str, evo_lbl: str,
                           ref_seq: str, evo_seq: str, start_pos: int, chunk: int) -> str:
        if not ref_seq and not evo_seq:
            return f"{header}\n{ref_lbl}: (empty)\n{evo_lbl}: (empty)"
        seq_len = max(len(ref_seq), len(evo_seq), 1)
        lines = [header]
        for offset in range(0, seq_len, chunk):
            chunk_ref = ref_seq[offset:offset + chunk]
            chunk_evo = evo_seq[offset:offset + chunk]
            base_index = start_pos + offset
            lines.append(f"{ref_lbl} {base_index:>7d}: {chunk_ref}")
            lines.append(f"{evo_lbl} {base_index:>7d}: {chunk_evo}")
        return "\n".join(lines)

    copy_blocks: list[str] = []
    st.caption(f"{len(selected_indices)} mutation site(s) selected")

    for pos, mut_idx in enumerate(selected_indices):
        orfs = sites[mut_idx]
        st.subheader(f"Mutation index {mut_idx}")

        ref_orfs = [o for o in orfs if ms.is_ref_row(o["row"]) and o["orf"] != -1]
        evo_orfs = [o for o in orfs if not ms.is_ref_row(o["row"]) and o["orf"] != -1]

        # ── pairing selector ─────────────────────────────────────────────────
        scores = pairing_scores.get(mut_idx)
        if scores:
            st.caption("Multiple pairings possible — ranked by AA overlap:")
            score_rows = []
            for s in scores:
                score_rows.append({
                    "ref":     f"row {s['ref_row']} F{s['ref_frame']} [{s['ref_start']},{s['ref_end']}]",
                    "evo":     f"row {s['evo_row']} F{s['evo_frame']} [{s['evo_start']},{s['evo_end']}]",
                    "matches": s["match_count"],
                    "overlap": s["overlap_len"],
                    "score":   f"{s['score']:.1%}",
                })
            st.dataframe(pd.DataFrame(score_rows), use_container_width=True, hide_index=True)
            pairing_options = [
                f"ref row {s['ref_row']} F{s['ref_frame']} ↔ evo row {s['evo_row']} F{s['evo_frame']}"
                for s in scores
            ]
            chosen_pairing_label = st.selectbox(
                "Choose pairing to view",
                options=pairing_options,
                key=f"view_pairing_sel_{mut_idx}",
            )
            chosen_score = scores[pairing_options.index(chosen_pairing_label)]
            ref_orf = {"row": chosen_score["ref_row"], "orf": chosen_score["ref_frame"],
                       "start": chosen_score["ref_start"], "end": chosen_score["ref_end"]}
            evo_orf = {"row": chosen_score["evo_row"], "orf": chosen_score["evo_frame"],
                       "start": chosen_score["evo_start"], "end": chosen_score["evo_end"]}
        elif ref_orfs and evo_orfs:
            fwd_ref = [o for o in ref_orfs if ms.is_fwd_row(o["row"])]
            fwd_evo = [o for o in evo_orfs if ms.is_fwd_row(o["row"])]
            rev_ref = [o for o in ref_orfs if not ms.is_fwd_row(o["row"])]
            rev_evo = [o for o in evo_orfs if not ms.is_fwd_row(o["row"])]
            if fwd_ref and fwd_evo:
                ref_orf, evo_orf = fwd_ref[0], fwd_evo[0]
            elif rev_ref and rev_evo:
                ref_orf, evo_orf = rev_ref[0], rev_evo[0]
            else:
                ref_orf = ref_orfs[0] if ref_orfs else None
                evo_orf = evo_orfs[0] if evo_orfs else None
        else:
            ref_orf = ref_orfs[0] if ref_orfs else None
            evo_orf = evo_orfs[0] if evo_orfs else None

        # ── sequence retrieval ───────────────────────────────────────────────
        if mode == "DNA sequence":
            ref_slice = _get_dna_slice(ref_orf) if ref_orf else ""
            evo_slice = _get_dna_slice(evo_orf) if evo_orf else ""
            ref_label, evo_label = "REF (DNA)", "EVO (DNA)"
        else:
            ref_slice = _get_aa_slice(ref_orf) if ref_orf else ""
            evo_slice = _get_aa_slice(evo_orf) if evo_orf else ""
            ref_label, evo_label = "REF (AA)", "EVO (AA)"

        # align by coordinate — pad the one that starts later
        if ref_orf and evo_orf:
            ref_start = ref_orf["start"]
            evo_start = evo_orf["start"]
            if mode == "DNA sequence":
                if evo_start < ref_start:
                    ref_slice = "-" * (ref_start - evo_start) + ref_slice
                elif ref_start < evo_start:
                    evo_slice = "-" * (evo_start - ref_start) + evo_slice
            else:
                if evo_start < ref_start:
                    ref_slice = "-" * ((ref_start - evo_start) // 3) + ref_slice
                elif ref_start < evo_start:
                    evo_slice = "-" * ((evo_start - ref_start) // 3) + evo_slice
            display_start = min(ref_start, evo_start)
        else:
            display_start = (ref_orf or evo_orf)["start"] if (ref_orf or evo_orf) else 0

        if not ref_slice:
            ref_slice = "(no ORF detected)"
        if not evo_slice:
            evo_slice = "(no ORF detected)"

        # ── render ───────────────────────────────────────────────────────────
        chunk = _chunk_size(max(len(ref_slice), len(evo_slice)))
        header = f"Mutation index {mut_idx}"
        highlight_html = _build_highlight_html(
            ref_slice, evo_slice, ref_label, evo_label, display_start, pdf_mode, chunk
        )

        blast_section = ""
        if show_blast:
            blast_key = (
                f"{mut_idx}_{ref_orf['row']}_{ref_orf['orf']}" if ref_orf
                else f"{mut_idx}_{evo_orf['row']}_{evo_orf['orf']}" if evo_orf
                else None
            )
            prot_text = blast_names.get(blast_key) if blast_key else None
            if prot_text:
                safe_prot = html.escape(str(prot_text)).replace("\n", "<br>")
                blast_section = (
                    f"<div class='blastbox'><div class='blasttitle'>BLAST result</div>"
                    f"<div class='blastbody'>{safe_prot}</div></div>"
                )
            else:
                blast_section = (
                    "<div class='blastbox'><div class='blasttitle'>BLAST result</div>"
                    "<div class='blastbody'>Submit to BLAST on page 2 to see protein results.</div></div>"
                )

        style_attr = "page-break-before: always;" if pdf_mode and pos > 0 else ""
        block_html = (
            f"<div class='pdf-block' style=\"{style_attr}\">"
            f"<div class='pdf-title'>{html.escape(header)}</div>"
            f"{highlight_html}"
            f"{blast_section if show_blast else ''}"
            f"</div>"
        )
        st.markdown(block_html, unsafe_allow_html=True)
        if not pdf_mode and pos < len(selected_indices) - 1:
            st.divider()

        copy_text = _build_copy_block(
            header, ref_label, evo_label, ref_slice, evo_slice, display_start, chunk
        )
        copy_blocks.append(copy_text)

    if copy_blocks:
        combined_joiner = "\n\f\n" if pdf_mode else "\n\n"
        combined_text = combined_joiner.join(copy_blocks)
        copy_html = f"""
        <div style="display:flex; justify-content:flex-end; margin: 8px 0 12px;">
            <button id="view-copy-all"
                    style="background:#1b73e8; color:#fff; border:0; border-radius:6px;
                           padding:6px 14px; font-weight:600; cursor:pointer;">
                Copy all
            </button>
        </div>
        <script>
        (function() {{
            const button = document.getElementById('view-copy-all');
            if (!button) return;
            const originalText = button.textContent;
            button.addEventListener('click', async () => {{
                try {{
                    await navigator.clipboard.writeText({json.dumps(combined_text)});
                    button.textContent = 'Copied!';
                    button.disabled = true;
                    setTimeout(() => {{
                        button.textContent = originalText;
                        button.disabled = false;
                    }}, 1500);
                }} catch (err) {{
                    console.error('Clipboard copy failed', err);
                }}
            }});
        }})();
        </script>
        """
        components.html(copy_html, height=60)


if page == "Debug":
    st.header("Debug")
    st.markdown("""
    <style>
    .seqwrap { 
        width:100%; max-width:100%; overflow-x:auto; overflow-y:hidden;
        border:1px solid rgba(0,0,0,0.08); padding:12px 14px; border-radius:10px; 
        background:#fff; display:block;
    }
    .rowline { 
        white-space:nowrap; 
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
        font-size:14px; line-height:1.7; color:#111; display:block; width:max-content;
    }
    .diffblock { display:inline; color:#111; }
    </style>
    """, unsafe_allow_html=True)
    res = st.session_state.get("results")
    if not res:
        st.info("Run the alignment in Tab 1 to populate debug data.")
    else:
        info_matrix = res.get("info_matrix")
        if info_matrix is None or getattr(info_matrix, "size", 0) == 0:
            st.info("info_matrix is not available yet.")
        else:
            info_matrix = np.asarray(info_matrix)
            n_rows, n_cols = info_matrix.shape

            st.session_state.setdefault("debug_start_idx", 0)
            st.session_state.setdefault("debug_end_idx", min(n_cols, 60))

            c_start, c_end = st.columns(2)
            max_start = max(n_cols - 1, 0)
            start_default = st.session_state.get("debug_start_idx", 0)
            start_idx = int(
                c_start.number_input(
                    "Start index (inclusive)",
                    min_value=0,
                    max_value=max_start,
                    value=min(max(start_default, 0), max_start),
                )
            )

            max_end = max(n_cols, 1)
            end_min_required = min(max_end, start_idx + 1)
            stored_end = st.session_state.get("debug_end_idx", end_min_required)
            end_default = min(max(stored_end, end_min_required), max_end)
            end_idx = int(
                c_end.number_input(
                    "End index (exclusive)",
                    min_value=end_min_required,
                    max_value=max_end,
                    value=end_default,
                )
            )

            if end_idx <= start_idx:
                end_idx = min(start_idx + 1, max_end)

            end_idx = max(end_idx, start_idx + 1)
            end_idx = min(end_idx, n_cols)

            st.session_state["debug_start_idx"] = start_idx
            st.session_state["debug_end_idx"] = end_idx

            if start_idx >= end_idx:
                st.warning("End index must be greater than start index to display a slice.")
            else:
                row_options = [
                    f"Row {idx} — {desc}" for idx, desc in INFO_MATRIX_CHEATSHEET
                ]

                selected_labels = st.multiselect(
                    "Select rows to view",
                    options=row_options,
                    default=[],
                    key="debug_row_multiselect",
                )

                selected_rows = [
                    idx for idx, desc in INFO_MATRIX_CHEATSHEET
                    if f"Row {idx} — {desc}" in selected_labels
                ]

                if not selected_rows:
                    st.info("Select at least one row to display its values.")
                else:
                    st.subheader(f"Selected row slices [{start_idx}:{end_idx})")
                    for idx in selected_rows:
                        row_slice = info_matrix[idx, start_idx:end_idx]
                        string_vals = ["" if v is None else str(v) for v in row_slice.tolist()]
                        non_empty = [s for s in string_vals if s]
                        if non_empty and all(len(s) == 1 for s in non_empty):
                            slice_values = "".join(string_vals)
                        else:
                            slice_values = " ".join(string_vals)
                        label_desc = next((desc for ridx, desc in INFO_MATRIX_CHEATSHEET if ridx == idx), "")
                        st.caption(f"Row {idx} — {label_desc}" if label_desc else f"Row {idx}")
                        st.markdown(
                            f"<div class='seqwrap'>"
                            f"<div class='rowline'><span class='diffblock' style=\"color:#111;\">{html.escape(slice_values)}</span></div>"
                            f"</div>",
                            unsafe_allow_html=True,
                        )
            st.markdown("---")
            st.subheader("REPL")
            st.caption(
                "Available variables: `info_matrix`, `results`, `np`, `pd`, `final_matrix`. "
                "Use `print()` or assign to `_out` to display output."
            )

            repl_code = st.text_area(
                "Python expression or block",
                height=160,
                placeholder="# examples:\nprint(info_matrix[1, 500:520])\nprint(np.nonzero(info_matrix[5, :]))\n_out = pd.Series(info_matrix[5,:]).value_counts()",
                key="repl_input",
            )

            if st.button("Run", key="repl_run"):
                _namespace = {
                    "info_matrix": info_matrix,
                    "final_matrix": res.get("final_matrix"),
                    "results": res,
                    "np": np,
                    "pd": pd,
                    "_out": None,
                }
                import io, contextlib
                _stdout = io.StringIO()
                try:
                    with contextlib.redirect_stdout(_stdout):
                        exec(repl_code, _namespace)
                    output = _stdout.getvalue()
                    if _namespace.get("_out") is not None:
                        _out = _namespace["_out"]
                        if isinstance(_out, pd.DataFrame):
                            st.dataframe(_out, use_container_width=True)
                        elif isinstance(_out, pd.Series):
                            st.dataframe(_out.to_frame(), use_container_width=True)
                        else:
                            st.write(_out)
                    if output:
                        st.code(output, language="text")
                    if _namespace.get("_out") is None and not output:
                        st.info("No output. Use print() or assign to `_out`.")
                except Exception as e:
                    st.exception(e)

            st.markdown("---")
            st.subheader("info_matrix cheatsheet")
            cheat_df = pd.DataFrame(INFO_MATRIX_CHEATSHEET, columns=["Row", "Description"])
            # Drop this temporarily anywhere after info_matrix is fully built

            st.dataframe(cheat_df, use_container_width=True, hide_index=True)

def _worker():
    try:
        res = run_pipeline(
            rotated_ref=rotated_ref,
            rotated_evolved=rotated_evolved,
            n_sections=int(n_sections),
            threshold=float(threshold),
            do_full=bool(do_full),
            on_progress=_on_progress,
            on_log=_on_log,
            do_blast=bool(st.session_state.get("do_blast", False)),
        )
        _result["val"] = res
        q.put(("done", None))
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        q.put(("error", e))