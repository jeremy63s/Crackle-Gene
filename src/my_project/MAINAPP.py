from __future__ import annotations
# --- Matplotlib: force headless backend for Streamlit ---
# --- HARD STOP if any code tries to import Tk / GUI backends ---
import os
import html
import json
import pandas as pd


import threading, queue, time
os.environ.setdefault("MPLBACKEND", "Agg")  # force headless BEFORE pyplot
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List
from my_project.pypeline import run_pipeline
#from my_project.pipeline import run_pipeline
from my_project.visualization import plot_interactive_matrix


# Your heavy deps (make sure they are installed in your venv)
import parasail
import string

# ---------------------------------------------------------------




os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
# --------------------------------------------------------
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
st.caption(f"Matplotlib backend: {matplotlib.get_backend()}")





import re
import streamlit as st
import streamlit.components.v1 as components
import microservices as ms  # <- our pure helpers (no Streamlit inside)
from typing import Optional, Tuple
import streamlit as st
import microservices as ms
#from io_helpers import get_sequence
#from sequence_utils import rotate_sequence




# MAINAPP.py
#from __future__ import annotations
import streamlit as st
import microservices as ms  # pure helpers
from my_project.blast import process_result_matrix
import streamlit as st
import microservices as ms  # pure helpers
from my_project.blast import process_result_matrix


INFO_MATRIX_CHEATSHEET = [
   (0, "Alignment data (chunkwise)"),
   (1, "Reference nucleotide sequence"),
   (2, "Reference nucleotide index"),
   (3, "Evolved nucleotide sequence"),
   (4, "Evolved nucleotide index"),
   (5, "Frameshift counter (del:+1, ins:-1)"),
   (6, "Mutation type (in/del:0, point:-1)"),
   (7, "Global index"),
   (8, "Amino acids in ref (naive ORF def)"),
   (9, "Amino acids in evo (naive ORF def)"),
   (10, "RevComp of reference seq"),
   (11, "AA in RevComp ref (naive ORF)"),
   (12, "RevComp of evolved seq"),
   (13, "AA in RevComp evolved (naive ORF)"),
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
tab1, tab2, tab3, tab4, tab5 = st.tabs(
   ["Upload sequences", "Mutations", "View", "Blast results", "Debug"]
)


with tab1:
   st.header("Upload / Paste Sequences")
   st.caption(f"BLAST flag = {bool(st.session_state.get('do_blast', False))}")
   st.caption(f"run_pipeline from {run_pipeline.__module__}")  # should show 'pypeline'


   ref_h, ref_seq = sequence_picker("Reference sequence", key_prefix="ref",   session_key="ref_seq")
   st.divider()
   evo_h, evolved_seq = sequence_picker("Query sequence", key_prefix="evo", session_key="evolved_seq")
   st.divider()
   start_h, start_subseq = sequence_picker("Starting subsequence", key_prefix="start", session_key="start_subseq")


   st.divider()
   # --- NEW: pipeline controls ---
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
   with tab1:
       # ... your inputs ...
       run_alignment_button = st.button("Run Alignment", key="run_alignment_btn", type="primary")
       if run_alignment_button:
           # gather inputs
           ref = st.session_state.get("ref_seq", "") or ""
           evo = st.session_state.get("evolved_seq", "") or ""
           start = st.session_state.get("start_subseq", "") or ""
           missing = [n for n, v in [("Reference", ref), ("Evolved", evo), ("Starting subsequence", start)] if not v]
           if missing:
               st.error(f"Please provide: {', '.join(missing)}.")
               st.stop()


           # normalize + rotate
           ref_seq = ms.get_sequence(ref)
           evolved_seq = ms.get_sequence(evo)
           start_subseq = ms.get_sequence(start)
           if start_subseq not in ref_seq:
               st.warning("Start subsequence not in **Reference**; using unrotated ref.")
           if start_subseq not in evolved_seq:
               st.warning("Start subsequence not in **Evolved**; using unrotated evolved.")
           rotated_ref = ms.rotate_sequence(ref_seq, start_subseq)
           rotated_evolved = ms.rotate_sequence(evolved_seq, start_subseq)


           # === UI progress bars (one per phase) ===
           prog_diag = st.progress(0, text="Diagonal check…")  # <-- NEW
           prog_compare = st.progress(0, text="Comparing chunks…")
           prog_matrices = st.progress(0, text="Building matrices…")


           import threading, queue, time


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
                       on_progress=_on_progress,  # <- critical!
                       on_log=_on_log,
                       do_blast=bool(st.session_state.get("do_blast", False)),
                   )
                   _result["val"] = res
                   q.put(("done", None))
               except Exception as e:
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
                   frac = min(max(int(i), 0) / total, 1.0)
                   if phase == "diagonal":  # <-- NEW
                       prog_diag.progress(frac, text=f"Diagonal check: {i}/{total}")
                   elif phase == "compare":






                       prog_compare.progress(frac, text=f"Comparing chunks: {i}/{total}")
                   elif phase == "matrices":
                       prog_matrices.progress(frac, text=f"Building matrices: {i}/{total}")
               elif kind == "log":
                   _, msg = item
                   with log_area.container():
                       st.write(msg)
               elif kind == "error":
                   _, err = item
                   prog_diag.empty();
                   prog_compare.empty();
                   prog_matrices.empty()
                   st.exception(err);
                   st.stop()
               elif kind == "done":
                   break
               time.sleep(0.01)
           prog_diag.empty()
           prog_compare.empty();
           prog_matrices.empty()


           results = _result.get("val")
           if results is None:
               st.error("Pipeline did not return results.");
               st.stop()


           st.session_state["results"] = results
           st.success("Alignment completed.")


   # ---------- RENDER SECTION (no computation) ----------
   results = st.session_state.get("results")
   if results:
       with st.expander("Sections (debug)"):
           st.write("Reference sections:",
                    [(i + 1, off, len(sec)) for i, (sec, off) in enumerate(results["ref_sections"])])
           st.write("Evolved sections:",
                    [(i + 1, off, len(sec)) for i, (sec, off) in enumerate(results["evo_sections"])])


       num_chunks = len(results["evo_sections"])
       _use_if = st.session_state.get("use_isolation_forest", True)

       if _use_if and results.get("isolation_forest") is not None:
           _if = results["isolation_forest"]
           _scores = _if["anomaly_scores"]   # shape (n_sections,)
           _labels = _if["anomaly_labels"]   # -1 anomaly, 1 normal
           _n_anom = int((_labels == -1).sum())

           st.subheader("Isolation Forest — Anomaly Scores per Chunk")
           st.caption(
               f"{_n_anom} of {len(_labels)} chunks flagged as anomalous (score closer to 1.0 = more anomalous)."
           )

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
           st.warning(
               "Isolation Forest results unavailable. "
               "Re-run the alignment — the error detail will appear in the log area during the run."
           )
           if results.get("isolation_forest_error"):
               st.error(results["isolation_forest_error"])

       st.write("final_matrix shape:", results["final_matrix"].shape)
       st.write("info_matrix shape:", results["info_matrix"].shape)
       if results["result_matrix"].size:
           st.write("result_matrix shape:", results["result_matrix"].shape)
       else:
           st.warning("No ORF boundaries computed.")


with tab2:
   st.header("Results")

   res = st.session_state.get("results")
   if not res:
       st.info("Run the alignment first to see results.")
       st.stop()

   rm = res.get("result_matrix")
   info_matrix = res.get("info_matrix")
   if rm is None or getattr(rm, "size", 0) == 0:
       st.warning("No ORF boundaries computed yet.")
       st.stop()

   import pandas as pd
   df = pd.DataFrame(rm, columns=["code", "start", "end"]).astype(int)

   st.session_state.setdefault("blast_logs", [])
   log_placeholder = None

   st.markdown("""
   <style>
   .seqwrap { width:100%; max-width:100%; overflow-x:auto; overflow-y:hidden;
              border:1px solid rgba(0,0,0,0.08); padding:12px 14px; border-radius:10px; background:#fff; }
   .rowline { white-space:nowrap; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
              font-size:16px; line-height:1.7; }
   .lab { color:#666; margin-right:12px; display:inline-block; min-width:160px; }
   .diffblock { display:inline; }
   .ok { color:#111; }
   .mm { background:#ffd24d; color:#000 !important; border-radius:2px; box-shadow: inset 0 0 0 1px rgba(0,0,0,.04); }
   .legend { color:#666; font-size:13px; margin-top:6px; }
   .blastbox { margin-top:10px; padding:12px 14px; border-radius:10px;
               border:1px solid rgba(12, 50, 135, 0.2); background:#f4f7ff; }
   .blasttitle { font-weight:600; color:#123d7a; margin-bottom:4px; }
   .blastbody { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
                font-size:14px; line-height:1.5; white-space:pre-wrap; word-break:break-word;
                color:#0b1533; }
   </style>
   """, unsafe_allow_html=True)

   def _select_row(i: int):
       st.session_state["selected_result_row"] = i

   top_left, top_right = st.columns([1, 1])
   with top_right:
       mode = st.radio("View mode", ["DNA sequence", "Amino acid sequence"],
                       index=0, horizontal=True, key="view_mode_tab2")

   details_container = st.container()

   st.subheader("Select mutated ORFs to view")

   fs_row = info_matrix[5, :]

   # --- multiselect replacing checkboxes ---
   mutation_options = []
   for idx, (code, start, end) in enumerate(df.itertuples(index=False, name=None)):
       try:
           diffs = ms.diff_indices(int(code), int(start), int(end), info_matrix)
       except Exception:
           diffs = []
       types = [ms.mut_type_from_frameshift(int(gi), fs_row) for gi in range(int(start), int(end))]
       non_point = [t for t in types if t != "point"]
       mut_type = max(set(non_point), key=non_point.count) if non_point else "point"
       length = int(end - start)
       label = f"Mutation {idx + 1} — code {code}, start {start}, end {end}, length {length}, diffs {len(diffs)}, type {mut_type}"
       mutation_options.append(label)

   selected_labels = st.multiselect(
       "Select mutations",
       options=mutation_options,
       default=[],
       key="mut_multiselect",
   )
   selected_rows = [mutation_options.index(l) for l in selected_labels]

   # view buttons — one per selected mutation
   for idx in selected_rows:
       st.button(f"View mutation {idx + 1}", key=f"rm_row_btn_{idx}",
                 on_click=_select_row, args=(idx,))

   # keep mut_checks in sync for tab 3/4 compatibility
   st.session_state["mut_checks"] = {i: (i in selected_rows) for i in range(len(df))}

   st.markdown("---")
   st.caption(f"{len(selected_rows)} mutation(s) selected")

   submit_pressed = st.button("Submit to BLAST", key="blast_submit", type="primary")
   log_container = st.container()
   with log_container:
       show_logs = st.checkbox(
           "Show BLAST log output",
           help="Display read-only streaming output while BLAST requests run.",
           key="blast_show_logs",
       )
       if show_logs:
           log_placeholder = st.empty()
           existing_log = "\n".join(st.session_state.get("blast_logs", [])) or "(No BLAST activity yet.)"
           log_placeholder.code(existing_log, language="text")
       else:
           log_placeholder = None

   if submit_pressed:
       if not selected_rows:
           st.warning("Select at least one mutation before submitting to BLAST.")
       else:
           subset_rows = df.iloc[selected_rows][["code", "start", "end"]].astype(int)
           subset_matrix = subset_rows.to_numpy()

           log_buffer: list[str] = []

           def _record_log(msg: str) -> None:
               log_buffer.append(msg)
               st.session_state["blast_logs"] = log_buffer.copy()
               if log_placeholder is not None:
                   log_placeholder.code("\n".join(st.session_state["blast_logs"]), language="text")

           st.session_state["blast_logs"] = []
           if log_placeholder is not None:
               log_placeholder.code("Submitting requests to NCBI BLAST…", language="text")

           with st.spinner("Submitting selected sequences to NCBI BLAST…"):
               try:
                   blast_results = process_result_matrix(
                       info_matrix,
                       subset_matrix,
                       do_blast=True,
                       on_log=_record_log,
                   )
               except Exception as exc:
                   st.error(f"BLAST request failed: {exc}")
                   if log_placeholder is not None:
                       log_placeholder.code(f"Error: {exc}", language="text")
               else:
                   st.session_state["blast_queue"] = subset_rows.apply(tuple, axis=1).tolist()
                   existing = dict(res.get("blast_names") or {})
                   for local_idx, global_idx in enumerate(selected_rows):
                       name = blast_results.get(local_idx)
                       if name:
                           existing[global_idx] = name
                       else:
                           existing.setdefault(global_idx, "No BLAST hits found.")
                   res["blast_names"] = existing
                   st.session_state["results"]["blast_names"] = existing
                   hit_count = sum(1 for name in blast_results.values() if name)
                   if hit_count:
                       st.success(f"Retrieved BLAST annotations for {hit_count} mutation(s).")
                   else:
                       st.info("BLAST completed but no hits were found for the selected mutations.")
                   if log_placeholder is not None:
                       final_log = "\n".join(st.session_state.get("blast_logs", []))
                       if final_log.strip():
                           log_placeholder.code(final_log, language="text")
                       else:
                           log_placeholder.code("BLAST completed.", language="text")

   with details_container:
       selected_idx = st.session_state.get("selected_result_row")
       if selected_idx is not None:
           try:
               code, start, end = map(int, df.iloc[selected_idx][["code", "start", "end"]].tolist())
           except Exception:
               code = start = end = None

           if code is not None:
               if mode == "DNA sequence":
                   r_row, e_row = (1, 3) if code in (1, 2, 3) else (10, 12)
                   ref_slice = "".join(info_matrix[r_row, start:end].tolist())
                   evo_slice = "".join(info_matrix[e_row, start:end].tolist())
                   st.subheader("DNA comparison")
                   st.markdown(
                       ms.highlighted_block(ref_slice, evo_slice,
                                            a_label="REF (DNA)", b_label="EVO (DNA)", start_pos=start),
                       unsafe_allow_html=True
                   )
               else:
                   aa_map = {1: (14, 17), 2: (15, 18), 3: (16, 19), 4: (20, 23), 5: (21, 24), 6: (22, 25)}
                   r_row, e_row = aa_map.get(code, (14, 17))
                   ref_slice = "".join(info_matrix[r_row, start:end:3].tolist())
                   evo_slice = "".join(info_matrix[e_row, start:end:3].tolist())
                   st.subheader("Amino acid comparison")
                   st.markdown(
                       ms.highlighted_block(ref_slice, evo_slice,
                                            a_label="REF (AA)", b_label="EVO (AA)", start_pos=start),
                       unsafe_allow_html=True
                   )

               st.markdown("---")
               st.subheader(f"Mutation {selected_idx + 1} details")

               diffs_nt = ms.diff_indices(code, start, end, info_matrix)
               fs_row = info_matrix[5, :]
               dna_r, dna_e = (1, 3) if code in (1, 2, 3) else (10, 12)
               aa_map = {1: (14, 17), 2: (15, 18), 3: (16, 19), 4: (20, 23), 5: (21, 24), 6: (22, 25)}
               aa_r, aa_e = aa_map.get(code, (14, 17))

               nt_detail = []
               for gi in diffs_nt:
                   try:
                       ref_base = str(info_matrix[dna_r, gi])
                       evo_base = str(info_matrix[dna_e, gi])
                   except Exception:
                       ref_base = evo_base = "?"
                   try:
                       ref_aa = str(info_matrix[aa_r, gi])
                       evo_aa = str(info_matrix[aa_e, gi])
                   except Exception:
                       ref_aa = evo_aa = "?"
                   nt_detail.append({
                       "nucleotide index": int(gi),
                       "ref→evo (nt)": f"{ref_base} → {evo_base}",
                       "ref→evo (aa)": f"{ref_aa} → {evo_aa}",
                       "mutation type": ms.mut_type_from_frameshift(int(gi), fs_row),
                   })

               if nt_detail:
                   st.dataframe(
                       pd.DataFrame(nt_detail)
                         .sort_values("nucleotide index")
                         [["nucleotide index", "ref→evo (nt)", "ref→evo (aa)", "mutation type"]],
                       use_container_width=True,
                       hide_index=True,
                   )
               else:
                   st.info("No nucleotide mismatches detected within this span.")

               st.markdown("---")
               st.write("**Protein annotation (from BLAST):** ")
               blast_names = res.get("blast_names") or {}
               prot_name = blast_names.get(selected_idx)
               if prot_name:
                   safe_prot = html.escape(prot_name).replace("\n", "<br>")
                   st.markdown(
                       f"<div class='blastbox'><div class='blasttitle'>Mutation {selected_idx + 1}</div>"
                       f"<div class='blastbody'>{safe_prot}</div></div>",
                       unsafe_allow_html=True,
                   )
               else:
                   st.info("Submit to BLAST to see protein results.")
       else:
           st.info("Select a mutation above and click View to see details.")

with tab4:
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


with tab3:
   st.header("Selected mutations summary")
   res = st.session_state.get("results")
   if not res:
       st.info("Run the alignment in Tab 1 to populate this view.")
   else:
       rm = res.get("result_matrix")
       info_matrix = res.get("info_matrix")
       if rm is None or getattr(rm, "size", 0) == 0 or info_matrix is None:
           st.info("No mutation data available yet.")
       else:
           df = pd.DataFrame(rm, columns=["code", "start", "end"]).astype(int)
           selected_rows = sorted(
               i for i, checked in st.session_state.get("mut_checks", {}).items() if checked and i < len(df)
           )
           if not selected_rows:
               st.info("Select at least one mutation in Tab 2 to see it here.")
           else:
               st.session_state.setdefault("tab4_show_blast", False)
               st.session_state.setdefault("tab4_view_mode", "DNA sequence")
               st.session_state.setdefault("tab4_view_all", False)
               st.session_state.setdefault("tab4_pdf_mode", False)

               controls_left, controls_mid, controls_right = st.columns([0.33, 0.27, 0.4])
               with controls_left:
                   show_blast = st.checkbox("Show BLAST results", key="tab4_show_blast")
               with controls_mid:
                   pdf_mode = st.checkbox("PDF friendly layout", key="tab4_pdf_mode")
               mode_options = ["DNA sequence", "Amino acid sequence"]
               saved_mode = st.session_state.get("tab4_view_mode", mode_options[0])
               default_idx = mode_options.index(saved_mode) if saved_mode in mode_options else 0
               with controls_right:
                   mode = st.radio(
                       "Sequence display",
                       mode_options,
                       index=default_idx,
                       horizontal=True,
                       key="tab4_view_mode",
                   )

               action_col1, action_col2 = st.columns([0.2, 0.2])

               def _set_tab4_view_all(value: bool) -> None:
                   st.session_state["tab4_view_all"] = value

               if action_col1.button("View all", key="tab4_view_all_btn"):
                   _set_tab4_view_all(True)
               if action_col2.button("Clear view", key="tab4_clear_view_btn"):
                   _set_tab4_view_all(False)

               blast_names = res.get("blast_names") or {}

               wrap_width = 70 if pdf_mode else None  # tuned width for tighter PDF output

               st.markdown(
                   """
                   <style>
                   .pdf-block {
                       margin-bottom: 16px;
                       padding: 14px 16px;
                       border-radius: 12px;
                       border: 1px solid rgba(0, 0, 0, 0.12);
                       background: #fff;
                   }
                   .pdf-title {
                       font-weight: 600;
                       font-size: 18px;
                       color: #111;
                       margin-bottom: 10px;
                    }
                   .seqwrap.pdf {
                       border-color: rgba(0, 0, 0, 0.25);
                       background: #fafafa;
                       padding: 14px 16px;
                    }
                   .seqwrap.pdf .pdf-pair {
                       margin-bottom: 14px;
                    }
                   .seqwrap.pdf .pdf-pair:last-child {
                       margin-bottom: 0;
                    }
                   .seqwrap.pdf .pdf-pair .rowline {
                       margin-bottom: 2px;
                    }
                   .seqwrap.pdf .pdf-pair .rowline:last-child {
                       margin-bottom: 0;
                    }
                   @media print {
                       .pdf-block {
                           page-break-before: always;
                           margin: 0 0 8px 0;
                           padding: 10px 12px;
                       }
                       .pdf-block:first-child {
                           page-break-before: auto;
                       }
                   }
                   </style>
                   """,
                   unsafe_allow_html=True,
               )

               def _retrieve_sequences(code: int, start: int, end: int) -> tuple[str, str, str, str]:
                   if mode == "DNA sequence":
                       r_row, e_row = (1, 3) if code in (1, 2, 3) else (10, 12)
                       ref_slice = "".join(info_matrix[r_row, start:end].tolist())
                       evo_slice = "".join(info_matrix[e_row, start:end].tolist())
                       return ref_slice, evo_slice, "REF (DNA)", "EVO (DNA)"
                   aa_map = {
                       1: (14, 17),
                       2: (15, 18),
                       3: (16, 19),
                       4: (20, 23),
                       5: (21, 24),
                       6: (22, 25),
                   }
                   r_row, e_row = aa_map.get(code, (14, 17))
                   ref_slice = "".join(info_matrix[r_row, start:end:3].tolist())
                   evo_slice = "".join(info_matrix[e_row, start:end:3].tolist())
                   return ref_slice, evo_slice, "REF (AA)", "EVO (AA)"

               def _chunk_size(length: int) -> int:
                   return wrap_width if wrap_width else max(length, 1)

               def _build_highlight_html(ref_seq: str, evo_seq: str, ref_label: str, evo_label: str,
                                         start_pos: int, pdf: bool, chunk: int) -> str:
                   if not ref_seq or not evo_seq:
                       return "<div class='seqwrap pdf'>No sequence data available.</div>" if pdf else "<div class='seqwrap'>No sequence data available.</div>"
                   container_class = "seqwrap pdf" if pdf else "seqwrap"
                   blocks = []
                   for offset in range(0, len(ref_seq), chunk):
                       chunk_ref = ref_seq[offset:offset + chunk]
                       chunk_evo = evo_seq[offset:offset + chunk]
                       base_index = start_pos + offset
                       def _spanify(a: str, b: str) -> str:
                           cls = "mm" if b and a != b else "ok"
                           return f"<span class='{cls}'>{html.escape(a)}</span>"
                       la = "".join(_spanify(a, chunk_evo[i] if i < len(chunk_evo) else "") for i, a in enumerate(chunk_ref))
                       lb = "".join(_spanify(b, chunk_ref[i] if i < len(chunk_ref) else "") for i, b in enumerate(chunk_evo))
                       lab_a = f"{ref_label} {base_index:>7d}: "
                       lab_b = f"{evo_label} {base_index:>7d}: "
                       pair_wrapper_start = "<div class='pdf-pair'>" if pdf else ""
                       pair_wrapper_end = "</div>" if pdf else ""
                       blocks.append(
                           f"{pair_wrapper_start}"
                           f"<div class='rowline'><span class='lab'>{html.escape(lab_a)}</span><span class='diffblock'>{la}</span></div>"
                           f"<div class='rowline'><span class='lab'>{html.escape(lab_b)}</span><span class='diffblock'>{lb}</span></div>"
                           f"{pair_wrapper_end}"
                       )
                   legend = "<div class='legend'></div>"
                   return f"<div class='{container_class}'>{''.join(blocks)}</div>{legend}"

               def _build_copy_block(header: str, ref_label: str, evo_label: str,
                                     ref_seq: str, evo_seq: str, start_pos: int, chunk: int) -> str:
                   if not ref_seq or not evo_seq:
                       return f"{header}\n{ref_label}: (empty)\n{evo_label}: (empty)"
                   lines = [header]
                   for offset in range(0, len(ref_seq), chunk):
                       chunk_ref = ref_seq[offset:offset + chunk]
                       chunk_evo = evo_seq[offset:offset + chunk]
                       base_index = start_pos + offset
                       lines.append(f"{ref_label} {base_index:>7d}: {chunk_ref}")
                       lines.append(f"{evo_label} {base_index:>7d}: {chunk_evo}")
                   return "\n".join(lines)

               st.caption(f"{len(selected_rows)} mutation(s) selected for export")

               copy_blocks: list[str] = []
               render_blocks: list[str] = []
               for pos, idx in enumerate(selected_rows):
                   row = df.iloc[idx]
                   code = int(row["code"])
                   start = int(row["start"])
                   end = int(row["end"])
                   ref_seq, evo_seq, ref_label, evo_label = _retrieve_sequences(code, start, end)
                   chunk = _chunk_size(len(ref_seq))
                   header = f"Mutation {idx + 1} — code {code}, start {start}, end {end}"
                   highlight_html = _build_highlight_html(ref_seq, evo_seq, ref_label, evo_label, start, pdf_mode, chunk)

                   blast_section = ""
                   if show_blast:
                       prot_text = blast_names.get(idx)
                       if prot_text:
                           safe_prot = html.escape(str(prot_text)).replace("\n", "<br>")
                           blast_section = (
                               f"<div class='blastbox'><div class='blasttitle'>BLAST result</div>"
                               f"<div class='blastbody'>{safe_prot}</div></div>"
                           )
                       else:
                           blast_section = "<div class='blastbox'><div class='blasttitle'>BLAST result</div><div class='blastbody'>Submit to BLAST to see protein results.</div></div>"

                   style_attr = "page-break-before: always;" if pdf_mode and pos > 0 else ""
                   block_html = (
                       f"<div class='pdf-block' style=\"{style_attr}\">"
                       f"<div class='pdf-title'>{html.escape(header)}</div>"
                       f"{highlight_html}"
                       f"{blast_section if show_blast else ''}"
                       f"</div>"
                   )
                   render_blocks.append(block_html)
                   copy_text = _build_copy_block(header, ref_label, evo_label, ref_seq, evo_seq, start, chunk)
                   copy_blocks.append(copy_text)
               if copy_blocks:
                   combined_joiner = "\n\f\n" if pdf_mode else "\n\n"
                   combined_text = combined_joiner.join(copy_blocks)
                   copy_html = f"""
                   <div style="display:flex; justify-content:flex-end; margin: 8px 0 12px;">
                       <button id="tab4-copy-all"
                               style="background:#1b73e8; color:#fff; border:0; border-radius:6px;
                                      padding:6px 14px; font-weight:600; cursor:pointer;">
                           Copy all
                       </button>
                   </div>
                   <script>
                   (function() {{
                       const button = document.getElementById('tab4-copy-all');
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
                   components.html(copy_html, height=100)

               for pos, block_html in enumerate(render_blocks):
                   st.markdown(block_html, unsafe_allow_html=True)
                   if not pdf_mode and pos < len(render_blocks) - 1:
                       st.divider()


with tab5:
    st.header("Debug")
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

