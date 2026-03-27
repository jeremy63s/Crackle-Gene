import sys, os

print("Python executable:", sys.executable)
print("Current working dir:", os.getcwd())


import parasail
from sequence_utils import substitution_matrix
from tqdm import tqdm            # ← add this!
from typing import Tuple, List
# alignment.py
from typing import Callable, Optional, Any, Dict, List, Tuple

ProgressPairCB = Optional[Callable[[int, int], None]]  # k, total

def compare_sections(
    evo_sections: List[Tuple[str, int]],
    ref_sections: List[Tuple[str, int]],
    threshold: float,
    rotated_evolved: str,
    rotated_ref: str,
    progress_cb: ProgressPairCB = None,   # <-- NEW
) -> List[Dict[str, Any]]:
    matches: List[Dict[str, Any]] = []

    total = max(1, len(evo_sections) * len(ref_sections))
    k = 0

    # If you already use tqdm here, keep it. Just ALSO call progress_cb per step.
    for i, (evo_sec, evo_off) in enumerate(evo_sections, start=1):
        for j, (ref_sec, ref_off) in enumerate(ref_sections, start=1):
            # --- your existing per-pair logic ---
            # Example (adjust to your real code):
            # aligned_evo, aligned_ref, midline, align_score, bq, eq, br, er = get_alignment(evo_sec, ref_sec)
            # norm_score = compute_norm_score(...)
            # if norm_score >= threshold: matches.append({...})

            # --- tick progress every pair ---
            k += 1
            if progress_cb:
                progress_cb(k, total)

    return matches

def extract_extended_context(full_sequence, section_offset, alignment_start, alignment_end, extension=200):
    """
    Given the full rotated sequence, the section offset, and the alignment start/end indices
    (relative to the section), extract an extended context (extension bp before and after).
    The sequence is treated as circular.
    """
    L = len(full_sequence)
    abs_start = section_offset + alignment_start
    abs_end = section_offset + alignment_end
    start_context = (abs_start - extension) % L
    end_context = (abs_end + extension) % L
    if start_context < end_context:
        return full_sequence[start_context:end_context]
    else:
        return full_sequence[start_context:] + full_sequence[:end_context]



def get_alignment(evolved, reference):
    """
    Retrieve detailed alignment information for two sequences using Parasail's sw_trace.

    Returns:
      aligned_evo, aligned_ref, midline, alignment_score,
      begin_query, end_query, begin_ref, end_ref
    """
    trace_result = parasail.sw_trace(evolved, reference, 1, 1, substitution_matrix)
    tb = trace_result.traceback
    aligned_evo = tb.query
    aligned_ref = tb.ref
    midline = "".join('|' if a == b else ' ' for a, b in zip(aligned_evo, aligned_ref))
    non_gap_query = len(aligned_evo.replace("-", ""))
    non_gap_ref = len(aligned_ref.replace("-", ""))
    begin_query = trace_result.end_query - non_gap_query
    begin_ref = trace_result.end_ref - non_gap_ref

    return (aligned_evo, aligned_ref, midline, trace_result.score,
            begin_query, trace_result.end_query,
            begin_ref, trace_result.end_ref)

def diagonal_check(
    evolved_sections,
    reference_sections,
    threshold,
    progress_cb: Optional[Callable[[int, int], None]] = None,  # <-- NEW
):
    """
    Fast check: compare section i vs section i for all i,
    returns (pct_met, max_consec_failed, raw_scores, norm_scores).
    Emits progress via progress_cb(i, total) if provided.
    """
    n = min(len(evolved_sections), len(reference_sections))
    total = max(1, n)

    raw_scores: List[float] = []
    norm_scores: List[float] = []

    # tqdm in terminal + enumerate so we know the step index
    for i, ((evo_sec, _), (ref_sec, _)) in enumerate(
        tqdm(zip(evolved_sections, reference_sections),
             total=n, desc="Diagonal check", leave=False),
        start=1
    ):
        # compute once and store both raw and normalized
        raw = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix).score
        raw_scores.append(raw)
        norm_scores.append(raw / max(1, len(evo_sec)))

        # 🔑 tell Streamlit
        if progress_cb:
            progress_cb(i, total)

    # compute percent met
    meets = [s >= threshold for s in norm_scores]
    pct_met = (sum(meets) / total) * 100.0

    # compute longest consecutive failures
    max_run = 0
    run = 0
    for ok in meets:
        if not ok:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0

    return pct_met, max_run, raw_scores, norm_scores


def compare_sections(evolved_sections, reference_sections, threshold, full_evolved, full_reference):
    """
    Compare each evolved section (tuple: (section, offset)) to every reference section.

    For each pair, compute the raw alignment score using sw_stats and normalize it by dividing
    by the length of the evolved section. If normalized score ≥ threshold, retrieve detailed alignment info.
    Returns a list of match dictionaries. Also adds a 'match_index' field for later reference.
    """
    results = []
    match_counter = 0
    total_comparisons = len(evolved_sections) * len(reference_sections)
    pbar = tqdm(total=total_comparisons, desc="Comparing sections", leave=False)
    for i, (evo_sec, evo_offset) in enumerate(evolved_sections, start=1):
        for j, (ref_sec, ref_offset) in enumerate(reference_sections, start=1):
            stats_result = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix)
            raw_score = stats_result.score
            normalized_score = raw_score / len(evo_sec)
            if normalized_score >= threshold:
                (aligned_evo, aligned_ref, midline, align_score,
                 begin_query, end_query, begin_ref, end_ref) = get_alignment(evo_sec, ref_sec)
                extended_evo = extract_extended_context(full_evolved, evo_offset, begin_query, end_query, extension=200)
                extended_ref = extract_extended_context(full_reference, ref_offset, begin_ref, end_ref, extension=200)
                match = {
                    "evo_index": i,
                    "ref_index": j,
                    "raw_alignment_score": raw_score,
                    "normalized_score": normalized_score,
                    "alignment_score": align_score,
                    "aligned_evo": aligned_evo,
                    "aligned_ref": aligned_ref,
                    "midline": midline,
                    "begin_query": begin_query,
                    "end_query": end_query,
                    "begin_ref": begin_ref,
                    "end_ref": end_ref,
                    "evo_offset": evo_offset,
                    "ref_offset": ref_offset,
                    "extended_evo": extended_evo,
                    "extended_ref": extended_ref,
                    "is_match": 1,
                    "match_index": match_counter
                }
                results.append(match)
                match_counter += 1
            pbar.update(1)
    pbar.close()
    return results