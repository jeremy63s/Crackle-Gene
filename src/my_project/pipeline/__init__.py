# no print, input(), plt.show(), saving or writing global files, or tkinter
#Pure
##Sequence Utils##
import streamlit as st
import pandas as pd
import re
import io
import parasail
from tqdm import tqdm
CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",  # Start codon (also codes for M)
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*",  # Stop codons
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGA": "*",  # Stop codon
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}
def rotate_sequence(seq, start_subseq):
    """
    Rotate a circular sequence so that it begins with the specified starting subsequence.
    The sequence is treated as circular by doubling it.
    """
    seq_len = len(seq)
    sub_len = len(start_subseq)
    extended_seq = seq + seq  # simulate circularity

    for i in range(seq_len):
        if extended_seq[i:i + sub_len] == start_subseq:
            rotated = extended_seq[i:i + seq_len]
            # Ensure the rotated sequence begins exactly with start_subseq.
            if not rotated.startswith(start_subseq):
                idx = rotated.find(start_subseq)
                if idx != -1:
                    rotated = rotated[idx:] + rotated[:idx]
            return rotated
    raise ValueError("Starting subsequence not found in the sequence (even when treated as circular).")


def divide_sequence_with_offset(seq, n):
    """
    Divide the sequence into n sections and return a list of tuples (section, offset),
    where offset is the starting index of that section in the full sequence.
    """
    if n <= 0:
        raise ValueError("The number of sections must be positive.")
    L = len(seq)
    base = L // n
    sections = []
    for i in range(n - 1):
        start = i * base
        end = (i + 1) * base
        sections.append((seq[start:end], start))
    start = (n - 1) * base
    sections.append((seq[start:], start))
    return sections

#io helpers

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

#formerly part of get fasta
def parse_fasta_content(content: str) -> str:
    """
    Pure function to parse FASTA content and return the sequence as a single string.

    Args:
        content (str): The raw content of a FASTA file as a string

    Returns:
        str: The concatenated sequence data (without headers)

    Raises:
        ValueError: If the content cannot be parsed
    """
    try:
        seq_lines = []
        for line in content.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue  # Skip header lines
            seq_lines.append(line)
        return "".join(seq_lines)
    except Exception as e:
        raise ValueError("Error parsing FASTA content: " + str(e))


def _parse_fasta_simple(content: str) -> str:
    """
    Simple built-in FASTA parser (private function).

    Args:
        content (str): FASTA file content

    Returns:
        str: Concatenated sequence data
    """
    seq_lines = []
    for line in content.split('\n'):
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        seq_lines.append(line)
    return "".join(seq_lines)

def process_sequence_input(sequence_text: str, is_fasta_content: bool = False,
                           fasta_parser=None) -> str:
    """
    Pure function to process sequence input and return uppercase sequence.

    Args:
        sequence_text (str): Either raw sequence or FASTA file content
        is_fasta_content (bool): True if sequence_text is FASTA content, False for raw sequence
        fasta_parser (callable, optional): Function to parse FASTA content. If None and
                                         is_fasta_content=True, uses built-in parser.

    Returns:
        str: Uppercase sequence string

    Raises:
        ValueError: If sequence cannot be processed
    """
    try:
        if is_fasta_content:
            if fasta_parser:
                # Use provided FASTA parser function
                sequence = fasta_parser(sequence_text)
            else:
                # Built-in simple FASTA parser
                sequence = _parse_fasta_simple(sequence_text)
        else:
            # Direct sequence input
            sequence = sequence_text.strip()

        return sequence.upper()
    except Exception as e:
        raise ValueError(f"Error processing sequence: {str(e)}")

#sequence utils

def translate_nuc_row_active(nuc_row, codon_table):
    """
    Translates a nucleotide row into an amino acid row—but only within active coding regions.
    An active coding region is defined as a series of complete codons that starts with the start codon (ATG)
    and continues until a stop codon (TAA, TAG, or TGA) is reached. Codons outside these regions (or incomplete codons)
    are marked with '0'.

    Each amino acid translation is repeated in the three positions that made up its codon.

    Parameters:
      nuc_row (list): List or array of single-character nucleotide strings (e.g., ['A','T','G',...]).
      codon_table (dict): Dictionary mapping codon strings to their corresponding single-letter amino acid codes.

    Returns:
      list: A list the same length as nuc_row containing the amino acid letters for coding positions and '0' elsewhere.
    """
    # Initialize output with '0'
    output = ['0'] * len(nuc_row)
    # Collect indices for valid nucleotides (ignoring gaps '-' and ambiguous 'N')
    valid_indices = [i for i, nt in enumerate(nuc_row) if nt.upper() not in ['-', 'N']]

    inside_coding = False
    # Process the valid indices in groups of three (complete codons)
    for i in range(len(valid_indices) // 3):
        indices = valid_indices[3 * i: 3 * i + 3]
        codon = ''.join(nuc_row[idx] for idx in indices).upper()

        if not inside_coding:
            # Check if the codon is a start codon (ATG)
            if codon == "ATG":
                inside_coding = True
                aa = codon_table.get(codon, '0')
                for idx in indices:
                    output[idx] = aa
            # If not in coding region, leave these positions as '0'
        else:
            # Inside an active coding region: translate this codon
            aa = codon_table.get(codon, '0')
            for idx in indices:
                output[idx] = aa
            # If a stop codon is encountered, mark the end of the coding region after translating it
            if codon in ["TAA", "TAG", "TGA"]:
                inside_coding = False
    return output


def reverse_complement(nuc_row):
    """
    Computes the reverse complement of a nucleotide row while preserving gap positions.

    Parameters:
      nuc_row (list): List or array of nucleotide characters.

    Returns:
      list: The reverse complement of nuc_row.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-', 'N': 'N'}
    # Reverse the nucleotide row and replace each nucleotide with its complement
    rev_comp = [complement.get(nt.upper(), nt) for nt in nuc_row[::-1]]
    return rev_comp

# ============================
# [PURE] core computation
# ============================
from typing import List, Dict, Any, Sequence, Optional
import numpy as np

def compute_matrix_slice_display(
    data: np.ndarray,
    click_idx: int,
    window: int = 10,
    row_label_width: int = 10,
    cell_width: int = 12,
    highlight_rows: Sequence[int] = (2, 5, 6),
) -> Dict[str, Any]:
    """
    [PURE] Given a data matrix and a clicked column index in the ORIGINAL data,
    compute the slice, fixed-width text, and an HTML <pre> with:
      - highlighted rows (yellow)
      - highlighted center column (light blue)

    Returns a dict with:
      start_idx, end_idx, center_cell_index, subset, text_plain, html_pre
    """
    # bounds for exclusive end
    start_idx = max(0, click_idx - window)
    end_idx   = min(data.shape[1], click_idx + window + 1)
    subset = data[:, start_idx:end_idx]

    # which column in the subset is the clicked column
    center_cell_index = click_idx - start_idx

    # --- build fixed-width plain text
    row_lines: List[str] = []
    for i, row in enumerate(subset):
        row_label = f"Row {i}:".ljust(row_label_width)
        cells = [f"{str(cell):{cell_width}}" for cell in row]
        line = row_label + " " + " ".join(cells)
        row_lines.append(line)
    text_plain = "\n".join(row_lines) + "\n"

    # --- build HTML <pre> with highlights
    # We’ll reconstruct the same layout, but wrap:
    #  * entire highlighted rows with a <span style="background:yellow">
    #  * the center column range per row with <span style="background:lightblue">
    # NOTE: Do NOT insert overlapping spans; nest the column span inside the row span if needed.
    def html_escape(s: str) -> str:
        return (s.replace("&", "&amp;")
                 .replace("<", "&lt;")
                 .replace(">", "&gt;"))

    html_lines: List[str] = []
    # offsets
    start_offset = row_label_width + 1  # label + single space
    gap = 1                             # the single space between columns in the join

    for i, row in enumerate(subset):
        row_label = f"Row {i}:".ljust(row_label_width)
        cells = [f"{str(cell):{cell_width}}" for cell in row]
        # Precompute character positions of the center cell
        # row_label + " " + (cell0 + " " + cell1 + " " + ... )
        # col_start = start_offset + center_cell_index * (cell_width + gap)
        col_start = start_offset + center_cell_index * (cell_width + gap)
        col_end = col_start + cell_width

        # Build the one-line string first, then splice with spans.
        base_line = row_label + " " + " ".join(cells)
        esc = html_escape(base_line)

        # Safeguard against out-of-range (can happen on edges)
        col_start = max(0, min(col_start, len(esc)))
        col_end   = max(col_start, min(col_end, len(esc)))

        # Inject the center column span
        center_span = (
            esc[:col_start]
            + '<span style="background:lightblue;">'
            + esc[col_start:col_end]
            + '</span>'
            + esc[col_end:]
        )

        # If row is in highlight_rows, wrap whole line with yellow (including the blue span inside)
        if i in highlight_rows:
            center_span = f'<span style="background:yellow;">{center_span}</span>'

        html_lines.append(center_span)

    html_pre = (
        '<div style="font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, \'Liberation Mono\', monospace; '
        'font-size: 14px; line-height: 1.3;">'
        '<pre style="margin:0;">' + "\n".join(html_lines) + "</pre></div>"
    )

    return {
        "start_idx": start_idx,
        "end_idx": end_idx,
        "center_cell_index": center_cell_index,
        "subset": subset,
        "text_plain": text_plain,
        "html_pre": html_pre,
    }

def create_annotation_array(seq_len):
    """Return a list of length seq_len filled with '0'."""
    return ['0'] * seq_len

def array_to_str(arr):
    """Convert an array of single characters to a string."""
    return "".join(arr)

def parse_orf_frame(description):
    """
    From the description string (e.g.
    "ID=Seq_ORF.1;ORF_type=complete;ORF_len=9;ORF_frame=1;Start:ATG;Stop:TAG")
    extract and return the ORF frame as an integer (1, 2, or 3).
    """
    for part in description.split(';'):
        if part.startswith("ORF_frame="):
            try:
                return int(part.split('=')[1])
            except ValueError:
                return None
    return None


def find_boundaries_for_row(info_matrix, n, row):
    """
    For a given index n and a specified row (of the global info_matrix), this function
    scans outward—decrementing from n until a cell with '0' is encountered and incrementing
        until a cell with '0' is encountered.

        Boundaries are:
          left_boundary  = n - i + 1  (first index after the encountered '0' when scanning left)
          right_boundary = n + j - 1  (last index before the encountered '0' when scanning right)

        Parameters:
          n   (int): The starting index.
          row (int): The row number of info_matrix to search.

        Returns:
          A tuple (left_boundary, right_boundary).

        Raises:
          ValueError if n is out of bounds or if info_matrix[row, n] is '0'.
        """
        # info_matrix = final_matrix
    num_columns = info_matrix.shape[1]
    if n < 0 or n >= num_columns:
        raise ValueError(f"Input index {n} is out of bounds (0 to {num_columns - 1}).")

    if info_matrix[row, n] == '0':
        raise ValueError(f"info_matrix[{row}, {n}] is '0'. No protein found at the given index.")

        # Scan left.
    i = 0
    while (n - i) >= 0 and info_matrix[row, n - i] != '0':
            i += 1
    left_boundary = n - i + 1

    # Scan right.
    j = 0
    while (n + j) < num_columns and info_matrix[row, n + j] != '0':
        j += 1
    right_boundary = n + j - 1

    return left_boundary, right_boundary

from typing import Iterable, List, Tuple, Any, Dict, Callable, Optional

# =========================
# [PURE] core computation
# =========================
def process_orf_array_pure(
    info_matrix: Any,
    orf_array: Iterable[int],
    primary_row: int,
    secondary_row: int,
    orf_code: Any,
    *,
    find_boundaries_for_row_fn: Callable[[Any, int, int], Tuple[int, int]],
) -> Dict[str, Any]:
    """
    [PURE] Compute unique boundary blocks for mutation indices in `orf_array`.
    - Does NOT print or touch UI.
    - Returns both the results and a list of log messages describing what happened.

    Returns dict:
      {
        "output": List[List[Any]],           # rows: [orf_code, left_boundary, right_boundary]
        "logs": List[Dict[str, Any]]         # structured messages with 'level' and 'text'
      }
    """
    unique_boundaries = set()
    output: List[List[Any]] = []
    logs: List[Dict[str, Any]] = []

    for n in orf_array:
        try:
            boundaries_primary = find_boundaries_for_row_fn(info_matrix, int(n), primary_row)
            boundaries_secondary = find_boundaries_for_row_fn(info_matrix, int(n), secondary_row)
        except ValueError as err:
            logs.append({
                "level": "warning",
                "text": f"Skipping index {n} for ORF code {orf_code}: {err}"
            })
            continue

        if boundaries_primary != boundaries_secondary:
            logs.append({
                "level": "warning",
                "text": (f"Warning for ORF code {orf_code} at mutation index {n}: "
                         f"boundaries differ between rows {primary_row} and {secondary_row}. "
                         f"Using row {primary_row} boundaries.")
            })

        boundaries = boundaries_primary

        if boundaries in unique_boundaries:
            logs.append({
                "level": "info",
                "text": (f"For ORF code {orf_code} at mutation index {n}: "
                         f"mutation already accounted for")
            })
        else:
            unique_boundaries.add(boundaries)
            output.append([orf_code, boundaries[0], boundaries[1]])
            logs.append({
                "level": "info",
                "text": (f"For ORF code {orf_code} at mutation index {n}: "
                         f"boundaries: Left = {boundaries[0]}, Right = {boundaries[1]}")
            })

    return {"output": output, "logs": logs}
from typing import Dict, Iterable, List, Optional, Tuple, Any

# ---------- Pure core ----------

def extract_protein_seq(info_matrix, a: int, b: int, c: int) -> str:
    """
    PURE.
    Given info_matrix and slice indices a (row), b (start col), c (end col),
    return the protein sequence by taking every 3rd symbol starting at b.
    Assumes info_matrix supports numpy-like slicing: info_matrix[row, start:end].
    """
    region = info_matrix[a, b:c]
    # Taking every 3rd symbol starting at 0 within the region
    seq = region[0::3]
    # Ensure elements are strings before join (in case they're numpy scalars)
    return ''.join(map(str, seq))


def compute_orf_protein_seqs(
    info_matrix,
    result_matrix: Iterable[Tuple[Any, Any, Any]],
    orf_map: Optional[Dict[int, int]] = None
) -> List[Dict[str, Any]]:
    """
    PURE.
    For each row [orf_code, left, right] in result_matrix:
      - map orf_code → a row index using orf_map
      - extract protein_seq
    Returns a list of dicts, one per input row, including any structured errors.
    No prints, no I/O, no BLAST calls.
    """
    if orf_map is None:
        orf_map = {1: 14, 2: 15, 3: 16, 4: 20, 5: 21, 6: 22}

    out: List[Dict[str, Any]] = []

    for idx, (orf_code_raw, left_raw, right_raw) in enumerate(result_matrix):
        item = {
            "index": idx,
            "orf_code_raw": orf_code_raw,
            "left_raw": left_raw,
            "right_raw": right_raw,
            "error": None,
            "orf_code": None,
            "matrix_row": None,
            "start": None,
            "end": None,
            "protein_seq": None,
        }

        try:
            orf_code = int(orf_code_raw)
        except (ValueError, TypeError):
            item["error"] = f"Non-integer ORF code: {orf_code_raw!r}"
            out.append(item)
            continue

        if orf_code not in orf_map:
            item["error"] = f"Unknown ORF code {orf_code}"
            out.append(item)
            continue

        try:
            start = int(left_raw)
            end = int(right_raw)
        except (ValueError, TypeError):
            item["error"] = f"Invalid start/end: {left_raw!r}, {right_raw!r}"
            out.append(item)
            continue

        a_row = orf_map[orf_code]

        try:
            protein_seq = extract_protein_seq(info_matrix, a_row, start, end)
        except Exception as e:
            item["error"] = f"Extraction failed: {e}"
            out.append(item)
            continue

        item.update(
            orf_code=orf_code,
            matrix_row=a_row,
            start=start,
            end=end,
            protein_seq=protein_seq,
        )
        out.append(item)

    return out

from typing import Any, Dict, List, Optional, Tuple
from io import StringIO

# ---------- Pure core (no I/O) ----------

def fasta_from_seq(seq: str, header: str = "query") -> str:
    """
    PURE.
    Build a FASTA string from a raw amino-acid sequence.
    """
    # Defensive cleanup: strip whitespace/newlines
    seq = "".join(str(seq).split())
    return f">{header}\n{seq}\n"


def parse_blast_xml(xml_text: str, top_k: Optional[int] = None) -> Dict[str, Any]:
    """
    PURE.
    Parse BLAST XML into a light-weight summary structure.
    Returns a dict with a list of 'hits' and basic metadata.
    """
    from Bio.Blast import NCBIXML  # import here to keep core pure of I/O
    record = NCBIXML.read(StringIO(xml_text))

    hits: List[Dict[str, Any]] = []
    alignments = getattr(record, "alignments", []) or []
    k = len(alignments) if top_k is None else min(top_k, len(alignments))

    for aln in alignments[:k]:
        # Take the best HSP (high-scoring pair) if present
        hsp = aln.hsps[0] if getattr(aln, "hsps", []) else None
        hits.append({
            "accession": getattr(aln, "accession", ""),
            "definition": getattr(aln, "hit_def", ""),
            "length": getattr(aln, "length", None),
            "e_value": getattr(hsp, "expect", None) if hsp else None,
            "bit_score": getattr(hsp, "bits", None) if hsp else None,
            "identities": getattr(hsp, "identities", None) if hsp else None,
            "align_len": getattr(hsp, "align_length", None) if hsp else None,
        })

    return {
        "program": getattr(record, "program", "blastp"),
        "database": getattr(record, "database", "nr"),
        "query_length": getattr(record, "query_length", None),
        "hits": hits,
        # Optionally also return the raw record if you want parity with the original API:
        "record": record,
    }

def compute_indices(aligned_seq, offset, begin):
    """
    Given an aligned sequence (with gaps), section offset, and beginning index,
    compute an array of absolute indices (starting at 1) for each column; gaps get 0.
    """
    indices = []
    current_index = offset + begin + 1
    for char in aligned_seq:
        if char == '-':
            indices.append(0)
        else:
            indices.append(current_index)
            current_index += 1
    return indices

def compute_indicator(aligned_evo, aligned_ref):
    """
    For each aligned column, compute an indicator:
      1 if both non-gap and identical,
     -1 if both non-gap but different,
      0 if either is a gap.
    Returns two identical lists.
    """
    indicator = []
    for a, r in zip(aligned_evo, aligned_ref):
        if a != '-' and r != '-':
            indicator.append(1 if a == r else -1)
        else:
            indicator.append(0)
    return indicator #, indicator

def build_chunk_matrix(match):
    """
    Given a match dictionary, build a 7 x (aligned_length) NumPy array:
      Row 1: Info ("R:<ref_index>, E:<evo_index>, Score:<normalized_score>") repeated across columns.
      Row 2: Reference nucleotides (aligned)
      Row 3: Reference indices (absolute, starting at 1; gaps=0)
      Row 4: Evolved nucleotides (aligned)
      Row 5: Evolved indices (absolute, starting at 1; gaps=0)
      Row 6: Evolved indicator (1 if match, -1 if substitution, 0 if gap)
      Row 7: Reference indicator (same as evolved)
    """
    aligned_evo = match["aligned_evo"]
    aligned_ref = match["aligned_ref"]
    N = len(aligned_evo)
    info_str = f"R:{match['ref_index']}, E:{match['evo_index']}, Score:{match['normalized_score']:.3f}"
    row1 = [info_str] * N
    row2 = list(aligned_ref)
    row3 = compute_indices(aligned_ref, match["ref_offset"], match["begin_ref"])
    row4 = list(aligned_evo)
    row5 = compute_indices(aligned_evo, match["evo_offset"], match["begin_query"])
    ind_evo, ind_ref = compute_indicator(aligned_evo, aligned_ref)
    row6 = ind_evo
    row7 = ind_ref
    chunk_matrix = np.array([row1, row2, row3, row4, row5, row6, row7], dtype=object)
    return chunk_matrix

def build_final_matrix(chunk_matrices):
    """
    Horizontally concatenate all chunk matrices (each 7 x aligned_length) into one big matrix,
    then append an 8th row that is a counter for each column.
    Returns an 8 x (total_columns) matrix.
    """
    if not chunk_matrices:
        return None
    concatenated = np.hstack(chunk_matrices)  # shape: (7, total_cols)
    total_cols = concatenated.shape[1]
    counter_row = np.arange(1, total_cols +  1).reshape(1, -1)
    final_matrix = np.vstack([concatenated, counter_row])
    #print("concatenated shape is", concatenated.shape)
    #print("counter_row shape is", counter_row.shape)
    return final_matrix


import parasail


def get_alignment_pure(evolved, reference, substitution_matrix):
    """
    Pure function version of get_alignment for Streamlit compatibility.

    Retrieve detailed alignment information for two sequences using Parasail's sw_trace.

    Args:
        evolved (str): Query sequence
        reference (str): Reference sequence
        substitution_matrix: Parasail substitution matrix object

    Returns:
        tuple: (aligned_evo, aligned_ref, midline, alignment_score,
                begin_query, end_query, begin_ref, end_ref)
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

    return (
        aligned_evo,
        aligned_ref,
        midline,
        trace_result.score,
        begin_query,
        trace_result.end_query,
        begin_ref,
        trace_result.end_ref
    )


def diagonal_check_pure(evolved_sections, reference_sections, threshold, substitution_matrix, progress_callback=None):
    """
    Pure function version of diagonal_check for Streamlit compatibility.

    Fast check: compare section i vs section i for all i,
    returns (pct_met, max_consec_failed, raw_scores, norm_scores).

    Args:
        evolved_sections (list): List of evolved sequence sections
        reference_sections (list): List of reference sequence sections
        threshold (float): Threshold for normalized scores
        substitution_matrix: Parasail substitution matrix object
        progress_callback (callable, optional): Function to call for progress updates

    Returns:
        tuple: (pct_met, max_consec_failed, raw_scores, norm_scores)
    """
    n = len(evolved_sections)
    raw_scores = []
    norm_scores = []

    for i, ((evo_sec, *_), (ref_sec, *_)) in enumerate(zip(evolved_sections, reference_sections)):
        # Update progress if callback provided
        if progress_callback:
            progress_callback(i / n, f"Processing section {i + 1}/{n}")

        # Compute once and store both raw and normalized
        raw = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix).score
        raw_scores.append(raw)
        norm_scores.append(raw / len(evo_sec))

    # Final progress update
    if progress_callback:
        progress_callback(1.0, "Complete")

    # Compute percent met
    meets = [s >= threshold for s in norm_scores]
    pct_met = sum(meets) / n * 100

    # Compute longest consecutive failures
    max_run = 0
    run = 0
    for ok in meets:
        if not ok:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0

    return pct_met, max_run, raw_scores, norm_scores


def compare_sections_pure(evolved_sections, reference_sections, threshold,
                          full_evolved, full_reference, substitution_matrix,
                          get_alignment_func, extract_extended_context_func,
                          progress_callback=None):
    """
    Pure function version of compare_sections for Streamlit compatibility.

    Compare each evolved section (tuple: (section, offset)) to every reference section.
    For each pair, compute the raw alignment score using sw_stats and normalize it by dividing
    by the length of the evolved section. If normalized score ≥ threshold, retrieve detailed alignment info.

    Args:
        evolved_sections (list): List of (section, offset) tuples for evolved sequences
        reference_sections (list): List of (section, offset) tuples for reference sequences
        threshold (float): Threshold for normalized scores
        full_evolved (str): Complete evolved sequence
        full_reference (str): Complete reference sequence
        substitution_matrix: Parasail substitution matrix object
        get_alignment_func (callable): Function to get detailed alignment
        extract_extended_context_func (callable): Function to extract extended context
        progress_callback (callable, optional): Function to call for progress updates

    Returns:
        list: List of match dictionaries with detailed alignment information
    """
    results = []
    match_counter = 0
    total_comparisons = len(evolved_sections) * len(reference_sections)

    comparison_count = 0
    for i, (evo_sec, evo_offset) in enumerate(evolved_sections, start=1):
        for j, (ref_sec, ref_offset) in enumerate(reference_sections, start=1):
            # Update progress if callback provided
            if progress_callback:
                progress = comparison_count / total_comparisons
                message = f"Comparing section {i} vs {j} ({comparison_count + 1}/{total_comparisons})"
                progress_callback(progress, message)

            stats_result = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix)
            raw_score = stats_result.score
            normalized_score = raw_score / len(evo_sec)

            if normalized_score >= threshold:
                (aligned_evo, aligned_ref, midline, align_score,
                 begin_query, end_query, begin_ref, end_ref) = get_alignment_func(evo_sec, ref_sec)

                extended_evo = extract_extended_context_func(
                    full_evolved, evo_offset, begin_query, end_query, extension=200
                )
                extended_ref = extract_extended_context_func(
                    full_reference, ref_offset, begin_ref, end_ref, extension=200
                )

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

            comparison_count += 1

    # Final progress update
    if progress_callback:
        progress_callback(1.0, f"Complete: Found {len(results)} matches")

    return results









# ---------- Streamlit adapter (impure by design: UI + optional BLAST) ----------



#####FINISHED ADAPTERS#######
#sequence_utils#

def read_fasta(filepath: str) -> str:
    """
    Adapter function to read a FASTA file and return the sequence as a single string.
    This function handles file I/O and delegates parsing to the pure function.

    Args:
        filepath (str): Path to the FASTA file

    Returns:
        str: The concatenated sequence data

    Raises:
        ValueError: If the file cannot be read or parsed
    """
    try:
        with open(filepath, 'r') as file:
            content = file.read()
        return parse_fasta_content(content)
    except Exception as e:
        raise ValueError("Error reading FASTA file: " + str(e))


def read_fasta_from_upload(uploaded_file) -> str:
    """
    Streamlit-specific adapter function to read FASTA content from an uploaded file.

    Args:
        uploaded_file: Streamlit UploadedFile object

    Returns:
        str: The concatenated sequence data

    Raises:
        ValueError: If the file content cannot be parsed
    """
    try:
        # Convert bytes to string if necessary
        if hasattr(uploaded_file, 'read'):
            content = uploaded_file.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
        else:
            content = str(uploaded_file)

        return parse_fasta_content(content)
    except Exception as e:
        raise ValueError("Error reading uploaded FASTA file: " + str(e))


def get_sequence_cli(seq_type: str, fasta_reader_func=None) -> str:
    """
    CLI adapter function to get sequence via command line input.

    Args:
        seq_type (str): Type of sequence (e.g., "DNA", "protein")
        fasta_reader_func: Function to read FASTA files (e.g., read_fasta)

    Returns:
        str: Uppercase sequence string
    """
    method = input(
        f"Do you want to input the {seq_type} sequence manually or via FASTA file? (M for manual, F for FASTA): "
    ).strip().upper()

    if method == 'F':
        filepath = input(f"Enter the path to the {seq_type} FASTA file: ").strip()
        if fasta_reader_func:
            sequence_content = fasta_reader_func(filepath)
        else:
            # Fallback to direct file reading
            with open(filepath, 'r') as f:
                sequence_content = f.read()
        return process_sequence_input(sequence_content, is_fasta_content=True)
    else:
        sequence_text = input(f"Enter the {seq_type} sequence: ").strip()
        return process_sequence_input(sequence_text, is_fasta_content=False)


def get_sequence_streamlit(seq_type: str, manual_input: str = None,
                           uploaded_file=None, fasta_content: str = None) -> str:
    """
    Streamlit adapter function to get sequence from UI components.

    Args:
        seq_type (str): Type of sequence (for error messages)
        manual_input (str, optional): Manually entered sequence
        uploaded_file (optional): Streamlit uploaded file object
        fasta_content (str, optional): Pre-read FASTA content

    Returns:
        str: Uppercase sequence string

    Raises:
        ValueError: If no valid input is provided
    """
    if uploaded_file is not None:
        # Handle uploaded file
        try:
            content = uploaded_file.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            return process_sequence_input(content, is_fasta_content=True)
        except Exception as e:
            raise ValueError(f"Error reading uploaded {seq_type} file: {str(e)}")

    elif fasta_content:
        # Handle pre-read FASTA content
        return process_sequence_input(fasta_content, is_fasta_content=True)

    elif manual_input:
        # Handle manual input
        return process_sequence_input(manual_input, is_fasta_content=False)

    else:
        raise ValueError(f"No {seq_type} sequence input provided")

def get_sequence_programmatic(sequence_data: str, input_type: str = "raw") -> str:
        """
        Programmatic adapter for direct function calls.

        Args:
            sequence_data (str): The sequence data (raw sequence or FASTA content)
            input_type (str): Either "raw" for direct sequence or "fasta" for FASTA content

        Returns:
            str: Uppercase sequence string
        """
        is_fasta = input_type.lower() == "fasta"
        return process_sequence_input(sequence_data, is_fasta_content=is_fasta)

# ============================================
# [ADAPTER] Matplotlib pick -> pure core call
# ============================================
def on_pick_mpl_adapter(event, scatter, plotted_indices, data: np.ndarray, *, window=10):
    """
    [ADAPTER] Minimal adapter that converts a Matplotlib pick event into a call
    to the pure function. Returns the computed payload; rendering is left to UI adapters.
    """
    if event.artist is not scatter:
        return None
    if not getattr(event, "ind", None):
        return None

    picked_index = event.ind[0]
    click_idx = plotted_indices[picked_index]

    return compute_matrix_slice_display(
        data=data,
        click_idx=int(click_idx),
        window=window,
        row_label_width=10,
        cell_width=12,
        highlight_rows=(2, 5, 6),
    )


# ==================================
# [ADAPTER] Streamlit render adapter
# ==================================
def render_slice_streamlit(payload, st, title_prefix: str = "Data Matrix Slice"):
    """
    [ADAPTER] Streamlit renderer for the payload produced by compute_matrix_slice_display.
    Uses HTML <pre> for highlighting.
    """
    if payload is None:
        return
    start_idx = payload["start_idx"]
    end_idx = payload["end_idx"]
    st.subheader(f"{title_prefix} (Columns {start_idx} to {end_idx - 1})")
    st.markdown(payload["html_pre"], unsafe_allow_html=True)
    # If you also want a raw copyable text block:
    with st.expander("Show plain text"):
        st.text(payload["text_plain"])

# ==================================
# [OPTIONAL ADAPTER] Tk renderer
# (Only if you still need desktop Tk)
# ==================================
import tkinter as tk
from tkinter import scrolledtext

def render_slice_tk(payload, root: Optional[tk.Tk] = None, title_prefix="Data Matrix Slice"):
    """
    [ADAPTER] Tkinter renderer. Behavior matches your original UI, but without
    any computation inside.
    """
    if payload is None:
        return

    if root is None:
        root = tk._default_root or tk.Tk()

    start_idx = payload["start_idx"]
    end_idx   = payload["end_idx"]
    center_cell_index = payload["center_cell_index"]
    subset = payload["subset"]
    text_plain = payload["text_plain"]

    row_label_width = 10
    cell_width = 12

    win = tk.Toplevel(root)
    win.title(f"{title_prefix} (Columns {start_idx} to {end_idx-1})")

    text_widget = scrolledtext.ScrolledText(win, wrap='none', font=("Courier", 14))
    text_widget.insert(tk.END, text_plain)
    text_widget.configure(state='normal')

    # highlight rows 2,5,6 in yellow
    highlight_rows = [2, 5, 6]
    for row_idx in highlight_rows:
        if row_idx < subset.shape[0]:
            line_no = row_idx + 1
            tag = f"highlight_row_{row_idx}"
            text_widget.tag_add(tag, f"{line_no}.0", f"{line_no}.end")
            text_widget.tag_config(tag, background="yellow")

    # highlight center column (light blue)
    start_offset = row_label_width + 1
    for i in range(subset.shape[0]):
        line_no = i + 1
        col_start = start_offset + center_cell_index * (cell_width + 1)
        col_end = col_start + cell_width
        text_widget.tag_add("center_col", f"{line_no}.{col_start}", f"{line_no}.{col_end}")
    text_widget.tag_config("center_col", background="lightblue")

    text_widget.configure(state='disabled')
    text_widget.pack(expand=True, fill='both')

# ==========================================
# [ADAPTER] Console/terminal print adapter
# ==========================================
def process_orf_array_console_adapter(
    info_matrix: Any,
    orf_array: Iterable[int],
    primary_row: int,
    secondary_row: int,
    orf_code: Any,
    *,
    find_boundaries_for_row_fn: Callable[[Any, int, int], Tuple[int, int]],
) -> List[List[Any]]:
    """
    [ADAPTER] Preserves original side-effect behavior by printing logs to stdout.
    """
    result = process_orf_array_pure(
        info_matrix, orf_array, primary_row, secondary_row, orf_code,
        find_boundaries_for_row_fn=find_boundaries_for_row_fn
    )
    for msg in result["logs"]:
        print(msg["text"])
    return result["output"]

# ==========================================
# [ADAPTER] Streamlit UI adapter
# ==========================================
def process_orf_array_streamlit_adapter(
    info_matrix: Any,
    orf_array: Iterable[int],
    primary_row: int,
    secondary_row: int,
    orf_code: Any,
    *,
    st,  # pass the streamlit module (import streamlit as st)
    find_boundaries_for_row_fn: Callable[[Any, int, int], Tuple[int, int]],
) -> List[List[Any]]:
    """
    [ADAPTER] Renders logs via Streamlit and returns the output table.
    """
    result = process_orf_array_pure(
        info_matrix, orf_array, primary_row, secondary_row, orf_code,
        find_boundaries_for_row_fn=find_boundaries_for_row_fn
    )

    # Show logs with appropriate levels
    for msg in result["logs"]:
        lvl = msg.get("level", "info")
        if   lvl == "warning": st.warning(msg["text"])
        elif lvl == "error":   st.error(msg["text"])
        elif lvl == "success": st.success(msg["text"])
        else:                  st.info(msg["text"])

    # Optionally show a table of unique blocks
    if result["output"]:
        st.write("**Unique ORF boundary blocks**")
        st.dataframe(result["output"], use_container_width=True)

    return result["output"]

def process_result_matrix_streamlit(
    info_matrix,
    result_matrix: Iterable[Tuple[Any, Any, Any]],
    do_blast: bool = True,
    orf_map: Optional[Dict[int, int]] = None,
    blast_top_k: int = 5,
):
    """
    IMPURE: renders to Streamlit and may call BLAST.
    - Uses the pure core to compute sequences
    - Optionally calls blast_protein_seq(protein_seq) for each non-empty sequence
      and shows the top K hits.

    Requirements:
      - `streamlit` must be installed
      - You must provide a `blast_protein_seq(seq: str)` function that returns an
        object with `.alignments`, where each alignment has `.accession`, `.hit_def`,
        and `.hsps[0].expect` (E-value). Adjust parsing if your BLAST wrapper differs.
    """


    results = compute_orf_protein_seqs(info_matrix, result_matrix, orf_map)

    st.subheader("Processed ORFs")
    for r in results:
        idx = r["index"]
        if r["error"]:
            st.warning(f"Row {idx}: {r['error']}")
            continue

        st.markdown(
            f"**Row {idx}** — ORF `{r['orf_code']}` → matrix row `{r['matrix_row']}`; "
            f"cols `{r['start']}`–`{r['end']}`"
        )
        protein_seq = r["protein_seq"] or ""
        st.code(protein_seq if protein_seq else "(empty sequence)", language="text")

        if do_blast:
            if protein_seq:
                with st.spinner("Submitting to BLAST…"):
                    try:
                        record = blast_protein_seq(protein_seq)  # you provide this
                        # Render top-K hits (defensively handle short lists)
                        hits = []
                        for aln in (record.alignments[:blast_top_k] if getattr(record, "alignments", None) else []):
                            try:
                                hsp = aln.hsps[0]
                                hits.append({
                                    "Accession": getattr(aln, "accession", ""),
                                    "Definition": getattr(aln, "hit_def", ""),
                                    "E-value": getattr(hsp, "expect", None),
                                })
                            except Exception:
                                continue
                        if hits:
                            import pandas as pd
                            st.dataframe(pd.DataFrame(hits))
                        else:
                            st.info("No BLAST hits found or record format unrecognized.")
                    except Exception as e:
                        st.error(f"BLAST failed: {e}")
            else:
                st.info("Empty sequence — skipping BLAST.")

###############################
##### IMPURE WORKING SPACE#####
###############################

# old


def blast_protein_seq(seq: str, email: Optional[str] = None, **qblast_kwargs):
    """
    IMPURE.
    Submit a BLASTP query to NCBI and return the parsed Record object
    (backwards-compatible with your original function).

    Args:
      seq: amino-acid sequence
      email: (recommended) contact email for NCBI usage policy
      **qblast_kwargs: forwarded to Bio.Blast.NCBIWWW.qblast (e.g., expect, hitlist_size)

    Returns:
      Bio.Blast.Record.Blast instance (same as NCBIXML.read(handle))
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    fasta = fasta_from_seq(seq, header="query")
    # Respect NCBI; supply email if given, and allow caller overrides
    kwargs = dict(program="blastp", database="nr", sequence=fasta, format_type="XML")
    kwargs.update(qblast_kwargs)

    # NCBIWWW.qblast signature typically: qblast(program, database, sequence, **kwargs)
    handle = NCBIWWW.qblast(kwargs.pop("program"), kwargs.pop("database"),
                            kwargs.pop("sequence"), **kwargs)
    # Read the XML once; then parse into a Record object (pure step)
    xml_text = handle.read()
    handle.close()
    record = NCBIXML.read(StringIO(xml_text))
    return record


# ---------- Streamlit adapter (UI + BLAST) ----------

def blast_protein_seq_streamlit(seq: str, email: Optional[str] = None, top_k: int = 5, **qblast_kwargs):
    """
    IMPURE (UI + network).
    Runs BLAST via NCBI, then displays a small table of the top-K hits in Streamlit.
    Returns the full Biopython Record object.
    """


    st.subheader("BLASTP Query (NCBI)")
    st.code(fasta_from_seq(seq), language="text")

    with st.spinner("Submitting to NCBI BLAST…"):
        record = blast_protein_seq(seq, email=email, **qblast_kwargs)

    # Convert the record to a light summary for display (pure parsing step)
    from Bio.Blast import NCBIXML
    # Re-serialize by walking the Record -> simple hits (avoid re-reading XML)
    # We'll mimic parse_blast_xml's hit extraction directly from 'record':
    hits = []
    alignments = getattr(record, "alignments", []) or []
    for aln in alignments[:top_k]:
        hsp = aln.hsps[0] if getattr(aln, "hsps", []) else None
        hits.append({
            "Accession": getattr(aln, "accession", ""),
            "Definition": getattr(aln, "hit_def", ""),
            "Length": getattr(aln, "length", None),
            "E-value": getattr(hsp, "expect", None) if hsp else None,
            "Bit score": getattr(hsp, "bits", None) if hsp else None,
            "Identities": getattr(hsp, "identities", None) if hsp else None,
            "Align len": getattr(hsp, "align_length", None) if hsp else None,
        })

    if hits:
        st.dataframe(pd.DataFrame(hits))
    else:
        st.info("No hits found or unexpected result format.")

    return record


def get_alignment(evolved, reference, substitution_matrix=None):
    """
    Adaptor function to maintain original interface.

    If substitution_matrix is not provided, creates a default one.
    This preserves the original function's behavior while making it compatible
    with the pure function approach.

    Args:
        evolved (str): Query sequence
        reference (str): Reference sequence
        substitution_matrix: Optional substitution matrix (will create default if None)

    Returns:
        tuple: Same as get_alignment_pure
    """
    if substitution_matrix is None:
        # Create a default substitution matrix if none provided
        substitution_matrix = parasail.blosum62  # or parasail.nuc44 for nucleotides

    return get_alignment_pure(evolved, reference, substitution_matrix)


# STREAMLIT USAGE EXAMPLE
def streamlit_alignment_interface():
    """
    Example of how to use this in a Streamlit app.
    This would be called from your main Streamlit script.
    """
    import streamlit as st

    st.title("Sequence Alignment Tool")

    # User inputs
    evolved_seq = st.text_area("Evolved Sequence", height=100)
    reference_seq = st.text_area("Reference Sequence", height=100)

    # Matrix selection
    matrix_option = st.selectbox(
        "Substitution Matrix",
        ["blosum62", "nuc44", "pam250"]
    )

    # Map selection to parasail matrix
    matrix_map = {
        "blosum62": parasail.blosum62,
        "nuc44": parasail.nuc44,
        "pam250": parasail.pam250
    }

    if st.button("Perform Alignment") and evolved_seq and reference_seq:
        try:
            # Call the pure function
            results = get_alignment_pure(
                evolved_seq.strip(),
                reference_seq.strip(),
                matrix_map[matrix_option]
            )

            aligned_evo, aligned_ref, midline, score, begin_q, end_q, begin_r, end_r = results

            # Display results
            st.subheader("Alignment Results")
            st.write(f"**Alignment Score:** {score}")
            st.write(f"**Query Position:** {begin_q} - {end_q}")
            st.write(f"**Reference Position:** {begin_r} - {end_r}")

            # Show alignment in monospace
            st.text("Alignment:")
            st.code(f"Query:     {aligned_evo}\nMidline:   {midline}\nReference: {aligned_ref}")

        except Exception as e:
            st.error(f"Error performing alignment: {str(e)}")


def diagonal_check(evolved_sections, reference_sections, threshold, substitution_matrix=None):
    """
    Adaptor function to maintain original interface with tqdm progress bar.

    Args:
        evolved_sections (list): List of evolved sequence sections
        reference_sections (list): List of reference sequence sections
        threshold (float): Threshold for normalized scores
        substitution_matrix: Optional substitution matrix (will create default if None)

    Returns:
        tuple: Same as diagonal_check_pure
    """
    if substitution_matrix is None:
        # Create a default substitution matrix if none provided
        substitution_matrix = parasail.blosum62  # or parasail.nuc44 for nucleotides

    n = len(evolved_sections)
    raw_scores = []
    norm_scores = []

    # Use tqdm for console progress (original behavior)
    for (evo_sec, *_), (ref_sec, *_) in tqdm(
            zip(evolved_sections, reference_sections),
            total=n,
            desc="Diagonal check",
            leave=False
    ):
        # Compute once and store both raw and normalized
        raw = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix).score
        raw_scores.append(raw)
        norm_scores.append(raw / len(evo_sec))

    # Compute percent met
    meets = [s >= threshold for s in norm_scores]
    pct_met = sum(meets) / n * 100

    # Compute longest consecutive failures
    max_run = 0
    run = 0
    for ok in meets:
        if not ok:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0

    return pct_met, max_run, raw_scores, norm_scores


# STREAMLIT USAGE EXAMPLE
def streamlit_diagonal_check_interface():
    """
    Example of how to use this in a Streamlit app.
    This would be called from your main Streamlit script.
    """
    import streamlit as st
    import pandas as pd
    import plotly.express as px

    st.title("Diagonal Check Analysis")

    # File upload or manual input
    upload_option = st.radio("Data Input Method", ["Upload Files", "Manual Input"])

    evolved_sections = []
    reference_sections = []

    if upload_option == "Upload Files":
        evo_file = st.file_uploader("Upload Evolved Sections", type=['txt', 'fasta'])
        ref_file = st.file_uploader("Upload Reference Sections", type=['txt', 'fasta'])

        if evo_file and ref_file:
            # Parse files (implementation depends on your file format)
            evolved_sections = parse_sections(evo_file)
            reference_sections = parse_sections(ref_file)
    else:
        # Manual input for testing
        num_sections = st.number_input("Number of sections", min_value=1, value=5)
        for i in range(num_sections):
            col1, col2 = st.columns(2)
            with col1:
                evo_seq = st.text_input(f"Evolved Section {i + 1}", key=f"evo_{i}")
                evolved_sections.append((evo_seq,))
            with col2:
                ref_seq = st.text_input(f"Reference Section {i + 1}", key=f"ref_{i}")
                reference_sections.append((ref_seq,))

    # Parameters
    threshold = st.slider("Threshold", min_value=0.0, max_value=100.0, value=10.0, step=0.1)

    matrix_option = st.selectbox(
        "Substitution Matrix",
        ["blosum62", "nuc44", "pam250"]
    )

    matrix_map = {
        "blosum62": parasail.blosum62,
        "nuc44": parasail.nuc44,
        "pam250": parasail.pam250
    }

    if st.button("Run Diagonal Check") and evolved_sections and reference_sections:
        if len(evolved_sections) != len(reference_sections):
            st.error("Number of evolved and reference sections must match!")
            return

        try:
            # Create progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()

            def update_progress(progress, message):
                progress_bar.progress(progress)
                status_text.text(message)

            # Call the pure function with progress callback
            pct_met, max_run, raw_scores, norm_scores = diagonal_check_pure(
                evolved_sections,
                reference_sections,
                threshold,
                matrix_map[matrix_option],
                progress_callback=update_progress
            )

            # Clear progress indicators
            progress_bar.empty()
            status_text.empty()

            # Display results
            st.subheader("Results Summary")
            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric("Sections Meeting Threshold", f"{pct_met:.1f}%")
            with col2:
                st.metric("Max Consecutive Failures", max_run)
            with col3:
                st.metric("Total Sections", len(evolved_sections))

            # Create detailed results DataFrame
            results_df = pd.DataFrame({
                'Section': range(1, len(raw_scores) + 1),
                'Raw Score': raw_scores,
                'Normalized Score': norm_scores,
                'Meets Threshold': [s >= threshold for s in norm_scores]
            })

            # Display detailed results
            st.subheader("Detailed Results")
            st.dataframe(results_df)

            # Visualization
            st.subheader("Score Distribution")
            fig = px.scatter(results_df, x='Section', y='Normalized Score',
                             color='Meets Threshold',
                             title="Normalized Scores by Section",
                             color_discrete_map={True: 'green', False: 'red'})
            fig.add_hline(y=threshold, line_dash="dash", line_color="blue",
                          annotation_text=f"Threshold: {threshold}")
            st.plotly_chart(fig)

        except Exception as e:
            st.error(f"Error running diagonal check: {str(e)}")


def parse_sections(file):
    """Helper function to parse uploaded section files"""
    # Implementation depends on your file format
    # This is a placeholder
    content = file.read().decode('utf-8')
    sections = content.strip().split('\n')
    return [(section,) for section in sections if section.strip()]


# CACHED VERSION FOR PERFORMANCE
@st.cache_data
def cached_diagonal_check(evolved_sections_tuple, reference_sections_tuple, threshold, matrix_name):
    """
    Cached version of diagonal check for better Streamlit performance.
    Note: Uses tuples instead of lists for hashability.
    """
    matrix_map = {
        "blosum62": parasail.blosum62,
        "nuc44": parasail.nuc44,
        "pam250": parasail.pam250
    }

    # Convert tuples back to lists
    evolved_sections = [list(section) for section in evolved_sections_tuple]
    reference_sections = [list(section) for section in reference_sections_tuple]

    return diagonal_check_pure(
        evolved_sections,
        reference_sections,
        threshold,
        matrix_map[matrix_name]
    )


import parasail
from tqdm import tqdm


def compare_sections_pure(evolved_sections, reference_sections, threshold,
                          full_evolved, full_reference, substitution_matrix,
                          get_alignment_func, extract_extended_context_func,
                          progress_callback=None):
    """
    Pure function version of compare_sections for Streamlit compatibility.

    Compare each evolved section (tuple: (section, offset)) to every reference section.
    For each pair, compute the raw alignment score using sw_stats and normalize it by dividing
    by the length of the evolved section. If normalized score ≥ threshold, retrieve detailed alignment info.

    Args:
        evolved_sections (list): List of (section, offset) tuples for evolved sequences
        reference_sections (list): List of (section, offset) tuples for reference sequences
        threshold (float): Threshold for normalized scores
        full_evolved (str): Complete evolved sequence
        full_reference (str): Complete reference sequence
        substitution_matrix: Parasail substitution matrix object
        get_alignment_func (callable): Function to get detailed alignment
        extract_extended_context_func (callable): Function to extract extended context
        progress_callback (callable, optional): Function to call for progress updates

    Returns:
        list: List of match dictionaries with detailed alignment information
    """
    results = []
    match_counter = 0
    total_comparisons = len(evolved_sections) * len(reference_sections)

    comparison_count = 0
    for i, (evo_sec, evo_offset) in enumerate(evolved_sections, start=1):
        for j, (ref_sec, ref_offset) in enumerate(reference_sections, start=1):
            # Update progress if callback provided
            if progress_callback:
                progress = comparison_count / total_comparisons
                message = f"Comparing section {i} vs {j} ({comparison_count + 1}/{total_comparisons})"
                progress_callback(progress, message)

            stats_result = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix)
            raw_score = stats_result.score
            normalized_score = raw_score / len(evo_sec)

            if normalized_score >= threshold:
                (aligned_evo, aligned_ref, midline, align_score,
                 begin_query, end_query, begin_ref, end_ref) = get_alignment_func(evo_sec, ref_sec)

                extended_evo = extract_extended_context_func(
                    full_evolved, evo_offset, begin_query, end_query, extension=200
                )
                extended_ref = extract_extended_context_func(
                    full_reference, ref_offset, begin_ref, end_ref, extension=200
                )

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

            comparison_count += 1

    # Final progress update
    if progress_callback:
        progress_callback(1.0, f"Complete: Found {len(results)} matches")

    return results


# ADAPTOR FUNCTION - Use this if you need to maintain the original interface
def compare_sections(evolved_sections, reference_sections, threshold,
                     full_evolved, full_reference, substitution_matrix=None,
                     get_alignment_func=None, extract_extended_context_func=None):
    """
    Adaptor function to maintain original interface with tqdm progress bar.

    Args:
        evolved_sections (list): List of (section, offset) tuples for evolved sequences
        reference_sections (list): List of (section, offset) tuples for reference sequences
        threshold (float): Threshold for normalized scores
        full_evolved (str): Complete evolved sequence
        full_reference (str): Complete reference sequence
        substitution_matrix: Optional substitution matrix (will create default if None)
        get_alignment_func (callable, optional): Function to get detailed alignment
        extract_extended_context_func (callable, optional): Function to extract extended context

    Returns:
        list: Same as compare_sections_pure
    """
    if substitution_matrix is None:
        substitution_matrix = parasail.blosum62

    # Default functions if not provided (assuming they exist in scope)
    if get_alignment_func is None:
        get_alignment_func = get_alignment  # Assumes original get_alignment exists
    if extract_extended_context_func is None:
        extract_extended_context_func = extract_extended_context  # Assumes original exists

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
                 begin_query, end_query, begin_ref, end_ref) = get_alignment_func(evo_sec, ref_sec)

                extended_evo = extract_extended_context_func(
                    full_evolved, evo_offset, begin_query, end_query, extension=200
                )
                extended_ref = extract_extended_context_func(
                    full_reference, ref_offset, begin_ref, end_ref, extension=200
                )

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


# STREAMLIT USAGE EXAMPLE
def streamlit_compare_sections_interface():
    """
    Example of how to use this in a Streamlit app.
    This would be called from your main Streamlit script.
    """
    import streamlit as st
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go

    st.title("Section Comparison Analysis")

    # Input sections
    st.subheader("Input Data")
    col1, col2 = st.columns(2)

    with col1:
        st.write("**Evolved Sections**")
        num_evo = st.number_input("Number of evolved sections", min_value=1, value=3, key="num_evo")
        evolved_sections = []
        full_evolved = ""

        for i in range(num_evo):
            section = st.text_area(f"Evolved Section {i + 1}", key=f"evo_sec_{i}")
            offset = st.number_input(f"Offset {i + 1}", min_value=0, value=i * 100, key=f"evo_off_{i}")
            if section:
                evolved_sections.append((section, offset))
                full_evolved += section

    with col2:
        st.write("**Reference Sections**")
        num_ref = st.number_input("Number of reference sections", min_value=1, value=3, key="num_ref")
        reference_sections = []
        full_reference = ""

        for i in range(num_ref):
            section = st.text_area(f"Reference Section {i + 1}", key=f"ref_sec_{i}")
            offset = st.number_input(f"Offset {i + 1}", min_value=0, value=i * 100, key=f"ref_off_{i}")
            if section:
                reference_sections.append((section, offset))
                full_reference += section

    # Parameters
    st.subheader("Parameters")
    col1, col2 = st.columns(2)

    with col1:
        threshold = st.slider("Threshold", min_value=0.0, max_value=100.0, value=10.0, step=0.1)

    with col2:
        matrix_option = st.selectbox("Substitution Matrix", ["blosum62", "nuc44", "pam250"])

    matrix_map = {
        "blosum62": parasail.blosum62,
        "nuc44": parasail.nuc44,
        "pam250": parasail.pam250
    }

    # Mock functions for the example (replace with your actual implementations)
    def mock_get_alignment(evo_sec, ref_sec):
        """Mock alignment function - replace with your actual get_alignment"""
        return (evo_sec, ref_sec, "|" * min(len(evo_sec), len(ref_sec)),
                100, 0, len(evo_sec), 0, len(ref_sec))

    def mock_extract_extended_context(full_seq, offset, begin, end, extension=200):
        """Mock context extraction - replace with your actual function"""
        start = max(0, offset + begin - extension)
        stop = min(len(full_seq), offset + end + extension)
        return full_seq[start:stop]

    if st.button("Run Comparison") and evolved_sections and reference_sections:
        try:
            # Create progress indicators
            progress_bar = st.progress(0)
            status_text = st.empty()

            def update_progress(progress, message):
                progress_bar.progress(progress)
                status_text.text(message)

            # Run the comparison
            matches = compare_sections_pure(
                evolved_sections,
                reference_sections,
                threshold,
                full_evolved,
                full_reference,
                matrix_map[matrix_option],
                mock_get_alignment,
                mock_extract_extended_context,
                progress_callback=update_progress
            )

            # Clear progress indicators
            progress_bar.empty()
            status_text.empty()

            # Display results
            st.subheader("Results Summary")
            total_comparisons = len(evolved_sections) * len(reference_sections)

            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Matches", len(matches))
            with col2:
                st.metric("Total Comparisons", total_comparisons)
            with col3:
                match_rate = (len(matches) / total_comparisons) * 100 if total_comparisons > 0 else 0
                st.metric("Match Rate", f"{match_rate:.1f}%")

            if matches:
                # Create results DataFrame
                results_df = pd.DataFrame([
                    {
                        'Match Index': match['match_index'],
                        'Evolved Index': match['evo_index'],
                        'Reference Index': match['ref_index'],
                        'Raw Score': match['raw_alignment_score'],
                        'Normalized Score': match['normalized_score'],
                        'Alignment Score': match['alignment_score']
                    }
                    for match in matches
                ])

                # Display matches table
                st.subheader("Match Details")
                st.dataframe(results_df)

                # Visualization
                st.subheader("Match Visualization")

                # Heatmap of matches
                fig = go.Figure(data=go.Heatmap(
                    z=[[0 for _ in reference_sections] for _ in evolved_sections],
                    x=[f"Ref {i + 1}" for i in range(len(reference_sections))],
                    y=[f"Evo {i + 1}" for i in range(len(evolved_sections))],
                    colorscale='Blues'
                ))

                # Add match points
                for match in matches:
                    fig.add_scatter(
                        x=[match['ref_index'] - 1],
                        y=[match['evo_index'] - 1],
                        mode='markers',
                        marker=dict(size=match['normalized_score'] * 2, color='red'),
                        name=f"Score: {match['normalized_score']:.2f}",
                        showlegend=False
                    )

                fig.update_layout(
                    title="Section Matches Heatmap",
                    xaxis_title="Reference Sections",
                    yaxis_title="Evolved Sections"
                )
                st.plotly_chart(fig)

                # Score distribution
                fig2 = px.histogram(results_df, x='Normalized Score',
                                    title="Distribution of Normalized Scores")
                st.plotly_chart(fig2)

                # Detailed match viewer
                st.subheader("Detailed Match Viewer")
                if matches:
                    selected_match = st.selectbox(
                        "Select match to view details",
                        range(len(matches)),
                        format_func=lambda
                            x: f"Match {x}: Evo{matches[x]['evo_index']} vs Ref{matches[x]['ref_index']} (Score: {matches[x]['normalized_score']:.2f})"
                    )

                    match = matches[selected_match]
                    st.write(f"**Alignment Score:** {match['alignment_score']}")
                    st.write(f"**Normalized Score:** {match['normalized_score']:.3f}")

                    st.text("Alignment:")
                    alignment_text = f"Query:     {match['aligned_evo']}\nMidline:   {match['midline']}\nReference: {match['aligned_ref']}"
                    st.code(alignment_text)

            else:
                st.info("No matches found above the specified threshold.")

        except Exception as e:
            st.error(f"Error running comparison: {str(e)}")


# CACHED VERSION FOR PERFORMANCE
@st.cache_data
def cached_compare_sections(evolved_sections_tuple, reference_sections_tuple,
                            threshold, full_evolved, full_reference, matrix_name):
    """
    Cached version for better Streamlit performance.
    Note: Uses tuples for hashability and mock functions for demo.
    """
    matrix_map = {
        "blosum62": parasail.blosum62,
        "nuc44": parasail.nuc44,
        "pam250": parasail.pam250
    }

    def mock_get_alignment(evo_sec, ref_sec):
        return (evo_sec, ref_sec, "|" * min(len(evo_sec), len(ref_sec)),
                100, 0, len(evo_sec), 0, len(ref_sec))

    def mock_extract_extended_context(full_seq, offset, begin, end, extension=200):
        start = max(0, offset + begin - extension)
        stop = min(len(full_seq), offset + end + extension)
        return full_seq[start:stop]

    # Convert tuples back to lists
    evolved_sections = list(evolved_sections_tuple)
    reference_sections = list(reference_sections_tuple)

    return compare_sections_pure(
        evolved_sections,
        reference_sections,
        threshold,
        full_evolved,
        full_reference,
        matrix_map[matrix_name],
        mock_get_alignment,
        mock_extract_extended_context
    )

def plot_interactive_matrix(num_chunks, matches, match_map, threshold):
    """
    Build and display an interactive square grid (num_chunks x num_chunks) for above-threshold matches.
    The color scale is relative (from threshold to 1) so differences above threshold are visible.
    A loading bar is shown while building the matrix.
    Clicking on a cell (with a match) opens a Tkinter window with a scrollable text widget showing:
      - The aligned reference and evolved sequences (with connecting midline if applicable)
      - The code to access the corresponding chunk matrix (e.g., "chunk_matrices[<match_index>]")
    """
    # Build the score matrix only for above-threshold matches.
    # Build the score matrix only for above-threshold matches, with a progress bar
    score_matrix = np.zeros((num_chunks, num_chunks))
    for match in tqdm.tqdm(matches,
                      total=len(matches),
                      desc="Building interactive matrix",
                      unit="match",
                      leave=False):
        i = match['evo_index'] - 1
        j = match['ref_index'] - 1
        score_matrix[i, j] = match['normalized_score']

    # Create a colormap scaled from threshold to 1.
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(score_matrix, cmap='viridis', interpolation='nearest', origin='upper', vmin=threshold, vmax=1)
    ax.set_title("Normalized Alignment Score Matrix (Matches Only)")
    ax.set_xlabel("Reference Chunk Index")
    ax.set_ylabel("Evolved Chunk Index")
    ax.set_xticks(np.arange(num_chunks))
    ax.set_yticks(np.arange(num_chunks))
    ax.set_xticklabels([str(x + 1) for x in range(num_chunks)])
    ax.set_yticklabels([str(x + 1) for x in range(num_chunks)])
    plt.colorbar(im, ax=ax, label="Normalized Score")

    def on_click(event):
        if event.xdata is None or event.ydata is None:
            return
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        if row < 0 or row >= num_chunks or col < 0 or col >= num_chunks:
            return
        key = (row + 1, col + 1)  # 1-indexed keys
        if key in match_map:
            match = match_map[key]

            # Create a new Toplevel window
            win = tk.Toplevel(root)
            win.title(f"Chunk Match: Evolved {key[0]} vs. Reference {key[1]}")

            # Create a horizontally scrollable text widget
            st = scrolledtext.ScrolledText(
                win,
                wrap="none",  # no wrapping => horizontal scrolling
                width=100,
                height=10,
                font=("Courier", 12)
            )
            st.pack(expand=True, fill="both")

            # Build exactly three lines: reference, midline, evolved
            lines = [
                match["aligned_ref"],
                match["midline"],
                match["aligned_evo"],
                "",  # blank line
                f"Access code: chunk_matrices[{match['match_index']}]"
            ]
            st.insert("1.0", "\n".join(lines))
        #else:
            #print(f"No match found for cell (Evolved {row + 1}, Reference {col + 1})")
            #return

    fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()
    return score_matrix

def print_alignment_summary(match, snippet_length=100):
    """Print a summary of a matching pair."""
    print(f"\nEvolved Section {match['evo_index']} vs. Reference Section {match['ref_index']}:")
    print(f"  Raw Alignment Score: {match['raw_alignment_score']}")
    print(f"  Normalized Alignment Score: {match['normalized_score']:.3f}")
    print("  Alignment (snippet):")
    if len(match["aligned_evo"]) > snippet_length:
        print(match["aligned_evo"][:snippet_length] + "...")
        print(match["midline"][:snippet_length] + "...")
        print(match["aligned_ref"][:snippet_length] + "...")
    else:
        print(match["aligned_evo"])
        print(match["midline"])
        print(match["aligned_ref"])
    print("  Extended Context (Evolved, snippet):", match["extended_evo"][:snippet_length] + "...")
    print("  Extended Context (Reference, snippet):", match["extended_ref"][:snippet_length] + "...")

#--------UI HELPERS------------------------------------------------------------------------------------

AMBIGUOUS_DNA = set("ACGTURYSWKMBDHVN-")  # includes RNA U and IUPAC symbols + gaps

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