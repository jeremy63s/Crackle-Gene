import matplotlib

matplotlib.use('TkAgg')  # Force TkAgg backend
from Bio.Blast import NCBIWWW, NCBIXML

import parasail
import string
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm  # Standard tqdm for PyCharm
import os
import tkinter as tk
from tkinter import scrolledtext

import os
import subprocess
import sys
import platform
import orfipy
import orfipy.findorfs as of
import orfipy_core
from Bio.Seq import Seq
print(os.getcwd())
# Create one root window and hide it.
root = tk.Tk()
root.withdraw()

# Create a ion matrix for uppercase letters.

substitution_matrix = parasail.matrix_create(ALPHABET, 1, -1)




def read_fasta(filepath):
    """Read a FASTA file and return the sequence as a single string."""
    try:
        with open(filepath, 'r') as file:
            seq_lines = []
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    continue  # Skip header lines
                seq_lines.append(line)
        return "".join(seq_lines)
    except Exception as e:
        raise ValueError("Error reading FASTA file: " + str(e))


def get_sequence(seq_type):
    """Get a sequence either by manual input or by reading a FASTA file."""
    method = input(
        f"Do you want to input the {seq_type} sequence manually or via FASTA file? (M for manual, F for FASTA): ").strip().upper()
    if method == 'F':
        filepath = input(f"Enter the path to the {seq_type} FASTA file: ").strip()
        seq = read_fasta(filepath)
    else:
        seq = input(f"Enter the {seq_type} sequence: ").strip()
    return seq.upper()


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



from tqdm import tqdm   # ← make sure this is imported at the top alongside your other imports

def diagonal_check(evolved_sections, reference_sections, threshold):
    """
    Fast check: compare section i vs section i for all i,
    returns (pct_met, max_consec_failed, raw_scores, norm_scores).
    """
    n = len(evolved_sections)
    raw_scores = []
    norm_scores = []

    for (evo_sec, _), (ref_sec, _) in tqdm(
            zip(evolved_sections, reference_sections),
            total=n,
            desc="Diagonal check",
            leave=False
        ):
        # compute once and store both raw and normalized
        raw = parasail.sw_stats(evo_sec, ref_sec, 1, 1, substitution_matrix).score
        raw_scores.append(raw)
        norm_scores.append(raw / len(evo_sec))

    # compute percent met
    meets = [s >= threshold for s in norm_scores]
    pct_met = sum(meets) / n * 100

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
    return indicator, indicator


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
    for match in tqdm(matches,
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
        else:
            print(f"No match found for cell (Evolved {row + 1}, Reference {col + 1})")

    fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()
    return score_matrix


# define final matrix:
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
    print("concatenated shape is", concatenated.shape)
    print("counter_row shape is", counter_row.shape)
    return final_matrix


def on_pick(event):
    # Ensure that the picked artist is our scatter plot.
    if event.artist != scatter:
        return

    # event.ind is a list of indices in the plotted (filtered) array.
    picked_index = event.ind[0]
    # Convert to the original data index.
    click_idx = plotted_indices[picked_index]

    # Define the slice: 10 columns to the left and right of the clicked column.
    start_idx = max(0, click_idx - 10)
    end_idx = min(data.shape[1], click_idx + 11)  # end index is exclusive
    subset = data[:, start_idx:end_idx]

    # Determine which column in the subset was clicked.
    center_cell_index = click_idx - start_idx

    # Build a fixed-width text display.
    row_label_width = 10  # Fixed width for row labels.
    cell_width = 12       # Fixed width for each cell.
    row_lines = []
    for i, row in enumerate(subset):
        # Create a fixed-width row label.
        row_label = f"Row {i}:".ljust(row_label_width)
        # Format each cell to a fixed width.
        cells = [f"{str(cell):{cell_width}}" for cell in row]
        # Join the cells with a space.
        line = row_label + " " + " ".join(cells)
        row_lines.append(line)
    subset_str = "\n".join(row_lines) + "\n"

    # Use the existing Tk root if available; otherwise create one.
    root = tk._default_root
    if root is None:
        root = tk.Tk()

    # Create a new window for the text display.
    win = tk.Toplevel(root)
    win.title(f"Data Matrix Slice (Columns {start_idx} to {end_idx-1})")

    # Create a ScrolledText widget with no wrapping and a larger monospaced font.
    text_widget = scrolledtext.ScrolledText(win, wrap='none', font=("Courier", 14))
    text_widget.insert(tk.END, subset_str)
    text_widget.configure(state='normal')  # Enable editing temporarily for tagging.

    # Highlight rows 2, 5, and 6 (if they exist in the subset) using yellow.
    highlight_rows = [2, 5, 6]
    for row_idx in highlight_rows:
        if row_idx < subset.shape[0]:
            line_no = row_idx + 1  # Text widget line numbers start at 1.
            tag_name = f"highlight_row_{row_idx}"
            text_widget.tag_add(tag_name, f"{line_no}.0", f"{line_no}.end")
            text_widget.tag_config(tag_name, background="yellow")

    # Now, highlight the center column (the clicked column) in light blue.
    # The row label occupies row_label_width characters plus one space.
    start_offset = row_label_width + 1
    for i in range(subset.shape[0]):
        line_no = i + 1
        # Calculate the start and end character positions for the center cell.
        col_start = start_offset + center_cell_index * (cell_width + 1)
        col_end = col_start + cell_width
        text_widget.tag_add("center_col", f"{line_no}.{col_start}", f"{line_no}.{col_end}")
    text_widget.tag_config("center_col", background="lightblue")

    text_widget.configure(state='disabled')  # Make the widget read-only.
    text_widget.pack(expand=True, fill='both')


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

# --- Helper Functions ---
def array_to_str(arr):
    """Convert an array of single characters to a string."""
    return "".join(arr)

def create_annotation_array(seq_len):
    """Return a list of length seq_len filled with '0'."""
    return ['0'] * seq_len

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



    # SAVE INFO MATRIX
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

    # NCBI BLAST:
def extract_protein_seq(info_matrix, a, b, c):
        """
        Given info_matrix and slice indices a (row), b (start col), c (end col),
        return the protein_seq from that ORF.
        """
        # slice out the triplet letters
        AA_triplocate = info_matrix[a, b:c]
        # take every 3rd position, starting at 0
        AA_seq = AA_triplocate[0::3]
        # join into a string
        return ''.join(AA_seq)

def process_orf_array(info_matrix, orf_array, primary_row, secondary_row, orf_code):
        """
        Given a 1D array (orf_array) of mutation location indices for an ORF, this function
        processes every value in that array. For each value n, it uses find_boundaries_for_row
        on both the primary and secondary rows (from info_matrix). If the two results differ,
        a warning is issued and the primary row’s boundaries are used.

        Duplicate boundary pairs (within the same ORF) are ignored (i.e. once a block is recorded,
        subsequent mutation indices within the same block print "mutation already accounted for").

        Returns:
          A list of rows, each row: [ ORF_code, left_boundary, right_boundary ]
        """
        unique_boundaries = set()
        output = []

        for n in orf_array:
            try:
                boundaries_primary = find_boundaries_for_row(info_matrix, n, primary_row)
                boundaries_secondary = find_boundaries_for_row(info_matrix, n, secondary_row)
            except ValueError as err:
                # Skip this index if a boundary cannot be computed.
                print(f"Skipping index {n} for ORF code {orf_code}: {err}")
                continue

            if boundaries_primary != boundaries_secondary:
                print(
                    f"Warning for ORF code {orf_code} at mutation index {n}: boundaries differ between rows {primary_row} and {secondary_row}. Using row {primary_row} boundaries.")
            boundaries = boundaries_primary

            if boundaries in unique_boundaries:
                print(f"For ORF code {orf_code} at mutation index {n}: mutation already accounted for")
            else:
                unique_boundaries.add(boundaries)
                output.append([orf_code, boundaries[0], boundaries[1]])
                print(
                    f"For ORF code {orf_code} at mutation index {n}: boundaries: Left = {boundaries[0]}, Right = {boundaries[1]}")

        return output

def blast_protein_seq(seq, email=None):
        """
        Run NCBIWWW.qblast on a single protein sequence (BLASTP against nr),
        returning the parsed record.
        """
        fasta = f">query\n{seq}"
        handle = NCBIWWW.qblast("blastp", "nr", fasta, format_type="XML")
        return NCBIXML.read(handle)

def process_result_matrix(info_matrix, result_matrix, do_blast=True):
        """
        For each row [orf_code, left, right] in result_matrix:
          - map orf_code → a row index
          - extract protein_seq
          - optionally BLAST and print top hits
        """
        a_map = {
            1: 14,
            2: 15,
            3: 16,
            4: 20,
            5: 21,
            6: 22
        }
        for idx, (orf_code, left, right) in enumerate(result_matrix):
            try:
                a_row = a_map[int(orf_code)]
            except KeyError:
                print(f"Row {idx}: unknown ORF code {orf_code}, skipping.")
                continue

            # extract
            protein_seq = extract_protein_seq(info_matrix, a_row, int(left), int(right))
            print(f"\nRow {idx} (ORF {orf_code}): cols {left}–{right}")
            print("Protein sequence:", protein_seq)

            if do_blast and protein_seq:
                print("Submitting to BLAST…")
                record = blast_protein_seq(protein_seq)
                print("Top 5 BLAST hits for this query:")
                for aln in record.alignments[:5]:
                    hsp = aln.hsps[0]
                    print(f" {aln.accession}\t{aln.hit_def}\tE={hsp.expect}")
            elif do_blast:
                print("Empty sequence—skipping BLAST.")
    # Save the matrix to an absolute path accessible to your system.


def main():
    print("Welcome to the Genome Sequence Comparator!")
    # Get sequences.
    ref_seq = get_sequence("reference")
    evolved_seq = get_sequence("evolved")
    start_subseq = get_sequence("starting subsequence")

    # Rotate sequences.
    try:
        rotated_ref = rotate_sequence(ref_seq, start_subseq)
        rotated_evolved = rotate_sequence(evolved_seq, start_subseq)
    except ValueError as e:
        print("\nError during rotation:", e)
        return None
    print("\nRotated Reference Sequence:")
    print(rotated_ref)
    print("\nRotated Evolved Sequence:")
    print(rotated_evolved)

    # Divide sequences.
    try:
        n = int(input("\nEnter the number of sections to divide the sequences into: ").strip())
        if n <= 0:
            raise ValueError("Number of sections must be positive.")
    except ValueError as e:
        print("Invalid input for number of sections:", e)
        return None
    ref_sections = divide_sequence_with_offset(rotated_ref, n)
    evolved_sections = divide_sequence_with_offset(rotated_evolved, n)
    print("\nReference Sequence Sections (with offsets):")
    for idx, (sec, offset) in enumerate(ref_sections, start=1):
        print(f"Section {idx} (offset {offset}): {sec}")
    print("\nEvolved Sequence Sections (with offsets):")
    for idx, (sec, offset) in enumerate(evolved_sections, start=1):
        print(f"Section {idx} (offset {offset}): {sec}")

    # Ask for a normalized alignment score threshold.
#    try:
#        threshold = float(input("\nEnter the normalized alignment score threshold (0 to 1, e.g., 0.9): ").strip())
#    except ValueError as e:
#        print("Invalid threshold:", e)
#        return None

        # Ask for a normalized alignment score threshold.
    try:
        threshold = float(input("\nEnter the normalized alignment score threshold (0 to 1, e.g., 0.9): ").strip())
    except ValueError as e:
        print("Invalid threshold:", e)
        return None

    # ----- NEW: diagonal pre-check -----
    # … after reading threshold from input …

    # NEW unpacking
    pct_met, max_consec_fail, raw_scores, diag_scores = diagonal_check(
        evolved_sections, ref_sections, threshold
    )

    print(f"\nDiagonal check:")
    print(f"  {pct_met:.1f}% of chunks on the diagonal met the alignment threshold")
    print(f"  Longest run of consecutive failures on diagonal: {max_consec_fail} chunks")
    if max_consec_fail > 2 and pct_met < 90:
        print("  Recommend running full analysis to check for translocations")

    run_full = input("\nRun full n×n comparison? (Y/N): ").strip().upper()
    do_full = (run_full == 'Y')

    if do_full:
        matches = compare_sections(
            evolved_sections, ref_sections,
            threshold, rotated_evolved, rotated_ref
        )
    else:
        matches = []
        for i, norm_score in enumerate(tqdm(diag_scores, desc="Checking chunks"), start=1):
            if norm_score < threshold:
                continue

            # reuse the precomputed raw_score
            raw = raw_scores[i - 1]

            evo_sec, evo_off = evolved_sections[i - 1]
            ref_sec, ref_off = ref_sections[i - 1]

            # detailed alignment (only done when needed)
            (aligned_evo, aligned_ref, midline,
             align_score, bq, eq, br, er) = get_alignment(evo_sec, ref_sec)

            # extended contexts
            ext_evo = extract_extended_context(rotated_evolved, evo_off, bq, eq)
            ext_ref = extract_extended_context(rotated_ref, ref_off, br, er)

            match = {
                "evo_index": i,
                "ref_index": i,
                "raw_alignment_score": raw,
                "normalized_score": norm_score,
                "alignment_score": align_score,
                "aligned_evo": aligned_evo,
                "aligned_ref": aligned_ref,
                "midline": midline,
                "begin_query": bq,
                "end_query": eq,
                "begin_ref": br,
                "end_ref": er,
                "evo_offset": evo_off,
                "ref_offset": ref_off,
                "extended_evo": ext_evo,
                "extended_ref": ext_ref,
                "is_match": 1,
                "match_index": len(matches),
            }
            matches.append(match)
    print("aligned_evo length is: ", len(aligned_evo), "aligned ref length is: ", len(aligned_ref))
    # 4) Build chunk_matrices and match_map as before
    chunk_matrices = [
        build_chunk_matrix(m)
        for m in tqdm(matches, desc="Building matrices")
    ]
    match_map = {(m["evo_index"], m["ref_index"]): m for m in matches}

    # 5) Interactive matrix (this will just show the diagonal if do_full is False)
    num_chunks = len(evolved_sections)
    print("\nBuilding interactive matrix visualization...")
    interactive_matrix = plot_interactive_matrix(
        num_chunks, matches, match_map, threshold
    )

    # 6) Final matrix + saving
    final_matrix = build_final_matrix(chunk_matrices)
    # Save Chunk Matrices to .npy file.
    # chunk_file_path = "/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/chunk_matrices.npy"
    # np.save(chunk_file_path, chunk_matrices)
    # if os.path.exists(chunk_file_path):
    #    print("Chunk matrices saved successfully at:", chunk_file_path)
    # else:
    #    print("Error: Chunk matrices file was not saved.")

    # Save Final Matrix to .npy file.
    final_matrix_path = "/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/final_matrix.npy"
    np.save(final_matrix_path, final_matrix)
    if os.path.exists(final_matrix_path):
        print("Final matrix saved successfully at:", final_matrix_path)
    else:
        print("Error: Final matrix file was not saved.")

    # Also, save the final matrix as a text file (optional).
    # final_matrix_txt_path = "/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/final_matrix.txt"
    # np.savetxt(final_matrix_txt_path, final_matrix, fmt='%s', delimiter='\t')
    # if os.path.exists(final_matrix_txt_path):
    #    print("Final matrix text file saved successfully at:", final_matrix_txt_path)
    # else:
    #    print("Error: Final matrix text file was not saved.")
    #final_matrix = results["final_matrix"]
    print(type(final_matrix), final_matrix.shape)
    final_matrix = np.array(final_matrix)
    print("final_matrix is a", type(final_matrix))
    print("final_matrix is shape", final_matrix.shape)

    # 1) Diagnose

    # if isinstance(final_matrix, dict):
    #   print(" keys:", final_matrix.keys())

    # 2) Convert to 2D array
    # if isinstance(final_matrix, dict):
    #    rows = [final_matrix[k] for k in sorted(final_matrix.keys())]
    #    data = np.array(rows)
    # else:
    #    data = np.array(final_matrix)

    # 3) Confirm
    # print("data.shape =", getattr(data, "shape", None))
    data = np.array(final_matrix)
    print("data is shape", data.shape)
    print("data is type", type(data))
    # 4) Find indices in 7th row (index 6)

    indices = np.where(data[6] == -1)[0]
    # print("–1’s at columns:", indices)
    # where are the indels and point mutations
    print("Indices where the value is -1:")
    for idx in indices:
        print(f"Index: {idx}, Value: {data[5, idx]}")

    # MAKING GRAPH

    # --- Assume your data matrix 'data' is already loaded ---
    # For example:
    # data = np.load('final_matrix.npy', allow_pickle=True)

    # Create a mask to filter points (only plot where data[6, :] != 1).
    mask = data[6, :] != 1
    # Get the original indices for the points we are plotting.
    plotted_indices = np.where(mask)[0]
    # Use these to obtain the x and y values.
    x_values = data[7, :][mask]  # x-axis from row 7
    y_values = data[6, :][mask]  # y-axis from row 6

    # Create the main scatter plot and enable picking.
    fig, scatter_ax = plt.subplots(figsize=(10, 6))
    scatter = scatter_ax.scatter(x_values, y_values, color='blue', picker=True)
    scatter_ax.set_xlabel('data[7, :]')
    scatter_ax.set_ylabel('Mutation Type')
    scatter_ax.set_title('Mutation Type vs. data[7, :] (Filtered)')
    scatter_ax.set_yticks([-1, 0])
    scatter_ax.set_yticklabels(['substitutions', 'indels'])

    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()

    # Extract the nucleotide sequences (as lists) for translation:
    ref_seq = list(final_matrix[1, :])  # Reference nucleotide sequence
    evolved_seq = list(final_matrix[3, :])  # Evolved nucleotide sequence

    # Create the six new rows using active coding region translation:
    aa_ref = translate_nuc_row_active(ref_seq, CODON_TABLE)  # Amino acids for reference sequence (active regions only)
    aa_evolved = translate_nuc_row_active(evolved_seq,
                                          CODON_TABLE)  # Amino acids for evolved sequence (active regions only)
    ref_revcomp = reverse_complement(ref_seq)  # Reverse complement of reference sequence
    aa_ref_revcomp = translate_nuc_row_active(ref_revcomp,
                                              CODON_TABLE)  # Amino acids for reverse complement of reference
    evolved_revcomp = reverse_complement(evolved_seq)  # Reverse complement of evolved sequence
    aa_evolved_revcomp = translate_nuc_row_active(evolved_revcomp,
                                                  CODON_TABLE)  # Amino acids for reverse complement of evolved

    # Build the new rows as a list of lists
    new_rows = [
        aa_ref,  # Row for amino acids from the reference sequence
        aa_evolved,  # Row for amino acids from the evolved sequence
        ref_revcomp,  # Row for reverse complement of the reference sequence
        aa_ref_revcomp,  # Row for amino acids from the reverse complement of the reference
        evolved_revcomp,  # Row for reverse complement of the evolved sequence
        aa_evolved_revcomp  # Row for amino acids from the reverse complement of the evolved sequence
    ]

    # Append the new rows to the existing final_matrix using vertical stacking
    final_matrix = np.vstack([final_matrix] + new_rows)

    # final_matrix now contains the original 8 rows followed by 6 new rows:
    #   Row 8: Amino acid translation for the reference (active coding regions only)
    #   Row 9: Amino acid translation for the evolved (active coding regions only)
    #   Row 10: Reverse complement of the reference sequence
    #   Row 11: Amino acid translation for the reverse complement of the reference (active coding regions only)
    #   Row 12: Reverse complement of the evolved sequence
    #   Row 13: Amino acid translation for the reverse complement of the evolved (active coding regions only)

    print(final_matrix[:, 3])
    # Now final_matrix has six additional rows at the bottom.
    # Initialize an empty list for storing the desired elements.


    ForwardAAMuts = []
    # ForwardAAMuts now contains every final_matrix[4, n] where final_matrix[8, n] and final_matrix[9, n] differ.
    # Iterate over each column index in the final_matrix.
    for n in range(final_matrix.shape[1]):
        # Check if the amino acid in row 8 (reference) is not equal to the one in row 9 (evolved)
        if final_matrix[8, n] != final_matrix[9, n]:
            ForwardAAMuts.append(final_matrix[7, n])

    # Optionally, convert ForwardAAMuts to a NumPy array
    ForwardAAMuts = np.array(ForwardAAMuts)

    print(ForwardAAMuts)

    aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_to_num = {aa: i + 1 for i, aa in enumerate(aa_order)}

    # Prepare lists for the x and y coordinates of points.
    x_coords = []
    y_coords = []

    # Loop over each column of final_matrix.
    for n in range(final_matrix.shape[1]):
        aa_val = final_matrix[8, n]  # Amino acid from row 8
        if aa_val != '0':
            # Get the relative position from row 7 (convert to float if needed)
            x_val = float(final_matrix[7, n])
            # Map the amino acid to its corresponding numeric value
            y_val = aa_to_num.get(aa_val, None)
            if y_val is not None:
                x_coords.append(x_val)
                y_coords.append(y_val)

    # Create the dot plot.
    # plt.figure(figsize=(12, 6))
    # plt.scatter(x_coords, y_coords, marker='o')
    # plt.xlabel('Relative Position (final_matrix[7, :])')
    # plt.ylabel('Amino Acid (alphabetical order)')
    # plt.title('Dot Plot of Amino Acid Translations (Reference Sequence)')

    # Set y-axis ticks to show amino acid letters.
    # plt.yticks(list(aa_to_num.values()), list(aa_to_num.keys()))
    # plt.grid(True)
    # plt.show()

    # Save the final_matrix to a .npy file named final_matrix14.npy
    #np.save('final_matrix14.npy', final_matrix)

    # Later, you can load it back using:
    # loaded_matrix = np.load('final_matrix14.npy')

    ###########################
    ##### END CODE BLOCK 3####
    ##########################

    loaded_matrix = final_matrix
    # print(loaded_matrix[:, 66120])

    from Bio import SeqIO

    fasta_file = '/Users/benumlauf/Desktop/Reference_Sequence.fa'

    # Parse the FASTA file and print each record's ID and sequence.

    # REPLACE ALL - WITH N

    # Assuming 'loaded_matrix' is already defined as a NumPy array.
    # Create a copy of loaded_matrix to preserve the original data.
    N_matrix = loaded_matrix.copy()

    # Replace '-' with 'N' in row index 1
    N_matrix[1, :] = np.where(N_matrix[1, :] == '-', 'N', N_matrix[1, :])

    # Replace '-' with 'N' in row index 3
    N_matrix[3, :] = np.where(N_matrix[3, :] == '-', 'N', N_matrix[3, :])

    N_matrix[10, :] = np.where(N_matrix[10, :] == '-', 'N', N_matrix[10, :])

    N_matrix[12, :] = np.where(N_matrix[12, :] == '-', 'N', N_matrix[12, :])

    # Now, N_matrix has the modifications, and all other elements remain unchanged.
    # print(N_matrix[:, 66120])

    # Parse the FASTA file and print each record's ID and sequence.
    # for record in SeqIO.parse(fasta_file, "fasta"):
    #    print("ID:", record.id)
    #    print("Sequence:", record.seq)

    np.save("/Users/benumlauf/Downloads/final_matrix14.npy", loaded_matrix)

    #############################
    # TURN .GBK INTO ARRAY     #
    ############################
    from Bio import SeqIO

    # Use the GenBank file exported from SnapGene that includes ORF annotations.
    file_path = "/Users/benumlauf/Desktop/RefSeqGB.gb"
    record = SeqIO.read(file_path, "genbank")

    # Iterate over annotated ORFs (or CDS features)
    for feature in record.features:
        if feature.type in ["ORF", "CDS"]:
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            # Extract the ORF sequence (handles reverse strands as needed)
            orf_seq = feature.extract(record.seq)
            # Translate the ORF sequence
            aa_seq = orf_seq.translate(to_stop=True)

            print(f"Detected ORF ({feature.type}) at positions {start} to {end} (strand {strand})")
            print("First 10 amino acids:", aa_seq[:10])
            print("Full translation (up to stop):", aa_seq)
            print("-" * 60)
    seq_indices = {
        "forward_ref": 1,
        "forward_evolved": 3,
        "revcomp_ref": 10,
        "revcomp_evolved": 12
    }

    # Assume loaded_matrix is already defined (e.g., loaded via np.loadtxt or similar)
    # For each sequence, convert the row (an array of chars) into a string.
    sequences = {}
    for name, idx in seq_indices.items():
        sequences[name] = array_to_str(N_matrix[idx, :])

    # --- Process Each Sequence and Create Annotation Arrays ---
    # We'll store the output in a dictionary:
    # results[name][frame] will be the annotation array for that sequence and reading frame.
    results = {}

    # For each sequence, we produce three annotation arrays (one for each reading frame).
    for name, seq in sequences.items():
        seq_len = len(seq)
        # Initialize annotation arrays for frames 0, 1, 2 (which correspond to ORF frames 1,2,3).
        frame_annotations = {
            0: create_annotation_array(seq_len),
            1: create_annotation_array(seq_len),
            2: create_annotation_array(seq_len)
        }

        # Run ORFipy on the sequence. We search only the forward strand ('f') because
        # the provided sequences are already oriented as desired.
        for orf in orfipy_core.orfs(seq, minlen=3, maxlen=1000000, strand='f'):
            start, stop, strand_val, description = orf
            # Extract the ORF frame from the description (values 1, 2, or 3).
            frame_val = parse_orf_frame(description)
            if frame_val is None:
                continue  # Skip if we can't parse the frame

            # Convert ORF frame (1,2,3) to index offset (0,1,2)
            frame_idx = frame_val - 1

            # Extract the nucleotide sequence for this ORF.
            orf_seq = seq[start:stop]
            # Translate the ORF sequence (stops translation at the first stop codon).
            protein_seq = str(Seq(orf_seq).translate(to_stop=True))

            # For each codon in the ORF, fill in the corresponding three positions
            # in the annotation array with the amino acid letter.
            for i, aa in enumerate(protein_seq):
                codon_start = start + i * 3
                codon_end = codon_start + 3
                # Ensure we don't exceed the sequence length.
                for pos in range(codon_start, min(codon_end, seq_len)):
                    frame_annotations[frame_idx][pos] = aa

        results[name] = frame_annotations

    # --- Output the Annotation Arrays as Strings ---
    # There will be 12 arrays total (4 sequences × 3 reading frames).
    for seq_name, frames in results.items():
        print(f"Annotations for sequence: {seq_name}")
        for frame in range(3):
            ann_str = "".join(frames[frame])
            print(f"  Reading Frame {frame} (ORF frame {frame + 1}):")
            print(ann_str)
            print("-" * 60)

    # --- Create a 12x(seq_length) Matrix for Protein Annotations ---
    # Define the order in which sequences should be added.
    seq_order = ["forward_ref", "forward_evolved", "revcomp_ref", "revcomp_evolved"]

    protein_matrix_rows = []
    for seq_name in seq_order:
        frame_annotations = results[seq_name]
        # Append rows for each of the three reading frames in order.
        for frame in range(3):
            protein_matrix_rows.append(frame_annotations[frame])

    # Convert the list of rows into a NumPy array.
    protein_matrix = np.array(protein_matrix_rows)
    print("Protein Matrix Shape:", protein_matrix.shape)
    print('N matrix shape:', N_matrix.shape)
    # The protein_matrix now has 12 rows and each row has the same number of columns as the sequence length.

    # Append protein_matrix to the bottom of N_matrix
    # Ensure that both matrices have the same number of columns
    info_matrix = np.vstack((N_matrix, protein_matrix))
    print("information matrix shape:", info_matrix.shape)
    # --- right after info_matrix is ready in main() ---
    a_map = {
        1: 14,
        2: 15,
        3: 16,
        4: 20,
        5: 21,
        6: 22
    }
    orf_mapping = {
        "ORF1_forward": {"code": 1, "rows": (14, 17)},
        "ORF2_forward": {"code": 2, "rows": (15, 18)},
        "ORF3_forward": {"code": 3, "rows": (16, 19)},
        "ORF1_reverse": {"code": 4, "rows": (20, 23)},
        "ORF2_reverse": {"code": 5, "rows": (21, 24)},
        "ORF3_reverse": {"code": 6, "rows": (22, 25)}
    }

    all_results = []
    for orf_name, props in orf_mapping.items():
        code = props["code"]
        r1, r2 = props["rows"]

        # compute the *integer* indices array right here
        arr = np.where(info_matrix[r1, :] != info_matrix[r2, :])[0]

        # just to be sure it's int64
        arr = arr.astype(int)

        # now call your helper — n will be ints
        results = process_orf_array(info_matrix, arr, r1, r2, code)
        all_results.extend(results)

    info_matrix_cheatsheet = np.array([['0), Alignment data (chunkwise)'],
                                       ['1), Reference nucleotide sequence'],
                                       ['2), Reference nucleotide index'],
                                       ['3), Evolvded nucleotide sequence'],
                                       ['4), Evolved nucleotide index'],
                                       ['5), Frameshift counter (del: +1, ins: -1)'],
                                       ['6), Mut Type (in/del: 0, point: -1)'],
                                       ['7), Global index'],
                                       ['8), Amino acids in ref (naive ORF def)'],
                                       ['9), Amino acids in evo (naive ORF def)'],
                                       ['10), RevComp of Reference seq'],
                                       ['11), Amino acids in RevComp of ref (naive ORF def)'],
                                       ['12), RevComp of Evolved seq'],
                                       ['13), Amino acids in RevComp of evo (naive ORF def)'],
                                       ['14), Reference ORF frame 1'],
                                       ['15), Reference ORF frame 2'],
                                       ['16), Reference ORF frame 3'],
                                       ['17), Evolved ORF frame 1'],
                                       ['18), Evolved ORF frame 2'],
                                       ['19), Evolved ORF frame 3'],
                                       ['20), Reverse complement reference ORF frame 1'],
                                       ['21), Reverse complement reference ORF frame 2'],
                                       ['22), Reverse complement reference ORF frame 3'],
                                       ['23), Reverse complement evolved ORF frame 1'],
                                       ['24), Reverse complement evolved ORF frame 2'],
                                       ['25), Reverse complement evolved ORF frame 3']])


    global_path_info = '/Users/benumlauf/Desktop/JS/pyfolder/info_matrix.npy'
    np.save(global_path_info, info_matrix)

    # SAVE MATRIX CHEATSHEET

    global_path_cheatsheet = '/Users/benumlauf/Desktop/JS/pyfolder/info_matrix_cheatsheet.npy'
    np.save(global_path_cheatsheet, info_matrix_cheatsheet)

    # Create arrays for the forward orientation mutations (at single index).
    #ORF1_forward = np.where(info_matrix[14, :] != info_matrix[17, :])[0]
    #ORF2_forward = np.where(info_matrix[15, :] != info_matrix[18, :])[0]
    #ORF3_forward = np.where(info_matrix[16, :] != info_matrix[19, :])[0]

    # Create arrays for the reverse (reverse complement) orientation mutations (at single index).
    #ORF1_reverse = np.where(info_matrix[20, :] != info_matrix[23, :])[0]
    #ORF2_reverse = np.where(info_matrix[21, :] != info_matrix[24, :])[0]
    #ORF3_reverse = np.where(info_matrix[22, :] != info_matrix[25, :])[0]

    result_matrix = None
    # --- User: Ensure that the following six arrays are defined in your global namespace ---
    # Example definitions (remove or replace these with your actual data):
    # ORF1_forward = [84838, 84839, 84840]
    # ORF2_forward = [64489, 64490, 64491, 66118, 66119, 66120, 66121, 66122, 66123, 201877, 201878, 201879, 201880, 201881, 201882, 201883, 201884, 201885]
    # ORF3_forward = [84806, 84807, 84808, 132305, 132306, 132307, 201878, 201879, 201880, 201881, 201882, 201883]
    # ORF1_reverse = [88851, 88852, 88853]
    # ORF2_reverse = [136321, 136322, 136323, 136351, 136352, 136353, 136374, 155032, 155033, 155034, 155040]
    # ORF3_reverse = [19274, 19275, 19276, 19277, 19278, 19279, 88850, 88851, 88852, 136322, 136323, 136324, 155033, 155034, 155041]

    # --- Mapping of ORF array names to info_matrix rows and ORF codes ---

    #global result_matrix
    try:
        info_matrix
    except NameError:
        print("Error: info_matrix is not defined. Please pre-define it before running this script.")
        return

    #all_results = []
    #for orf_name, props in orf_mapping.items():
    #    arr = globals()[orf_name]
    #    code, (r1, r2) = props["code"], props["rows"]
    #    results = process_orf_array(info_matrix, arr, r1, r2, code)
    #    all_results.extend(results)

    for row_code, left, right in all_results:
        seq = extract_protein_seq(info_matrix, a_map[row_code], left, right)
        print("Protein:", seq)

    if all_results:
        result_matrix = np.array(all_results)
        print("\nUnique boundaries matrix (columns: ORF_code, left boundary, right boundary):")
        print(result_matrix)
        result_matrix_path = '/Users/benumlauf/Desktop/JS/pyfolder/mutated_ORFs'
        np.save(result_matrix_path, result_matrix)
        print("mutated ORFs successfully saved at", result_matrix_path)
    else:
        print("No boundaries computed for any ORF.")


    process_result_matrix(info_matrix, result_matrix, do_blast=True)
    return {
        "matches": matches,
        "chunk_matrices": chunk_matrices,
        "final_matrix": final_matrix,
        "interactive_matrix": interactive_matrix,
        "match_map": match_map,
        "result_matrix": result_matrix,
        "info_matrix": info_matrix,
        # and preserve diag info if you like:
        "diag_pct_met": pct_met,
        "diag_max_consec_fail": max_consec_fail,
    }

if __name__ =="__main__":
    #info_matrix = main()


    results = main()
#    final_matrix = main()


# /Users/benumlauf/Downloads/contig_1_Genome_R0_FLIPPED_Red_FASTA(1).fa
# /Users/benumlauf/Downloads/contig_1_R75_Genome_FASTA(1).fa
# /Users/benumlauf/Downloads/DNA_pol_R75(2).fa
#result = main()                       # result is your dict
      # this is the actual NumPy array
    # e.g. (<class 'numpy.ndarray'>, (8, 178))














# Connect the pick event (instead of a general click event).


###########################
#### END CODE BLOCK 2#####
##########################





#final_matrix = np.load('/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/final_matrix.npy', allow_pickle=True)
# Define the standard genetic code (using single-letter amino acid codes)



# --- Main processing ---

# Assuming final_matrix is a NumPy array where each row is defined as follows:
# Rows 0-7 (original data):
#   0: Alignment scores for each chunk (from Smith-Waterman)
#   1: Reference nucleotide sequence (list of characters)
#   2: Nucleotide numbering for the reference sequence
#   3: Evolved nucleotide sequence (list of characters)
#   4: Nucleotide numbering for the evolved sequence
#   5: Frame shift counter (increments for deletions, decrements for insertions)
#   6: Labels (1 for aligned, 0 for insertion/deletion, -1 for point mutation)
#   7: Count of total elements per row



########################################################
#looking for AA mutations in forward sequence###########
########################################################



# Initialize an empty list for storing the desired elements.


# ForwardAAMuts now contains every final_matrix[4, n] where final_matrix[8, n] and final_matrix[9, n] differ.


########################################
#######GRAPHING AA SEQUENCES############
########################################


# Define the mapping from amino acid letters to numeric values in alphabetical order.

#####################################
#       ORFIPY LFGGGGGGGGGG        #
####################################







# --- Define Sequences of Interest ---
# Mapping: name -> row index in loaded_matrix






# 0), Alignment data (chunkwise)
# 1), Reference nucleotide sequence
# 2), Reference nucleotide index
# 3), Evolvded nucleotide sequence
# 4), Evolved nucleotide index
# 5), Frameshift counter (del: +1, ins: -1)
# 6), Mut Type (in/del: 0, point: -1)
# 7), Global index
# 8), Amino acids in ref (naive ORF def)
# 9), Amino acids in evo (naive ORF def)
# 10), RevComp of Reference seq
# 11), Amino acids in RevComp of ref (naive ORF def)
# 12), RevComp of Evolved seq
# 13), Amino acids in RevComp of evo (naive ORF def)
# 14), Reference ORF frame 1
# 15), Reference ORF frame 2
# 16), Reference ORF frame 3
# 17), Evolved ORF frame 1
# 18), Evolved ORF frame 2
# 19), Evolved ORF frame 3
# 20), Reverse complement reference ORF frame 1
# 21), Reverse complement reference ORF frame 2
# 22), Reverse complement reference ORF frame 3
# 23), Reverse complement evolved ORF frame 1
# 24), Reverse complement evolved ORF frame 2
# 25), Reverse complement evolved ORF frame 3



# /Users/benumlauf/Downloads/contig_1_Genome_R0_FLIPPED_Red_FASTA(1).fa
# /Users/benumlauf/Downloads/contig_1_R75_Genome_FASTA(1).fa
# /Users/benumlauf/Downloads/DNA_pol_R75(2).fa

#actgatcggggctagcttagcgattagcgtactgactatacacgagccacactctctcttttgaggggaagcgatggattgcgatcgcgcagtcgatgcagggggacgatatctatagcggcgatgctgattctgatgcttagtgctagctgcgcggatggctatagcgcgtaatttaggcgaaagcgcgcgcgcgctagc

#actcatcggggcttagcttagcgattagcgtactgactatacacgagccacactctctcttttgaagtggaatcgatggattgcgatcgcgcaaatggtcgatgcagggggacgatatctatagcggcgagattctgatgcttagtgctagctgcgcggatggctatagcgcgtaatttaggcgaaagcgcgcgcgcgctagc