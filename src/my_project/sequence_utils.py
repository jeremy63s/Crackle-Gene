import parasail
import string
#import tkinter as tk
ALPHABET = string.ascii_uppercase
substitution_matrix = parasail.matrix_create(ALPHABET, 1, -1)
scatter = None
plotted_indices = None
data = None
#from tkinter import scrolledtext

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
    #root = tk._default_root
   #if root is None:
       # root = tk.Tk()

    # Create a new window for the text display.
   # win = tk.Toplevel(root)
    win.title(f"Data Matrix Slice (Columns {start_idx} to {end_idx-1})")

    # Create a ScrolledText widget with no wrapping and a larger monospaced font.
    text_widget = scrolledtext.ScrolledText(win, wrap='none', font=("Courier", 14))
   # text_widget.insert(tk.END, subset_str)
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