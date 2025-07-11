import matplotlib

matplotlib.use('TkAgg')  # Force TkAgg backend
from Bio.Blast import NCBIWWW, NCBIXML

import parasail
import string
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm  # Standard tqdm for PyCharm
import os


import os
import subprocess
import sys
import platform
import orfipy
import orfipy.findorfs as of
import orfipy_core
from Bio.Seq import Seq

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