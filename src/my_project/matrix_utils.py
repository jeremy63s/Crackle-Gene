import numpy as np

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