import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys
import platform
import orfipy
import orfipy.findorfs as of
import orfipy_core
from Bio.Seq import Seq
import parasail
from sequence_utils import substitution_matrix
from tqdm import tqdm            # ← add this!
from typing import Tuple, List
import parasail
import string
ALPHABET = string.ascii_uppercase
substitution_matrix = parasail.matrix_create(ALPHABET, 1, -1)

from my_project.io_helpers import read_fasta, get_sequence, extract_extended_context
from my_project.sequence_utils import rotate_sequence, divide_sequence_with_offset
from my_project.alignment import diagonal_check, compare_sections, get_alignment
from my_project.matrix_utils import build_chunk_matrix, build_final_matrix
from my_project.blast import blast_protein_seq, process_result_matrix,  extract_protein_seq
from my_project.visualization import plot_interactive_matrix
from sequence_utils import on_pick, CODON_TABLE, translate_nuc_row_active, reverse_complement, create_annotation_array
from sequence_utils import array_to_str, parse_orf_frame,  process_orf_array

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
    #print("aligned_evo length is: ", len(aligned_evo), "aligned ref length is: ", len(aligned_ref))
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
    #final_matrix_path = "/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/final_matrix.npy"
    #np.save(final_matrix_path, final_matrix)
    #if os.path.exists(final_matrix_path):
    #    print("Final matrix saved successfully at:", final_matrix_path)
    #else:
    #    print("Error: Final matrix file was not saved.")

    # Also, save the final matrix as a text file (optional).
    # final_matrix_txt_path = "/Users/benumlauf/PycharmProjects/GENOMEANALYSIS/final_matrix.txt"
    # np.savetxt(final_matrix_txt_path, final_matrix, fmt='%s', delimiter='\t')
    # if os.path.exists(final_matrix_txt_path):
    #    print("Final matrix text file saved successfully at:", final_matrix_txt_path)
    # else:
    #    print("Error: Final matrix text file was not saved.")
    #final_matrix = results["final_matrix"]
    #print(type(final_matrix), final_matrix.shape)
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

    #np.save("/Users/benumlauf/Downloads/final_matrix14.npy", loaded_matrix)

    #############################
    # TURN .GBK INTO ARRAY     #
    ############################
    from Bio import SeqIO

    # Use the GenBank file exported from SnapGene that includes ORF annotations.
    file_path = "/Users/siegelmanfamily/Downloads/crackle_gene/GBK.gb"
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


    global_path_info = '/Users/siegelmanfamily/Desktop/JS/pyproject/info_matrix.npy'
    np.save(global_path_info, info_matrix)

    # SAVE MATRIX CHEATSHEET

    #global_path_cheatsheet = '/Users/benumlauf/Desktop/JS/pyfolder/info_matrix_cheatsheet.npy'
    #np.save(global_path_cheatsheet, info_matrix_cheatsheet)

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

    result_matrix_path = "/Users/siegelmanfamily/Desktop/JS/pyproject/result_matrix"
#chatgpt told me to do this######
    # --- after you build all_results ---
    if all_results:
        result_matrix = np.array(all_results)
        np.save(result_matrix_path, result_matrix)
        print("mutated ORFs successfully saved at", result_matrix_path)
    else:
        print("No boundaries computed for any ORF.")
        result_matrix = np.empty((0, 3), dtype=int)

    process_result_matrix(info_matrix, result_matrix, do_blast=False)
    #return { … }

#    if all_results:
#        result_matrix = np.array(all_results)
#        print("\nUnique boundaries matrix (columns: ORF_code, left boundary, right boundary):")
#        print(result_matrix)
#        result_matrix_path = '/Users/siegelmanfamily/Desktop/JS/pyproject'
#        np.save(result_matrix_path, result_matrix)
#        print("mutated ORFs successfully saved at", result_matrix_path)
#    else:
#        print("No boundaries computed for any ORF.")


    process_result_matrix(info_matrix, result_matrix, do_blast=False)
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
#actgatcggggctagcttagcgattagcgtactgactatacacgagccacactctctcttttgaggggaagcgatggattgcgatcgcgcagtcgatgcagggggacgatatctatagcggcgatgctgattctgatgcttagtgctagctgcgcggatggctatagcgcgtaatttaggcgaaagcgcgcgcgcgctagc

# /Users/siegelmanfamily/Downloads/contig_1_Genome_R0_FLIPPED_Red_FASTA.fa
#/Users/siegelmanfamily/Downloads/contig_1_R75_Genome_FASTA.fa
#/Users/siegelmanfamily/Desktop/JS/genome_stuff/Sf_e293f_R52_genome_FASTA.fa
#/Users/siegelmanfamily/Downloads/contig_1_Genome_R0_FASTA.fa

# RNA polymerase subunit alpha:
# ttactcgtcagcgatgcttgccggtggccagttttccaagcgcatgcccagagacagtccacgggaagccagcacgtctttaatctcagtaagagattttttaccaaggttaggcgttttaaggagctcaacctcggtacgctgtaccagatcaccgatatagtggatagcttctgctttaaggcagttagcagagcggacagtcaattccagatcgtcaacagggcgcagcaggatcggatcgaactctggtttctcttctttcacttcaggctgacgtacatcacgtaagtcaacgaaagcttccagttgttcagccagaatggttgccgcacgacgaatcgc