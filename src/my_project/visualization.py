import numpy as np
import tqdm
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import scrolledtext
root = tk.Tk()
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
        else:
            print(f"No match found for cell (Evolved {row + 1}, Reference {col + 1})")

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


