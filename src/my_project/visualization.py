import numpy as np
import tqdm
import matplotlib.pyplot as plt
#import tkinter as tk
#from tkinter import scrolledtext
#root = tk.Tk()
# src/my_project/visualization.py
# Streamlit-safe plotting (no tkinter, no plt.show)
import numpy as np
import matplotlib.pyplot as plt

def plot_interactive_matrix(num_chunks, matches, match_map, threshold):
    """
    Build a simple heatmap of normalized scores.
    - No tkinter
    - No plt.show()
    - Returns a Matplotlib Figure (call st.pyplot(fig) in Streamlit)
    """
    # initialize with NaNs so "missing" cells render light/white
    mat = np.full((num_chunks, num_chunks), np.nan, dtype=float)
    for m in matches:
        i = int(m["evo_index"]) - 1
        j = int(m["ref_index"]) - 1
        # normalized_score expected in [0,1]
        mat[i, j] = float(m.get("normalized_score", np.nan))

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(mat, origin="lower", aspect="auto", vmin=0.0, vmax=1.0)
    ax.set_xlabel("Reference chunk")
    ax.set_ylabel("Evolved chunk")
    ax.set_title(f"Match matrix (threshold ≥ {threshold:.2f})")
    cbar = fig.colorbar(im, ax=ax, label="Normalized score")
    return fig

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


