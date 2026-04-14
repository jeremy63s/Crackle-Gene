# Crackle-Gene

A bioinformatics tool for comparing two DNA sequences — a reference and an evolved variant — to identify and characterize mutations at the nucleotide and amino acid level.

The pipeline aligns chunked sequence sections using [parasail](https://github.com/jeffdaily/parasail), predicts open reading frames (ORFs) using [Prodigal](https://github.com/hyattpd/Prodigal), tracks active ORF frames across the alignment, and surfaces mutations with their amino-acid consequences. Results are presented in an interactive [Streamlit](https://streamlit.io) web UI that supports BLAST submission for each affected ORF.

## Features

- Chunk-based pairwise alignment with diagonal pre-screening
- Forward and reverse-complement ORF tracking across all three reading frames
- Nucleotide → amino acid translation for reference and evolved sequences
- Mutation site detection with ORF-aware pairing and scoring
- BLAST integration: submit individual ORF sequences directly from the UI
- PDF-friendly sequence comparison view with syntax highlighting
- Copy-all button for exporting selected mutation summaries

## Requirements

- Python 3.11+
- [Prodigal](https://github.com/hyattpd/Prodigal) installed and on `PATH` (e.g. `brew install prodigal` on macOS)
- Python packages listed in `requirements.txt`

## Quickstart (macOS / Linux / WSL)

```bash
git clone https://github.com/jeremy63s/Crackle-Gene.git
cd Crackle-Gene
bash scripts/setup_and_run.sh
```

The script creates a virtual environment, installs dependencies, and launches the Streamlit app at `http://localhost:8501`.


```

## Usage

1. **Upload sequences** — provide a reference FASTA and an evolved FASTA (upload or paste).
2. **Run alignment** — the pipeline aligns and builds the info matrix.
3. **Submit to BLAST** — select mutation sites and submit individual ORF sequences to NCBI BLAST.
4. **View results** — browse sequence comparisons and BLAST results per mutation site.

## Project structure

```
src/my_project/
  MAINAPP.py        # Streamlit UI and page routing
  pypeline.py       # Pipeline orchestration
  microservices.py  # ORF tracking, diff detection, result matrix construction
  alignment.py      # Parasail-based pairwise alignment
  sequence_utils.py # Sequence manipulation utilities
  blast.py          # NCBI BLAST submission and result parsing
  io_helpers.py     # FASTA I/O helpers
test_data/          # Example E. coli reference and query sequences
tests/              # Unit tests
```

## License

MIT — see [LICENSE](LICENSE).
