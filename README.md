# crackle_gene

A Genome Sequence Comparator with Parasail, Biopython, ORFipy, and a Streamlit UI.


## Quickstart (macOS/Linux)
```bash
git clone https://github.com/jeremy63s/Crackle-Gene.git
cd Crackle-Gene-main
bash scripts/setup_and_run.sh

## Powershell
git clone https://github.com/jeremy63s/Crackle-Gene.git
cd Crackle-Gene-main
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\scripts\setup_and_run.ps1


## Docker
docker build -t crackle-gene .
docker run --rm -p 8501:8501 crackle-gene
# Visit http://localhost:8501
