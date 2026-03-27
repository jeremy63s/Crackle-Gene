$ErrorActionPreference = "Stop"

# repo root (this script is in scripts/)
Set-Location -Path (Join-Path $PSScriptRoot "..")

# Prefer Python 3.11
$py = "py -3.11"
try {
  $v = (& py -3.11 --version) 2>$null
} catch {
  $py = "python"
}

# Create venv if missing
if (-not (Test-Path ".venv")) {
  & $py -m venv .venv
}

# Activate venv
& .\.venv\Scripts\Activate.ps1

# Upgrade tooling and install deps
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install .

# Headless-friendly default
$env:MPLBACKEND = $env:MPLBACKEND -ne $null ? $env:MPLBACKEND : "Agg"

$env:PYTHONPATH = (Join-Path (Get-Location) "src")

# Run Streamlit
python -m streamlit run src\my_project\MAINAPP.py
