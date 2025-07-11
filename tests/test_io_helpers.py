import tempfile
from my_project.io_helpers import read_fasta

def test_read_fasta(tmp_path):
    f = tmp_path / "test.fa"
    f.write_text(">hdr\nATG\nCGA\n")
    assert read_fasta(str(f)) == "ATGCGA"
