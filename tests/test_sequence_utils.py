import pytest
from my_project.sequence_utils import rotate_sequence

def test_rotate_sequence_simple():
    assert rotate_sequence("ABCDE", "CD") == "CDEAB"
