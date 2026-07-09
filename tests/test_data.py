from __future__ import annotations

import os
import pickle
import tempfile
from unittest.mock import MagicMock, patch
import pytest

from genomeprep.data import get_data_dir, ReferenceData


class TestGetDataDir:
    def test_returns_explicit_arg(self, tmp_path):
        result = get_data_dir(str(tmp_path))
        assert result == str(tmp_path)

    def test_reads_env_var(self, tmp_path, monkeypatch):
        monkeypatch.setenv("GENOMEPREP_DATA_DIR", str(tmp_path))
        result = get_data_dir(None)
        assert result == str(tmp_path)

    def test_explicit_arg_beats_env_var(self, tmp_path, monkeypatch):
        other = tmp_path / "other"
        other.mkdir()
        monkeypatch.setenv("GENOMEPREP_DATA_DIR", str(tmp_path))
        result = get_data_dir(str(other))
        assert result == str(other)

    def test_falls_back_to_cache(self, monkeypatch):
        monkeypatch.delenv("GENOMEPREP_DATA_DIR", raising=False)
        result = get_data_dir(None)
        assert "genomeprep" in result


class TestReferenceDataFields:
    def test_dataclass_has_expected_fields(self):
        import dataclasses
        fields = {f.name for f in dataclasses.fields(ReferenceData)}
        assert fields == {"api_23andme", "the_list", "reverse_snps", "bad_allele_clusters", "faidx"}


class TestLoadReferenceData:
    def test_raises_file_not_found_for_missing_the_list(self, tmp_path):
        from genomeprep.data import load_reference_data
        with pytest.raises(FileNotFoundError, match="THE_LIST.dat"):
            load_reference_data(str(tmp_path))

    def test_raises_file_not_found_for_missing_fasta(self, tmp_path):
        from genomeprep.data import load_reference_data
        # Create all non-FASTA files so error points to FASTA
        for fname in ["THE_LIST.dat", "badalleles.dat", "RS2GRCh37Orien_1.dat"]:
            with open(tmp_path / fname, "wb") as fh:
                pickle.dump({}, fh)
        (tmp_path / "api.23andme.com").write_text("# header\n")
        with pytest.raises(FileNotFoundError, match="GRCh37"):
            load_reference_data(str(tmp_path))
