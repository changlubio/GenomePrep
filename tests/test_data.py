from __future__ import annotations

import os
import pickle
from unittest.mock import MagicMock, patch

import pytest

from genomeprep.data import ReferenceData, get_data_dir


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


class TestDownload:
    def test_download_creates_dir_and_fetches_files(self, tmp_path):
        from genomeprep.data import _DATA_FILES, download

        dest = str(tmp_path / "data")
        with patch("urllib.request.urlretrieve") as mock_fetch:
            result = download(dest=dest)

        assert result == dest
        assert os.path.isdir(dest)
        assert mock_fetch.call_count == len(_DATA_FILES)

    def test_download_skips_existing_files(self, tmp_path):
        from genomeprep.data import _DATA_FILES, download

        dest = str(tmp_path)
        # Pre-create all files
        for fname in _DATA_FILES:
            (tmp_path / fname).write_text("existing")

        with patch("urllib.request.urlretrieve") as mock_fetch:
            download(dest=dest)

        assert mock_fetch.call_count == 0  # nothing downloaded

    def test_download_with_custom_base_url(self, tmp_path):
        from genomeprep.data import download

        with patch("urllib.request.urlretrieve") as mock_fetch:
            download(dest=str(tmp_path), base_url="https://example.com/data/")

        first_call_url = mock_fetch.call_args_list[0][0][0]
        assert first_call_url.startswith("https://example.com/data/")


class TestReadApi:
    def test_parses_api_file(self, tmp_path):
        from genomeprep.data import _read_23andme_api

        api_file = tmp_path / "api.23andme.com"
        api_file.write_text("# comment\nsome data 1 100\nsome data 2 200\n")
        result = _read_23andme_api(str(api_file))
        assert "1:100" in result
        assert "2:200" in result

    def test_skips_short_lines(self, tmp_path):
        from genomeprep.data import _read_23andme_api

        api_file = tmp_path / "api.23andme.com"
        api_file.write_text("# comment\nonly_one_part\n")
        result = _read_23andme_api(str(api_file))
        assert len(result) == 0


class TestLoadReferenceDataSuccess:
    def test_loads_all_files(self, tmp_path):
        from genomeprep.data import load_reference_data

        # Create all required files
        the_list = {"1:100": ["cluster1"]}
        bad_alleles = {"cluster1": {"1:200"}}
        rev_snps = {"1:300": True}

        with open(tmp_path / "THE_LIST.dat", "wb") as fh:
            pickle.dump(the_list, fh)
        with open(tmp_path / "badalleles.dat", "wb") as fh:
            pickle.dump(bad_alleles, fh)
        with open(tmp_path / "RS2GRCh37Orien_1.dat", "wb") as fh:
            pickle.dump(rev_snps, fh)
        (tmp_path / "api.23andme.com").write_text("# comment\ndata col 1 100\n")
        (tmp_path / "Homo_sapiens.GRCh37.75.dna.toplevel.fa").write_text(">1\nACGT\n")

        with patch("pyfaidx.Fasta") as mock_fasta_cls:
            mock_fasta = MagicMock()
            mock_fasta_cls.return_value = mock_fasta
            ref = load_reference_data(str(tmp_path))

        assert ref.the_list == the_list
        assert ref.bad_allele_clusters == bad_alleles
        assert ref.reverse_snps == rev_snps
        assert "1:100" in ref.api_23andme
