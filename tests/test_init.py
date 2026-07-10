from __future__ import annotations

import io
from unittest.mock import MagicMock, patch


class TestProcess:
    """Tests for genomeprep.process() — all heavy dependencies are mocked."""

    def _mock_ref(self):
        ref = MagicMock()
        ref.api_23andme = set()
        ref.the_list = {}
        ref.reverse_snps = {}
        ref.bad_allele_clusters = {}
        ref.faidx = MagicMock()
        return ref

    def _make_reject_tuple(self, reject_type: str, count_snps: int = 5000):
        return (
            True,
            reject_type,
            [],
            count_snps,
            [],
            0,
            0,
            -1.0,
            -1.0,
            -1,
            -1.0,
            -1.0,
            -1.0,
            -1,
            -1,
        )

    def test_file_not_resolved(self, tmp_path):
        from genomeprep import process

        with patch("genomeprep.resolve_file_for_genome", return_value=(False, None, False)):
            result = process(str(tmp_path / "x.txt"))
        assert result.status == "Invalid"
        assert "Cannot identify" in result.info

    def test_file_not_resolved_from_zip(self, tmp_path):
        from genomeprep import process

        with patch("genomeprep.resolve_file_for_genome", return_value=(False, None, True)):
            result = process(str(tmp_path / "x.zip"))
        assert result.status == "Invalid"
        assert "ZIP" in result.info

    def test_handle_none_after_resolve(self, tmp_path):
        from genomeprep import process

        with patch("genomeprep.resolve_file_for_genome", return_value=(True, None, False)):
            result = process(str(tmp_path / "x.txt"))
        assert result.status == "error"

    def test_parser_rejects_grch(self, tmp_path):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = True
        mock_parser.reject_type = "grch38"
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                result = process("x.txt")
        assert result.status == "Invalid"
        assert "GRCh37" in result.info

    def test_parser_rejects_other(self, tmp_path):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = True
        mock_parser.reject_type = "invalid"
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                result = process("x.txt")
        assert result.status == "Invalid"

    def test_dtc_read_records_reject_vcf(self, tmp_path):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = False
        mock_parser.total_badsnp = 0
        ref = self._mock_ref()
        reject_tuple = self._make_reject_tuple("vcf", count_snps=100)
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                with patch("genomeprep.get_data_dir", return_value="/fake"):
                    with patch("genomeprep.load_reference_data", return_value=ref):
                        with patch("genomeprep.dtc_read_records", return_value=reject_tuple):
                            result = process("x.txt")
        assert result.status == "vcf"

    def test_dtc_read_records_reject_invalid_zero_snps(self):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = False
        mock_parser.total_badsnp = 0
        ref = self._mock_ref()
        reject_tuple = self._make_reject_tuple("invalid", count_snps=0)
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                with patch("genomeprep.get_data_dir", return_value="/fake"):
                    with patch("genomeprep.load_reference_data", return_value=ref):
                        with patch("genomeprep.dtc_read_records", return_value=reject_tuple):
                            result = process("x.txt")
        assert result.status == "Invalid"

    def test_dtc_read_records_reject_other_type(self):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = False
        mock_parser.total_badsnp = 5
        ref = self._mock_ref()
        reject_tuple = self._make_reject_tuple("incomplete", count_snps=5000)
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                with patch("genomeprep.get_data_dir", return_value="/fake"):
                    with patch("genomeprep.load_reference_data", return_value=ref):
                        with patch("genomeprep.dtc_read_records", return_value=reject_tuple):
                            result = process("x.txt")
        assert result.status == "incomplete"

    def test_happy_path_returns_processed(self, tmp_path):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = False
        mock_parser.total_badsnp = 10
        mock_parser.company = "23andme"
        mock_parser.build = "GRCh37"
        ref = self._mock_ref()
        store_snps = [
            {
                "rsid": "rs1",
                "chro": "1",
                "pos": "100",
                "call": "AA",
                "ref": "A",
                "gav": "c1",
            }
        ]
        ok_tuple = (
            False,
            "",
            store_snps,
            500000,
            [],
            100,
            450000,
            0.95,
            0.02,
            5,
            0.9,
            0.03,
            0.01,
            "cluster1",
            0.85,
        )
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                with patch("genomeprep.get_data_dir", return_value="/fake"):
                    with patch("genomeprep.load_reference_data", return_value=ref):
                        with patch("genomeprep.dtc_read_records", return_value=ok_tuple):
                            with patch("genomeprep.filter_by_badalleles", return_value=store_snps):
                                with patch("genomeprep.output_as_23andme", return_value=450000):
                                    with patch(
                                        "genomeprep.directconvertSNPtoVCFformat",
                                        return_value=(450000, 100, 50, 20),
                                    ):
                                        result = process("x.txt", outdir=str(tmp_path))
        assert result.status == "Processed"
        assert isinstance(result.info, dict)
        assert "final_snps" in result.info

    def test_reject_with_integer_cluster_id(self):
        from genomeprep import process

        mock_handle = io.StringIO("data")
        mock_parser = MagicMock()
        mock_parser.reject = False
        mock_parser.total_badsnp = 3
        ref = self._mock_ref()
        # cluster_id is integer (not str), to exercise the str(cluster_id) branch
        reject_tuple = self._make_reject_tuple("incomplete", count_snps=5000)
        with patch("genomeprep.resolve_file_for_genome", return_value=(True, mock_handle, False)):
            with patch("genomeprep.DTCparser", return_value=mock_parser):
                with patch("genomeprep.get_data_dir", return_value="/fake"):
                    with patch("genomeprep.load_reference_data", return_value=ref):
                        with patch("genomeprep.dtc_read_records", return_value=reject_tuple):
                            result = process("x.txt")
        assert result.status == "incomplete"
