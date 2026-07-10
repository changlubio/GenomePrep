from __future__ import annotations

import argparse
import subprocess
import sys
from unittest.mock import patch

import pytest

from genomeprep._types import ProcessResult


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "genomeprep.cli", *args],
        capture_output=True,
        text=True,
    )


class TestCliHelp:
    def test_help_exits_cleanly(self):
        result = run_cli("--help")
        assert result.returncode == 0
        assert "genomeprep" in result.stdout.lower()

    def test_process_subcommand_help(self):
        result = run_cli("process", "--help")
        assert result.returncode == 0
        assert "infile" in result.stdout

    def test_download_data_subcommand_help(self):
        result = run_cli("download-data", "--help")
        assert result.returncode == 0

    def test_no_args_prints_help(self):
        result = run_cli()
        assert result.returncode != 0 or "usage" in result.stdout.lower()

    def test_process_missing_file_exits_nonzero(self):
        result = run_cli("process", "/nonexistent/file.txt")
        assert result.returncode != 0


class TestProcessCmdDirect:
    def test_file_not_found_returns_1(self, tmp_path):
        from genomeprep.cli import _process_cmd

        args = argparse.Namespace(
            infile=str(tmp_path / "missing.txt"),
            datadir=None,
            outdir=str(tmp_path),
            outindex="out",
        )
        assert _process_cmd(args) == 1

    def test_processed_prints_success_and_returns_0(self, tmp_path):
        from genomeprep.cli import _process_cmd

        f = tmp_path / "genome.txt"
        f.write_text("x")
        args = argparse.Namespace(
            infile=str(f),
            datadir=None,
            outdir=str(tmp_path),
            outindex="out",
        )
        mock_result = ProcessResult(status="Processed", info={"snps": 500000, "assembly": "GRCh37"})
        with patch("genomeprep.process", return_value=mock_result):
            assert _process_cmd(args) == 0

    def test_invalid_status_returns_1(self, tmp_path):
        from genomeprep.cli import _process_cmd

        f = tmp_path / "genome.txt"
        f.write_text("x")
        args = argparse.Namespace(
            infile=str(f),
            datadir=None,
            outdir=str(tmp_path),
            outindex="out",
        )
        mock_result = ProcessResult(status="Invalid", info="bad file")
        with patch("genomeprep.process", return_value=mock_result):
            assert _process_cmd(args) == 1

    def test_processed_with_string_info(self, tmp_path):
        from genomeprep.cli import _process_cmd

        f = tmp_path / "genome.txt"
        f.write_text("x")
        args = argparse.Namespace(
            infile=str(f),
            datadir=None,
            outdir=str(tmp_path),
            outindex="out",
        )
        # result.status == "Processed" but info is a string (edge case)
        mock_result = ProcessResult(status="Processed", info="ok")
        with patch("genomeprep.process", return_value=mock_result):
            assert _process_cmd(args) == 0


class TestDownloadCmdDirect:
    def test_calls_download_and_returns_0(self, tmp_path):
        from genomeprep.cli import _download_cmd

        args = argparse.Namespace(dest=str(tmp_path))
        with patch("genomeprep.data.download", return_value=str(tmp_path)) as mock_dl:
            result = _download_cmd(args)
        assert result == 0
        mock_dl.assert_called_once_with(dest=str(tmp_path))

    def test_download_with_none_dest(self):
        from genomeprep.cli import _download_cmd

        args = argparse.Namespace(dest=None)
        with patch("genomeprep.data.download", return_value="/tmp/genomeprep"):
            result = _download_cmd(args)
        assert result == 0


class TestMainDirect:
    def test_main_process_subcommand_missing_file(self, tmp_path):
        from genomeprep.cli import main

        with patch("sys.argv", ["genomeprep", "process", str(tmp_path / "missing.txt")]):
            with pytest.raises(SystemExit) as exc:
                main()
            assert exc.value.code == 1

    def test_main_no_command_exits_nonzero(self):
        from genomeprep.cli import main

        with patch("sys.argv", ["genomeprep"]):
            with pytest.raises(SystemExit) as exc:
                main()
            assert exc.value.code != 0

    def test_main_download_data(self):
        from genomeprep.cli import main

        with patch("sys.argv", ["genomeprep", "download-data", "--dest", "/tmp/gp"]):
            with patch("genomeprep.data.download", return_value="/tmp/gp"):
                with pytest.raises(SystemExit) as exc:
                    main()
                assert exc.value.code == 0
