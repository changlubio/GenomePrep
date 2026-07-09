import subprocess
import sys
import pytest


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
