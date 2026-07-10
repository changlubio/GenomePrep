import gzip
import io
import os
import zipfile

import pytest

from genomeprep.io import check_for_vcf, resolve_file_for_genome

TUTORIAL = os.path.join(os.path.dirname(__file__), "..", "tutorial")


class TestResolveFileForGenome:
    def test_resolves_plain_text_file(self, tmp_path):
        f = tmp_path / "genome.txt"
        f.write_text("# header\nrs1\t1\t100\tAA\n")
        resolved, handle, iszip = resolve_file_for_genome(str(f))
        assert resolved
        assert handle is not None
        assert not iszip
        handle.close()

    def test_returns_false_for_nonexistent_file(self, tmp_path):
        resolved, handle, iszip = resolve_file_for_genome(str(tmp_path / "missing.txt"))
        assert not resolved
        assert handle is None

    def test_resolves_gzip_file(self, tmp_path):
        f = tmp_path / "genome.txt.gz"
        with gzip.open(str(f), "wt") as gh:
            gh.write("# header\nrs1\t1\t100\tAA\n")
        resolved, handle, iszip = resolve_file_for_genome(str(f))
        assert resolved
        assert handle is not None
        assert not iszip
        content = handle.read()
        assert "rs1" in content
        handle.close()

    def test_resolves_zip_containing_text(self, tmp_path):
        inner = tmp_path / "genome.txt"
        inner.write_text("# header\nrs1\t1\t100\tAA\n")
        zf = tmp_path / "genome.zip"
        with zipfile.ZipFile(str(zf), "w") as z:
            z.write(str(inner), arcname="genome.txt")
        resolved, handle, iszip = resolve_file_for_genome(str(zf))
        assert resolved
        assert iszip
        assert handle is not None
        handle.close()

    def test_tutorial_testgenome_zip(self):
        zpath = os.path.join(TUTORIAL, "testgenome.zip")
        if not os.path.exists(zpath):
            pytest.skip("tutorial/testgenome.zip not present")
        resolved, handle, iszip = resolve_file_for_genome(zpath)
        assert resolved
        assert iszip
        if handle:
            handle.close()


class TestCheckForVcf:
    def test_returns_false_for_dtc_format(self):
        content = "# header\nrs1\t1\t100\tAA\n" * 100
        assert not check_for_vcf(io.StringIO(content))

    def test_returns_true_for_vcf_format(self):
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        row = "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0|1\n"
        content = header + row * 6000
        assert check_for_vcf(io.StringIO(content))

    def test_returns_false_for_none(self):
        assert not check_for_vcf(None)
