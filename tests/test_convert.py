from genomeprep._types import SNPRecord
from genomeprep.convert import directconvertSNPtoVCFformat, flip_genotype, output_as_23andme


def _snp(
    rsid: str = "rs1",
    chro: str = "1",
    pos: str = "100",
    call: str = "AG",
    ref: str = "A",
    gav: str = "",
) -> SNPRecord:
    return SNPRecord(rsid=rsid, chro=chro, pos=pos, call=call, ref=ref, gav=gav, ba=False)


class TestFlipGenotype:
    def test_flips_all_bases(self):
        assert flip_genotype("ATCG") == "TAGC"

    def test_empty_string(self):
        assert flip_genotype("") == ""

    def test_single_base(self):
        assert flip_genotype("A") == "T"
        assert flip_genotype("C") == "G"


class TestOutputAs23andme:
    def test_writes_snps_and_returns_count(self, tmp_path):
        snps: list[SNPRecord] = [
            _snp(rsid="rs1", chro="1", pos="100", call="AA"),
            _snp(rsid="rs2", chro="1", pos="200", call="GG"),
        ]
        out = str(tmp_path / "out.23andme")
        count = output_as_23andme(snps, "##header", out)
        assert count == 2
        lines = open(out).readlines()
        assert lines[0].strip() == "##header"
        assert "rs1\t1\t100\tAA" in lines[1]

    def test_empty_snps(self, tmp_path):
        out = str(tmp_path / "out.23andme")
        count = output_as_23andme([], "##header", out)
        assert count == 0


class TestDirectConvertSNPtoVCF:
    def test_writes_vcf_and_returns_counts(self, tmp_path):
        snps: list[SNPRecord] = [
            _snp(rsid="rs1", chro="1", pos="100", call="AG", ref="A"),
        ]
        out = str(tmp_path / "out.vcf")
        final, indels, refref, multi = directconvertSNPtoVCFformat(snps, "##header", "sample1", out)
        assert final == 1
        assert indels == 0
        assert refref == 0
        assert multi == 0
        content = open(out).read()
        assert "rs1" in content
        assert "0/1" in content or "1/0" in content
        assert "GAV=none" in content
        assert "PASS" in content

    def test_emits_indels_with_tag(self, tmp_path):
        snps: list[SNPRecord] = [_snp(call="ID", ref="A")]
        out = str(tmp_path / "out.vcf")
        final, indels, refref, multi = directconvertSNPtoVCFformat(snps, "##header", "s", out)
        assert final == 1  # indels ARE emitted
        assert indels == 1
        content = open(out).read()
        assert "INDEL=ID" in content
        assert "indel" in content

    def test_emits_refref_with_dot_alt(self, tmp_path):
        snps: list[SNPRecord] = [_snp(call="AA", ref="A")]
        out = str(tmp_path / "out.vcf")
        final, indels, refref, multi = directconvertSNPtoVCFformat(snps, "##header", "s", out)
        assert final == 1  # refref IS emitted
        assert refref == 1
        content = open(out).read()
        assert "\t.\t" in content  # alt is "."

    def test_badallele_adds_bz3_filter(self, tmp_path):
        snp = _snp(call="AG", ref="A")
        snp["ba"] = True
        out = str(tmp_path / "out.vcf")
        final, indels, refref, multi = directconvertSNPtoVCFformat([snp], "##header", "s", out)
        assert final == 1
        content = open(out).read()
        assert "bz3" in content

    def test_gav_annotation_in_info(self, tmp_path):
        snp_with_gav = _snp(call="AG", ref="A", gav="cluster1")
        snp_no_gav = _snp(rsid="rs2", pos="200", call="AG", ref="A", gav="")
        out = str(tmp_path / "out.vcf")
        directconvertSNPtoVCFformat([snp_with_gav, snp_no_gav], "##header", "s", out)
        content = open(out).read()
        assert "GAV=cluster1" in content
        assert "GAV=none" in content

    def test_gt_separator_is_slash(self, tmp_path):
        snp = _snp(call="AG", ref="A")
        out = str(tmp_path / "out.vcf")
        directconvertSNPtoVCFformat([snp], "##header", "s", out)
        content = open(out).read()
        assert "0/1" in content or "1/0" in content
