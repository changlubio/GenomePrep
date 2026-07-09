from __future__ import annotations

from collections.abc import Iterator
from typing import IO

from genomeprep._types import Company
from genomeprep.convert import flip_genotype

valid_call: frozenset[str] = frozenset({"A", "C", "G", "T", "I", "D", "_", "-", "0"})
valid_chromo: frozenset[str] = frozenset(
    {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "MT",
    }
)
COMPANY: list[str] = [
    "23andme",
    "ancestrydna",
    "ftdna",
    "myheritage",
    "decodeme",
    "livingdnacsv",
    "all",
]
reasons: dict[str, str] = {
    "invalid": (
        "We are unable to properly process your upload. Please check your file and try again."
    ),
    "small": "Your genome is too small (<10000 lines).",
    "incomplete": (
        "Your genome is incomplete. We only process complete genomes "
        "where all chromosomes are sequenced."
    ),
    "grch36": "Your DNA file is generated using human genome assembly GRCh36. ",
    "grch38": "Your DNA file is generated using human genome assembly GRCh38.",
    "notgrch37": ("save4check Your DNA file contains positions not from assembly GRCh37."),
    "vcf": "save4check This could be because it is an vcf file.",
    "zip": "save4check This could be because it is compressed/zipped.",
    "ref": "save4check This could be because it has unspecified reverse strand SNPs.",
    "imputation": "save4check This could be because it is imputation.",
}


class DTCparser:
    """Parse a direct-to-consumer genome file. Iterable; yields DTCrecord objects."""

    def __init__(self, handle: IO[str], filename: str = "") -> None:
        self.path = filename
        self.company: Company | None = None
        self.build = "GRCh37"
        self.meta_data: list[str] = []
        self.total_badsnp = 0
        self.reject = False
        self.reject_type: str = ""
        self.handle = handle
        self.norsid = 0
        self._buffered_line: str = ""
        try:
            self._read_headers()
        except Exception as exc:
            raise RuntimeError("Cannot read the file by handle") from exc

    def _read_headers(self) -> None:
        line: str = ""
        while True:
            try:
                raw = self.handle.readline()
            except UnicodeDecodeError:
                self.reject = True
                break
            line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
            if not line or not line.startswith("#"):
                break
            self.meta_data.append(line.rstrip())

        if not line:
            self.reject = True
            self.reject_type = "invalid"
            self.handle.close()
            return

        if not self.meta_data:
            self.meta_data.append(line.rstrip())
        else:
            # Buffer the first non-comment line so it is not lost during iteration
            self._buffered_line = line

        first = self.meta_data[0]
        if "23andme" in first or "23andMe" in first:
            self.company = "23andme"
            self._check_build_in_meta()
        elif "AncestryDNA" in first:
            self.company = "ancestrydna"
            self._check_build_in_meta()
        elif "RSID,CHROMOSOME,POSITION,RESULT" in first:
            self.company = "ftdna"
        elif "MyHeritage" in first:
            self.company = "myheritage"
        elif "Name,Variation,Chromosome,Position,Strand,YourCode" in first:
            self.company = "decodeme"
        elif "Living DNA customer genotype data download file" in first:
            self.company = "livingdnacsv"
            self._check_build_in_meta()
        else:
            self.company = "all"

    def _check_build_in_meta(self) -> None:
        for x in self.meta_data:
            if "build" in x:
                if "36" in x:
                    self.reject, self.reject_type, self.build = True, "grch36", "GRCh36"
                    return
                if "38" in x:
                    self.reject, self.reject_type, self.build = True, "grch38", "GRCh38"
                    return

    def __iter__(self) -> Iterator[DTCrecord]:
        return self

    def __next__(self) -> DTCrecord:
        if self._buffered_line:
            line = self._buffered_line
            self._buffered_line = ""
        else:
            try:
                raw = self.handle.readline()
            except Exception:
                raw = ""
            line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
        if not line:
            self.handle.close()
            raise StopIteration
        record = DTCrecord(line, self.company or "all")
        if not record.goodrsid:
            self.norsid += 1
        return record


class DTCrecord:
    """A single SNP record parsed from one line of a DTC file."""

    def __init__(self, line: str, company: Company) -> None:
        self.line = line.rstrip()
        self.rsid = self.chro = self.pos = self.genotype = ""
        self.badsnp = False
        getattr(self, f"_parse_{company}")()
        self._check_badsnp()
        self.loci = f"{self.chro}:{self.pos}"

    @property
    def goodrsid(self) -> bool:
        return self.rsid.startswith("rs")

    def _check_badsnp(self) -> None:
        if self.chro not in valid_chromo:
            self.badsnp = True
        try:
            int(self.pos)
        except ValueError:
            self.badsnp = True
        genotypes = list(self.genotype)
        if not (1 <= len(genotypes) <= 2):
            self.badsnp = True
        if not all(g in valid_call for g in genotypes):
            self.badsnp = True

    def _parse_23andme(self) -> None:
        try:
            self.rsid, self.chro, self.pos, self.genotype = self.line.rstrip().split()
        except ValueError:
            self.badsnp = True

    def _parse_ancestrydna(self) -> None:
        try:
            self.rsid, chro, self.pos, a1, a2 = self.line.rstrip().split()
        except ValueError:
            try:
                self.rsid, chro, self.pos, a1 = self.line.rstrip().split()
                a2 = "_"
            except ValueError:
                self.badsnp = True
                return
        self.chro = {"23": "X", "25": "X", "24": "Y", "26": "MT"}.get(chro, chro)
        self.genotype = a1 + a2

    def _parse_ftdna(self) -> None:
        try:
            self.rsid, self.chro, self.pos, self.genotype = self.line.replace('"', "").split(",")
        except ValueError:
            try:
                self.rsid, self.chro, self.pos, self.genotype = self.line.split()
            except ValueError:
                self.badsnp = True
                return
        if not self.rsid.startswith("rs"):
            self.badsnp = True

    def _parse_myheritage(self) -> None:
        try:
            self.rsid, self.chro, self.pos, self.genotype = self.line.replace('"', "").split(",")
        except ValueError:
            try:
                self.rsid, self.chro, self.pos, self.genotype = self.line.split()
            except ValueError:
                self.badsnp = True
                return
        if not self.rsid.startswith("rs"):
            self.badsnp = True

    def _parse_decodeme(self) -> None:
        try:
            self.rsid, _var, self.chro, self.pos, strand, genotype = self.line.replace(
                " ", ""
            ).split(",")
        except ValueError:
            self.badsnp = True
            return
        self.genotype = flip_genotype(genotype) if strand == "-" else genotype
        if self.chro == "M":
            self.chro = "MT"

    def _parse_livingdnacsv(self) -> None:
        try:
            self.rsid, self.chro, self.pos, self.genotype = self.line.split()
        except ValueError:
            self.badsnp = True
        if not self.rsid.startswith("rs"):
            self.badsnp = True

    def _parse_all(self) -> None:
        a1: str | None = None
        a2: str | None = None
        try:
            self.rsid, self.chro, self.pos, self.genotype = self.line.rstrip().split()
        except ValueError:
            try:
                self.rsid, self.chro, self.pos, a1, a2 = self.line.rstrip().split()
            except ValueError:
                try:
                    self.rsid, self.chro, self.pos, self.genotype = self.line.replace(
                        '"', ""
                    ).split(",")
                except ValueError:
                    try:
                        self.rsid, _var, self.chro, self.pos, strand, self.genotype = (
                            self.line.replace(" ", "").split(",")
                        )
                        if strand == "-":
                            self.genotype = flip_genotype(self.genotype)
                    except ValueError:
                        self.badsnp = True
                        return

        self.chro = {"23": "X", "25": "X", "24": "Y", "26": "MT"}.get(self.chro, self.chro)
        if self.chro == "M":
            self.chro = "MT"
        if a1 is not None and a2 is not None:
            self.genotype = a1 + a2
