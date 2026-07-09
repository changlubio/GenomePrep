from __future__ import annotations

from genomeprep._types import SNPRecord

_FLIP: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}
_VALID_CALL: frozenset[str] = frozenset("ACGT")


def flip_genotype(genotype: str) -> str:
    return "".join(_FLIP.get(g, g) for g in genotype)


def output_as_23andme(store_snps: list[SNPRecord], genome_info: str, output_file: str) -> int:
    count = 0
    with open(output_file, "w") as fh:
        print(genome_info, file=fh)
        for snp in store_snps:
            count += 1
            print("\t".join([snp["rsid"], snp["chro"], snp["pos"], snp["call"]]), file=fh)
    return count


def directconvertSNPtoVCFformat(
    store_snps: list[SNPRecord],
    vcf_info: str,
    vcf_name: str,
    out_vcf: str,
) -> tuple[int, int, int, int]:
    """Returns (final_snps, indels, refref, multiallelic)."""
    final_snps = indels = refref = multiallelic = 0

    with open(out_vcf, "w") as fh:
        print(vcf_info, file=fh)
        fields_header = [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            vcf_name,
        ]
        print("#" + "\t".join(fields_header), file=fh)

        for snp in store_snps:
            genotypes = list(snp["call"])
            if not (1 <= len(genotypes) <= 2):
                continue

            ref = snp["ref"]
            gset = set(genotypes)

            # INFO: GAV
            info_parts: list[str] = []
            info_parts.append("GAV=none" if snp["gav"] == "" else f"GAV={snp['gav']}")

            # FILTER: badallele
            filter_parts: list[str] = []
            if snp.get("ba", False):
                filter_parts.append("bz3")

            if not all(g in _VALID_CALL for g in gset):
                # INDEL: emit with INDEL tag
                indels += 1
                alt = "."
                info_parts.append(f"INDEL={snp['call']}")
                filter_parts.append("indel")
            else:
                if len(gset) == 1:
                    if ref in gset:
                        # refref: emit with alt="."
                        refref += 1
                        alt = "."
                    else:
                        alt = genotypes[0]
                else:
                    if ref in gset:
                        alt = (gset - {ref}).pop()
                    else:
                        # multiallelic: skip
                        multiallelic += 1
                        continue

            calls = sorted("0" if g == ref else "1" for g in genotypes)
            filter_str = ";".join(filter_parts) if filter_parts else "PASS"
            info_str = ";".join(info_parts)
            row = [
                snp["chro"],
                snp["pos"],
                snp["rsid"],
                ref,
                alt,
                ".",
                filter_str,
                info_str,
                "GT",
                "/".join(calls),
            ]
            print("\t".join(row), file=fh)
            final_snps += 1

    return final_snps, indels, refref, multiallelic
