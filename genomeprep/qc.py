from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING, Any

from genomeprep._types import SNPRecord
from genomeprep.convert import flip_genotype

if TYPE_CHECKING:
    import pyfaidx

    from genomeprep.parser import DTCparser

MIN_CLUSTER_COVERAGE: float = 0.8
MIN_COUNT_SNP: int = 10000
_USABLE: frozenset[str] = frozenset("ACGTDI")


def get_snp_reference(loci: str, faidx: pyfaidx.Fasta) -> str | None:
    try:
        chro, pos_str = loci.split(":")
    except ValueError:
        return None
    try:
        return str(faidx[chro][int(pos_str) - 1].seq)
    except Exception:
        return None


def filter_by_badalleles(store_snps: list[SNPRecord], badalleles: set[str]) -> list[SNPRecord]:
    for snp in store_snps:
        snp["ba"] = f"{snp['chro']}:{snp['pos']}" in badalleles
    return store_snps


def dtc_read_records(
    dtc_parser: DTCparser,
    api_23andme: set[str],
    faidx: pyfaidx.Fasta,
    the_list: dict[str, list[str]],
    reverse_snps: dict[str, Any],
) -> tuple[
    bool,
    str,
    list[SNPRecord],
    int,
    list[str],
    int,
    int,
    float,
    float,
    int,
    float,
    float,
    float,
    str | int,
    float | int,
]:
    """
    Read all records from a DTCparser and perform quality checks.

    Returns:
        reject, reject_type, store_snps, count_snps, duplicate_snps,
        count_not_found_in_ref, valid_snps, in_ref, revd_in_ref,
        in_rev_strand, rev_in_ref, rev_revd_in_ref, not_in_list,
        cluster_id, cluster_mcov
    """
    reject = False
    reject_type = ""

    count_snps = 0
    count_not_found_in_ref = 0

    vcf_lines = 0
    non_vcf_lines = 0

    # valid_snps: exclude no calls;
    # also make sure that a reference genotype can be found,
    # otherwise it may be from a different build
    valid_snps = 0
    store_snps: list[SNPRecord] = []
    snps_seen: set[str] = set()
    duplicate_snps: list[str] = []
    count_genotype_contains_reference = 0

    count_not_in_list = 0
    cluster_count: Counter[str] = Counter()

    try:
        for record in dtc_parser:
            if record.badsnp:
                dtc_parser.total_badsnp += 1
                # check vcf
                if len(record.line.split("\t")) > 9:
                    vcf_lines += 1
                else:
                    non_vcf_lines += 1
                continue

            # special treatment for 23andme: filter to known loci when api set is provided
            if dtc_parser.company == "23andme" and record.loci not in api_23andme:
                record.badsnp = True
                dtc_parser.total_badsnp += 1
                continue

            count_snps += 1

            # resolve duplicates
            snpid = f"{record.chro}_{record.pos}"
            if snpid in snps_seen:
                duplicate_snps.append(f"{record.chro}_{record.pos}_{record.genotype}")
                continue
            snps_seen.add(snpid)

            if _USABLE.intersection(record.genotype):
                ref = get_snp_reference(record.loci, faidx)

                if ref and len(ref) == 1:
                    valid_snps += 1
                    # compare to reference
                    if ref in record.genotype:
                        count_genotype_contains_reference += 1

                    # cluster check
                    gavs: list[str] = []
                    if record.loci in the_list:
                        for clust in the_list[record.loci]:
                            cluster_count[clust] += 1
                            gavs.append(clust)
                    else:
                        count_not_in_list += 1

                    store_snps.append(
                        {
                            "rsid": record.rsid,
                            "chro": record.chro,
                            "pos": record.pos,
                            "call": record.genotype,
                            "ref": ref,
                            "gav": ",".join(gavs),
                        }
                    )
                else:
                    count_not_found_in_ref += 1

    except RecursionError:
        reject = True
        reject_type = "invalid"

    # invalid if parsed lines < 10000
    if not reject:
        if count_snps < MIN_COUNT_SNP:
            reject = True
            reject_type = "invalid"

    # reject if <22 chromosomes or fewer than 10000 valid SNPs
    if not reject:
        chromo_set = {s["chro"] for s in store_snps}
        if len(chromo_set) < 22 or valid_snps < 10000:
            reject = True
            reject_type = "incomplete"

    # assembly check
    if not reject:
        in_ref = round(count_genotype_contains_reference / valid_snps, 4)
        if in_ref < 0.6 or count_not_found_in_ref > 0:
            reject = True
            reject_type = "notgrch37"
    else:
        in_ref = -1.0

    # imputation check
    if not reject:
        if count_snps > 2_000_000:
            reject = True
            reject_type = "imputation"

    # reverse strand check
    in_rev_strand: int = -1
    rev_in_ref: float = -1.0
    rev_revd_in_ref: float = -1.0
    revd_in_ref: float = -1.0
    if not reject:
        count_rev_in_ref = 0
        count_rev_revd_in_ref = 0
        count_revd_in_ref = 0
        in_rev_strand = 0

        for snp in store_snps:
            newcall = flip_genotype(snp["call"])
            if snp["ref"] in newcall:
                count_revd_in_ref += 1

            snppos = f"{snp['chro']}:{snp['pos']}"
            if snppos in reverse_snps:
                in_rev_strand += 1
                if snp["ref"] in snp["call"]:
                    count_rev_in_ref += 1
                if snp["ref"] in newcall:
                    count_rev_revd_in_ref += 1

        rev_revd_in_ref = round(count_rev_revd_in_ref / (in_rev_strand + 1), 4)
        rev_in_ref = round(count_rev_in_ref / (in_rev_strand + 1), 4)
        revd_in_ref = round(count_revd_in_ref / valid_snps, 4)

        if rev_revd_in_ref > 0.1 or in_ref < 0.75 or revd_in_ref > 0.08 or in_ref > 0.98:
            reject = True
            reject_type = "ref"

    # cluster info
    not_in_list: float = -1.0
    cluster_id: str | int = -1
    cluster_mcov: float | int = -1
    if not reject:
        not_in_list = round(count_not_in_list / valid_snps, 3)
        try:
            cluster_id, count = cluster_count.most_common(1)[0]
            cluster_mcov = round(count / valid_snps, 2)
            if cluster_mcov < MIN_CLUSTER_COVERAGE:
                cluster_id = "all"
        except IndexError:
            cluster_id = cluster_mcov = -1

    if reject:
        total = vcf_lines + non_vcf_lines
        if vcf_lines > 5000 and total > 0 and non_vcf_lines / total < 0.1:
            reject_type = "vcf"

    return (
        reject,
        reject_type,
        store_snps,
        count_snps,
        duplicate_snps,
        count_not_found_in_ref,
        valid_snps,
        in_ref,
        revd_in_ref,
        in_rev_strand,
        rev_in_ref,
        rev_revd_in_ref,
        not_in_list,
        cluster_id,
        cluster_mcov,
    )
