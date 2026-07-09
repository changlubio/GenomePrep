from __future__ import annotations

import os
from datetime import datetime

from genomeprep._types import ProcessResult, SNPRecord
from genomeprep.convert import directconvertSNPtoVCFformat, output_as_23andme
from genomeprep.data import get_data_dir, load_reference_data
from genomeprep.io import resolve_file_for_genome
from genomeprep.parser import DTCparser, DTCrecord
from genomeprep.qc import dtc_read_records, filter_by_badalleles

__version__ = "1.2.0"
__all__ = ["process", "DTCparser", "DTCrecord", "ProcessResult", "__version__"]

_REF_GENOMES = {"GRCh37": "ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz"}


def process(
    infile: str,
    datadir: str | None = None,
    outdir: str = "./outputs",
    outindex: str = "outindex",
) -> ProcessResult:
    """Parse, QC, and convert a DTC genome file. Outputs .23andme.checked and .vcf files."""
    os.makedirs(outdir, exist_ok=True)
    output_23andme = os.path.join(outdir, f"{outindex}.23andme.checked")
    output_vcf = os.path.join(outdir, f"{outindex}.vcf")

    resolved, handle, iszip = resolve_file_for_genome(infile)
    if not resolved:
        msg = "Cannot identify genetic data in your uploaded ZIP file" if iszip else "Cannot identify genetic data in your upload"
        return ProcessResult(status="Invalid", info=msg)
    if handle is None:
        return ProcessResult(status="error", info="Contact us — data found but handle unavailable")

    parser = DTCparser(handle, filename=infile)
    if parser.reject:
        if "grch" in parser.reject_type:
            return ProcessResult(status="Invalid", info="Only GRCh37 assembly is supported.")
        return ProcessResult(status="Invalid", info="Error reading file headers.")

    resolved_datadir = get_data_dir(datadir)
    ref = load_reference_data(resolved_datadir)

    (reject, reject_type, store_snps, count_snps, duplicate_snps,
     count_not_found_in_ref, valid_snps, in_ref, revd_in_ref,
     in_rev_strand, rev_in_ref, rev_revd_in_ref, not_in_list,
     cluster_id, cluster_mcov) = dtc_read_records(
        parser, ref.api_23andme, ref.faidx, ref.the_list, ref.reverse_snps
    )

    excl_duplicates = count_snps - len(duplicate_snps)

    if reject:
        base_info: dict[str, int | float | str] = {
            "badalleles": parser.total_badsnp,
            "count_snps": count_snps,
            "excl_duplicates": excl_duplicates,
            "count_not_found_in_ref": count_not_found_in_ref,
            "valid_snps": valid_snps,
            "in_ref": in_ref,
            "cluster_id": cluster_id if isinstance(cluster_id, str) else str(cluster_id),
            "cluster_mcov": cluster_mcov if isinstance(cluster_mcov, float) else float(cluster_mcov),
        }
        if reject_type == "vcf":
            return ProcessResult(status="vcf", info="VCF files are not processed by this server.")
        if reject_type == "invalid" and count_snps == 0:
            return ProcessResult(status="Invalid", info="Does not look like a genetic data file.")
        return ProcessResult(status=reject_type, info=base_info)

    # cluster_id can be str | int; bad_allele_clusters is keyed by str
    cluster_key = cluster_id if isinstance(cluster_id, str) else str(cluster_id)
    badalleles: set[str] = ref.bad_allele_clusters.get(cluster_key, set())
    store_snps = filter_by_badalleles(store_snps, badalleles)

    genome_info = (
        "##fileformat=VCFv4.2\n"
        f"##fileDate={datetime.now().strftime('%y-%m-%d')}\n"
        f"##source=GenomePrepV1.2(dtc_format={parser.company},clusterversion={cluster_id},"
        f"clustermaxcov={cluster_mcov},from_valid_snps={valid_snps})\n"
        f"##reference={_REF_GENOMES['GRCh37']}\n"
        "##phasing=unphased\n"
        "##INFO=<ID=INDEL,Number=1,Type=String,Description=\"Insert or Deletion\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )

    afterbadallele = output_as_23andme(store_snps, genome_info, output_23andme)
    final_snps, indels, refref, multialleic = directconvertSNPtoVCFformat(
        store_snps, genome_info, outindex, output_vcf
    )

    return ProcessResult(
        status="Processed",
        info={
            "badalleles": parser.total_badsnp,
            "assembly": parser.build,
            "count_snps": count_snps,
            "excl_duplicates": excl_duplicates,
            "count_not_found_in_ref": count_not_found_in_ref,
            "valid_snps": valid_snps,
            "in_ref": in_ref,
            "cluster_id": cluster_id if isinstance(cluster_id, str) else str(cluster_id),
            "cluster_mcov": cluster_mcov if isinstance(cluster_mcov, float) else float(cluster_mcov),
            "afterbadallele": afterbadallele,
            "indels": indels,
            "refref": refref,
            "multialleic": multialleic,
            "final_snps": final_snps,
        },
    )
