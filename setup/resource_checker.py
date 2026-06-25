"""
resource_checker.py — Validate reference genome and annotation resource paths.

Checks existence, indexing, and basic format expectations.
Does not download anything — only validates what is already present.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Resources marked required=True will cause setup to fail if missing.
# Resources marked required=False will be reported as warnings.
RESOURCE_SPEC: Dict[str, dict] = {
    "reference_fasta": {
        "required": True,
        "extensions": [".fa", ".fasta", ".fna"],
        "companion_suffixes": [".fai"],        # .dict is also expected but not enforced here
        "description": "GRCh38 reference FASTA (with .fai index)",
    },
    "vep_cache_dir": {
        "required": True,
        "is_dir": True,
        "description": "VEP offline cache directory (GRCh38)",
    },
    "vep_plugins_dir": {
        "required": False,
        "is_dir": True,
        "description": "VEP plugins directory",
        "user_hint": (
            "Directory containing VEP plugin .pm files (e.g. ~/.vep/Plugins). "
            "Download plugins from: github.com/Ensembl/VEP_plugins"
        ),
    },
    "clinvar_vcf": {
        "required": False,
        "extensions": [".vcf.gz"],
        "companion_suffixes": [".tbi"],
        "description": "ClinVar VCF (GRCh38, bgzipped + tabix indexed)",
        "user_hint": (
            "File: clinvar.vcf.gz  (with clinvar.vcf.gz.tbi)\n"
            "  Download: ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/\n"
            "  Prepare:  bgzip clinvar.vcf && tabix -p vcf clinvar.vcf.gz"
        ),
    },
    "gnomad_vcf": {
        "required": False,
        "extensions": [".vcf.gz"],
        "companion_suffixes": [".tbi"],
        "description": "gnomAD sites VCF (GRCh38, bgzipped + tabix indexed)",
        "user_hint": (
            "File: gnomad.genomes.vX.vcf.bgz  (with .tbi)\n"
            "  Download: gnomad.broadinstitute.org/downloads\n"
            "  Note: gnomAD AF is already available from the VEP cache (--af_gnomad);\n"
            "  this file is only needed for additional fields not in the cache."
        ),
    },
    "cadd_snv_file": {
        "required": False,
        "extensions": [".tsv.gz"],
        "companion_suffixes": [".tbi"],
        "description": "CADD SNV scores (whole_genome_SNVs.tsv.gz + .tbi)",
        "user_hint": (
            "File: whole_genome_SNVs.tsv.gz  (with .tbi)\n"
            "  Download: cadd.gs.washington.edu/download  (GRCh38-v1.7 or later)"
        ),
    },
    "cadd_indel_file": {
        "required": False,
        "extensions": [".tsv.gz"],
        "companion_suffixes": [".tbi"],
        "description": "CADD indel scores (.tsv.gz + .tbi)",
        "user_hint": (
            "File: gnomad.genomes.rX.indel.tsv.gz  (with .tbi)\n"
            "  Download: cadd.gs.washington.edu/download  (GRCh38 indel file)"
        ),
    },
    "revel_file": {
        "required": False,
        "extensions": [".tsv.gz", ".csv.gz", ".txt.gz"],
        "description": "REVEL missense pathogenicity scores for VEP plugin",
        "user_hint": (
            "File: revel_all_chromosomes.tsv.gz  (bgzipped + tabix indexed)\n"
            "  Download: sites.google.com/site/revel-tool  (Whole genome REVEL scores)\n"
            "  Prepare:  bgzip revel_all_chromosomes.tsv\n"
            "            tabix -f -s 1 -b 2 -e 2 revel_all_chromosomes.tsv.gz"
        ),
    },
    "spliceai_snv_file": {
        "required": False,
        "extensions": [".vcf.gz"],
        "companion_suffixes": [".tbi"],
        "description": "SpliceAI SNV scores (hg38, bgzipped + tabix indexed)",
        "user_hint": (
            "File: spliceai_scores.raw.snv.hg38.vcf.gz  (with .tbi)\n"
            "  Download: Illumina BaseSpace (search 'SpliceAI lookup') or Zenodo\n"
            "  Requires spliceai_indel_file as well (both SNV and indel are needed)"
        ),
    },
    "spliceai_indel_file": {
        "required": False,
        "extensions": [".vcf.gz"],
        "companion_suffixes": [".tbi"],
        "description": "SpliceAI indel scores (hg38, bgzipped + tabix indexed)",
        "user_hint": (
            "File: spliceai_scores.raw.indel.hg38.vcf.gz  (with .tbi)\n"
            "  Download: Illumina BaseSpace (search 'SpliceAI lookup') or Zenodo\n"
            "  Requires spliceai_snv_file as well (both SNV and indel are needed)"
        ),
    },
    "alphamissense_file": {
        "required": False,
        "extensions": [".tsv.gz", ".tsv"],
        "description": "AlphaMissense pathogenicity scores for VEP plugin",
        "user_hint": (
            "File: AlphaMissense_hg38.tsv.gz\n"
            "  Download: Zenodo (search 'AlphaMissense hg38') or\n"
            "            storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz\n"
            "  Prepare:  bgzip AlphaMissense_hg38.tsv && tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz"
        ),
    },
    "loftee_path": {
        "required": False,
        "is_dir": True,
        "description": "LOFTEE GRCh38 source directory (contains LoF.pm, gerp_dist.pl, splice_data/, maxEntScan/)",
        "user_hint": (
            "Clone: git clone --branch grch38 https://github.com/konradjk/loftee.git ~/.vep/loftee_grch38\n"
            "  This is the full LOFTEE source, separate from vep_plugins_dir."
        ),
    },
    "loftee_human_ancestor_fa": {
        "required": False,
        "extensions": [".fa.gz"],
        "companion_suffixes": [".fai", ".gzi"],
        "description": "LOFTEE human ancestor FASTA for GRCh38 (bgzipped + indexed)",
        "user_hint": (
            "File: human_ancestor_GRCh38.fa.gz  (with .fai and .gzi)\n"
            "  Download: github.com/konradjk/loftee — GRCh38 data files\n"
            "  Required for LOFTEE ancestral allele filter"
        ),
    },
    "loftee_conservation_file": {
        "required": False,
        "extensions": [".sql"],
        "description": "LOFTEE GERP conservation SQLite database (loftee.sql)",
        "user_hint": (
            "File: loftee.sql\n"
            "  Download: github.com/konradjk/loftee — GRCh38 data files\n"
            "  GERP conservation scores database used by LOFTEE"
        ),
    },
    "loftee_gerp_bigwig": {
        "required": False,
        "extensions": [".bw"],
        "description": "LOFTEE GERP conservation BigWig (GRCh38)",
        "user_hint": (
            "File: gerp_conservation_scores.homo_sapiens.GRCh38.bw\n"
            "  Download: github.com/konradjk/loftee — GRCh38 data files\n"
            "  Used by LOFTEE for BigWig-based conservation filtering"
        ),
    },
    "clingen_gene_validity_csv": {
        "required": False,
        "extensions": [".csv"],
        "description": "ClinGen Gene Disease Validity Curations CSV",
        "user_hint": (
            "File: ClinGen_gene_validity.csv\n"
            "  Download: search.clinicalgenome.org/kb/gene-validity/download\n"
            "  Provides gene-level disease validity classification "
            "(Definitive/Strong/Moderate/Limited/Disputed/Refuted)"
        ),
    },
    "bwa_index_prefix": {
        "required": False,
        "is_prefix": True,
        "description": (
            "Legacy bwa mem 0.7.17 index prefix — superseded by minibwa_index_prefix. "
            "Kept for reference; pipeline now uses minibwa."
        ),
        "user_hint": (
            "This resource is no longer used by the pipeline (minibwa replaced bwa mem).\n"
            "  Configure minibwa_index_prefix instead."
        ),
    },
    "minibwa_exe": {
        "required": False,
        "is_exe": True,
        "description": "minibwa executable (v0.2+). Required for FASTQ input.",
        "user_hint": (
            "Path to the minibwa binary, e.g. /home/user/tools/minibwa/minibwa\n"
            "  Build: git clone --depth=1 https://github.com/lh3/minibwa.git && make -C minibwa -j8\n"
            "  Or add minibwa to PATH — setup will find it automatically."
        ),
    },
    "minibwa_index_prefix": {
        "required": False,
        "is_prefix": True,
        "description": (
            "minibwa index prefix (path without extension, e.g. /data/minibwa/genome). "
            "Build with: minibwa index -l <reference.fa>"
        ),
        "user_hint": (
            "Prefix path (no extension), e.g. /data/minibwa/genome\n"
            "  The index files genome.mbw and genome.l2b must exist at that prefix.\n"
            "  Build: minibwa index -l /path/to/genome.fa  (takes ~40 min for GRCh38)\n"
            "  The -l flag is required for genomes larger than 2 GB (e.g. hg38).\n"
            "  Only needed for FASTQ input — VCF input does not use minibwa."
        ),
    },
}

# All five files that bwa index produces alongside the reference FASTA (legacy).
_BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

# Minimum acceptable size for the .bwt file.  A real GRCh38 .bwt is ~3 GB.
_BWA_BWT_MIN_BYTES = 1_048_576  # 1 MiB

# minibwa index files: .mbw (main index, ~5.8 GB for GRCh38) and .l2b (lookup table, ~800 MB)
_MINIBWA_INDEX_SUFFIXES = (".mbw", ".l2b")
_MINIBWA_MBW_MIN_BYTES = 1_048_576  # 1 MiB


class ResourceCheckResult:
    """Result for a single resource check."""

    def __init__(self, resource: str):
        self.resource = resource
        self.path: Optional[str] = None
        self.status: str = "not_checked"   # found | missing | bad_format | not_provided
        self.required: bool = RESOURCE_SPEC.get(resource, {}).get("required", False)
        self.messages: List[str] = []

    def as_dict(self) -> dict:
        return {
            "path": self.path,
            "status": self.status,
        }


def _check_exe_resource(resource: str, path: str, spec: dict) -> ResourceCheckResult:
    """Validate an executable resource — must exist and be executable."""
    result = ResourceCheckResult(resource)
    result.path = path
    p = Path(path)
    if not p.exists():
        result.status = "missing"
        result.messages.append(
            f"Executable not found: {path}\n"
            f"  Expected: {spec['description']}"
        )
        return result
    if not p.is_file():
        result.status = "bad_format"
        result.messages.append(f"Path exists but is not a file: {path}")
        return result
    import os as _os
    if not _os.access(path, _os.X_OK):
        result.status = "bad_format"
        result.messages.append(f"File is not executable: {path}\n  Run: chmod +x {path}")
        return result
    result.status = "found"
    result.messages.append(f"OK — {path}")
    return result


def _check_prefix_resource(resource: str, path: str, spec: dict) -> ResourceCheckResult:
    """
    Validate an index prefix resource.

    For bwa_index_prefix: checks .amb/.ann/.bwt/.pac/.sa files exist.
    For minibwa_index_prefix: checks .mbw and .l2b files exist.
    Falls back to bwa validation for any unrecognised prefix resource.
    """
    result = ResourceCheckResult(resource)
    result.path = path

    if resource == "minibwa_index_prefix":
        missing = [s for s in _MINIBWA_INDEX_SUFFIXES if not Path(path + s).exists()]
        if missing:
            result.status = "missing"
            result.messages.append(
                f"minibwa index incomplete at prefix: {path}\n"
                f"  Missing files: {', '.join(path + s for s in missing)}\n"
                f"  Build the index (takes ~40 min for GRCh38; -l required for large genomes):\n"
                f"    minibwa index -l <genome.fa> {path}\n"
                f"  Or point resources.minibwa_index_prefix.path at an existing prefix."
            )
            return result
        mbw_path = Path(path + ".mbw")
        mbw_size = mbw_path.stat().st_size
        if mbw_size <= _MINIBWA_MBW_MIN_BYTES:
            result.status = "bad_format"
            result.messages.append(
                f"minibwa index appears to be a placeholder at prefix: {path}\n"
                f"  {mbw_path.name} is {mbw_size:,} bytes "
                f"(minimum expected: {_MINIBWA_MBW_MIN_BYTES:,} bytes for a real index)\n"
                f"  A real GRCh38 .mbw is approximately 5.8 GB.\n"
                f"  Rebuild: minibwa index -l <genome.fa> {path}"
            )
            return result
        result.status = "found"
        result.messages.append(f"OK — {path}")
        return result

    # bwa_index_prefix (legacy) — check five standard bwa files
    missing = [s for s in _BWA_INDEX_SUFFIXES if not Path(path + s).exists()]
    if missing:
        result.status = "missing"
        result.messages.append(
            f"BWA index incomplete at prefix: {path}\n"
            f"  Missing files: {', '.join(path + s for s in missing)}\n"
            f"  Note: pipeline now uses minibwa — configure minibwa_index_prefix instead."
        )
        return result

    bwt_path = Path(path + ".bwt")
    bwt_size = bwt_path.stat().st_size
    if bwt_size <= _BWA_BWT_MIN_BYTES:
        result.status = "bad_format"
        result.messages.append(
            f"BWA index appears to be a placeholder at prefix: {path}\n"
            f"  {bwt_path.name} is {bwt_size:,} bytes "
            f"(minimum expected: {_BWA_BWT_MIN_BYTES:,} bytes for a real index)"
        )
        return result

    result.status = "found"
    result.messages.append(f"OK (legacy bwa index) — {path}")
    return result


def _check_file_resource(resource: str, path: str, spec: dict) -> ResourceCheckResult:
    result = ResourceCheckResult(resource)
    result.path = path
    p = Path(path)

    if not p.exists():
        result.status = "missing"
        result.messages.append(
            f"File not found: {path}\n"
            f"  Expected: {spec['description']}"
        )
        return result

    if not p.is_file():
        result.status = "bad_format"
        result.messages.append(f"Path exists but is not a file: {path}")
        return result

    # Check companion index files (e.g. .tbi)
    for suffix in spec.get("companion_suffixes", []):
        companion = Path(str(path) + suffix)
        if not companion.exists():
            result.messages.append(
                f"WARNING: companion index not found: {companion}\n"
                f"  Run: tabix -p vcf {path}"
            )

    result.status = "found"
    result.messages.append(f"OK — {path}")
    return result


def _check_dir_resource(resource: str, path: str, spec: dict) -> ResourceCheckResult:
    result = ResourceCheckResult(resource)
    result.path = path
    p = Path(path)

    if not p.exists():
        result.status = "missing"
        result.messages.append(
            f"Directory not found: {path}\n"
            f"  Expected: {spec['description']}"
        )
        return result

    if not p.is_dir():
        result.status = "bad_format"
        result.messages.append(f"Path exists but is not a directory: {path}")
        return result

    result.status = "found"
    result.messages.append(f"OK — {path}")
    return result


def check_resource(resource: str, path: Optional[str]) -> ResourceCheckResult:
    """
    Validate a single resource.

    Parameters
    ----------
    resource:  Resource key (must be in RESOURCE_SPEC).
    path:      Path string supplied by user/config, or None if not provided.
    """
    spec = RESOURCE_SPEC.get(resource, {})
    result = ResourceCheckResult(resource)

    if not path:
        result.status = "not_provided"
        if spec.get("required", False):
            result.messages.append(
                f"REQUIRED resource not provided: {resource}\n"
                f"  Expected: {spec.get('description', resource)}\n"
                f"  Provide the path in config/paths.yaml under resources.{resource}.path"
            )
        else:
            hint = spec.get("user_hint", "")
            hint_block = f"\n  File needed:\n    {hint.replace(chr(10), chr(10) + '    ')}" if hint else ""
            result.messages.append(
                f"Optional resource not provided: {resource}\n"
                f"  {spec.get('description', resource)} — skipped (some annotations will be unavailable)"
                f"{hint_block}\n"
                f"  To configure: set resources.{resource}.path in config/paths.yaml,\n"
                f"                then re-run: python cli.py --setup --force"
            )
        return result

    if spec.get("is_exe", False):
        return _check_exe_resource(resource, path, spec)
    elif spec.get("is_prefix", False):
        return _check_prefix_resource(resource, path, spec)
    elif spec.get("is_dir", False):
        return _check_dir_resource(resource, path, spec)
    else:
        return _check_file_resource(resource, path, spec)


def check_all_resources(provided_paths: Dict[str, Optional[str]]) -> Dict[str, ResourceCheckResult]:
    """
    Check all resources listed in RESOURCE_SPEC.

    Parameters
    ----------
    provided_paths: Mapping of resource_name → path string (or None).

    Returns
    -------
    Dict of resource_name → ResourceCheckResult.
    """
    results = {}
    for resource in RESOURCE_SPEC:
        results[resource] = check_resource(resource, provided_paths.get(resource))
    return results
