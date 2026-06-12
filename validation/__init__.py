"""
validation/ — Internal QC validation framework for VI_Pipeline.

Run a validation:
    python validate.py --truth-set test_data/truth_set.tsv \\
                       --vcf runs/pt1/pt1.norm.vep.vcf \\
                       --output-dir validation_out/

Produces:
    validation_out/validation_summary.pdf
    validation_out/discordance_log.xlsx
    validation_out/tool_concordance.xlsx
    validation_out/raw_results.tsv
    validation_out/raw_results.json
"""
