"""Python port of the TargetScan target-prediction pipeline.

Ported from the original Perl scripts in ``tsh_orig/`` (TargetScan Release
7 code, which already contains the TargetScan 8 context++ scoring model).
The algorithms are genome-agnostic: they operate on aligned UTR/ORF
sequences and miRNA seed tables, so the same code works for hg19 (vert70)
or hg38 (vert80) input data -- only the downloaded data files differ.
"""
