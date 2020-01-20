# Custom scripts

This repository contains custom R scripts used in the following manuscript:
Sendell-Price, A.T., Ruegg, K. C. & Clegg, S.M. (in review) Rapid morphological divergence following a human-mediated introduction: The role of drift and selection.

1) process_satsuma_output.R - R script used to process satsuma synteny output. This script assigns zosterops lateralis scaffolds to zebra finch scaffolds and determines mean postion and orientation. Requires the following files: satsuma_output_AssembledChroms.txt (satsuma synteny output - available here: https://www.dropbox.com/s/yocojw8q4b1zpt9/satsuma_output_AssembledChroms.txt?dl=0).

2) Reorder_VCF.R - R Script used to reorder VCF file based on output from satsuma synteny. Requires the following files: ZL_FP_inc_Rurutu_NZ_75ind_64K.recode.vcf (Zosterops lateralis mapped VCF file), ZLat_scaffold_order_from_ZFinch.csv (tab delimited file specifying lateralis scaffold order and orientation relative to zebra finch genome), Zost_Scaffold_Lengths.txt (tab delimited file specifying lateralis scaffold lengths).
