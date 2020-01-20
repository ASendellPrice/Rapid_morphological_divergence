# 20 Jan 2020
# R script used to reorder vcf file based on satsuma output

setwd("~/Dropbox/DPhil/IncipientDivergenceFP/ReorderingVCF/")

# ------------------------------------------------------------
# STEP 1: Read in and massage the satsuma and scaffold data
# ------------------------------------------------------------

library(tidyverse)
scaff_ord <- read_csv("ZLat_scaffold_order_from_ZFinch.csv") %>%
  rename(scaffold = sca)
scaff_lengths <- read_tsv("Zost_Scaffold_Lengths.txt")

missing_lengths <- anti_join(scaff_ord, scaff_lengths)

missing_ords <- anti_join(scaff_lengths, scaff_ord)

scaffs <- scaff_ord %>%
  left_join(scaff_lengths) %>%
  rename(ZFCHROM = chr, CHROM = scaffold)

scaffs


# ------------------------------------------------------------
# STEP 2: Read in zosterops lateralis mapped vcf file
# ------------------------------------------------------------

vcf <- read_tsv("../VCFs/FP_inc_Rurutu/ZL_FP_inc_Rurutu_NZ_75ind_64K.recode.vcf", 
                comment = "##", progress = FALSE) %>%
  rename(CHROM = `#CHROM`)

# ------------------------------------------------------------
# STEP 3: Join VCF and scaffs dataframes
# ------------------------------------------------------------

combo <- scaffs %>%
  left_join(vcf)

# ------------------------------------------------------------
# STEP 4: Get Zebra finch positions
# ------------------------------------------------------------

zf_ified <- combo %>%
  mutate(ZFPOS = {
    ML = floor(mean.loc)  # some temp variables to make it easier to express
    Lo2 = floor(length/2)
    L = length
    ifelse(sca.ori == 1,
           ML - Lo2 + POS,           # forward orientation
           ML - Lo2 + (L - POS))     # reverse orientation
  }) %>%
  select(ZFCHROM, ZFPOS, everything()) %>%
  mutate(ZFCHROM = factor(ZFCHROM, levels = unique(ZFCHROM))) %>%   # this is to get them to sort correctly
  arrange(ZFCHROM, ZFPOS)

head(zf_ified)

# ------------------------------------------------------------
# STEP 4: Remove un-needed columns
# ------------------------------------------------------------

zf_ified$CHROM <- NULL
zf_ified$POS <- NULL
zf_ified$mean.loc <- NULL
zf_ified$sca.ori <- NULL
zf_ified$length <- NULL

# ------------------------------------------------------------
# STEP 5: Save file
# ------------------------------------------------------------
write.table(zf_ified, "../VCFs/ZF_ordered/FP_inc_Rurutu.ZF.ordered.vcf", quote = FALSE, sep = "\t", row.names = FALSE)


# ------------------------------------------------------------
# STEP 6: Add header to make a workable VCF file
# ------------------------------------------------------------
#To make into a functional vcf file need to add header. This is done by copying
# the header from the original vcf file and adding to the new file
# This was done in xcode

# Also need to remove scaffolds not mapped to ZFinch Chroms (NAs)
# sed '/NA/d' ZFified_ZOLA_NonVariant_plus_Biallelic_SNPs_indv163_No_Un.vcf

# and remove lines mapped to negative positions (-)
# sed '/-/d' ZFified_ZOLA_NonVariant_plus_Biallelic_SNPs_indv163_No_Un.vcf