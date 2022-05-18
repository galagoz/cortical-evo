## CONFIGURATION FILE FOR VARIANT-LEVEL QC ##

#-------------------------------------------#
# Define input data and directories         #
#-------------------------------------------#

# Name of the data subset (prefix of SNPstats files), provided through the command line.
# The SNPstats files are required to be in a format <subset_name>_chr<chr>.snpstats.txt
subset_name="ukb43760_enigmaEvo_replication"

# Directory containing the SNP statistics files
stats_dir= "/data/clusterfs/lag/users/barmol/enigma_evo/output/subset_replication"

# Directory where output is written to
working_dir= "/data/clusterfs/lag/users/barmol/enigma_evo/output/vqc_replication/"

# Do you want the script to generate only compact output (only variant identifiers of SNPs to exclude and to keep)? Then compact_output=T.
# Otherwise, when compact_output=F, the script generates a full merged SNPstats table and more detailed lists of SNPs to keep and exclude.
compact_output=F

# Do you want the script to generate plots of MAF and INFO scores, before and after exclusion of filtered SNPs? Then plots=T.
# Otherwise, when plots=F, no plots are generated.
plots=T

#-------------------------------------------#
# Define filtering criteria and thresholds  #
#-------------------------------------------#

# If you do not want to filter on one or more of these criteria, set it as NULL

## MINOR ALLELE FREQUENCY
maf_thr = 0.001 # MAF threshold (based on the MAF in the subset)
maf_diff = 0.2 # MAF difference (based on the absolute difference in MAF between the subset and the total UKB sample, and the subset and the HRC reference panel)

## IMPUTATION QUALITY/INFO SCORE
info_thr = 0.7 # INFO score threshold (based on the INFO score in the subset)
info_diff = NULL  # INFO score difference (based on the absolute difference in INFO score between the subset and the total UKB sample)

## HARDY-WEINBERG EQUILIBRIUM
hwe_thr =1e-6  # Based on the HWE p-value in the subset

## GENOTYPE MISSINGNESS RATE
geno_miss = NULL

## ADDITIONAL FILTERS
remove_indel = F # Logical whether to mark variants other than SNPs (e.g. indertions, deletions) for removal
remove_multiallelic = T # Logical whether to remove multiallelic variants