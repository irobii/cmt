#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# Basic argument should be assigned with your directory and files
PROJECT
REF
AlignedPath
INPUT=$4
ID=${INPUT%%.RGadded*}
DIR=/data/project/$PROJECT

# fixed annotation data path
# Single nucleotide polymorphism in dbSNP
dbSNP=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris.vcf.gz
# All structural variations
SV=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_structural_variations.vcf.gz
# All consequences of the variations on the Ensembl transcriptome, as called by the variation consequence pipeline
VAR=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_incl_consequences.vcf.gz

# step 1. Variant call by GATK
gatk HaplotypeCaller \
-R $REF \
-I $AlignedPath$ID"/"$INPUT \
--dbsnp $dbSNP \
-stand-call-conf 30 \
--min-pruning 3 \
-O $AnalysisPath$ID.Haplotype.raw.vcf

# step 2. Extract the SNPs and INDELs each from the call set
# step 2-1. Extract the SNPs from the call set
gatk SelectVariants \
-R $REF \
-V $AnalysisPath$ID.Haplotype.raw.vcf \
-select-type SNP \
-O $AnalysisPath$ID.raw_snps.vcf

# step 2-2. Filter to the SNP call set
gatk VariantFiltration \
-R $REF \
-V $AnalysisPath$ID.raw_snps.vcf \
-O $AnalysisPath$ID.filtered_snps.vcf \
-filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filter-name "HARD_TO_VALIDATE" \
-filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "GATK_snp" \
-filter "DP < 5 " \
--filter-name "LowCoverage" \
-filter "QUAL < 30.0 " \
--filter-name "VeryLowQual" \
-filter "QUAL > 30.0 && QUAL < 50.0 " \
--filter-name "LowQual" \
-filter "QD < 1.5 " \
--filter-name "LowQD" \
-filter "SB > -10.0 " \
--filter-name "StrandBias"

bgzip -c $AnalysisPath$ID.filtered_snps.vcf >$AnalysisPath$ID.filtered_snps.vcf.gz
tabix -p vcf $AnalysisPath$ID.filtered_snps.vcf.gz

# step 2-3. Extract the Indels from the call set
gatk SelectVariants \
-R $REF \
-V $AnalysisPath$ID.Haplotype.raw.vcf \
-select-type INDEL \
-O $AnalysisPath$ID.raw_indels.vcf

# step 2-4. Filter to the Indel call set
gatk VariantFiltration \
-R $REF \
-V $AnalysisPath$ID.raw_indels.vcf \
-O $AnalysisPath$ID.filtered_indels.vcf \
-filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filter-name "HARD_TO_VALIDATE" \
-filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "GATK_indel" \
-filter "DP < 5 " \
--filter-name "LowCoverage" \
-filter "QUAL < 30.0 " \
--filter-name "VeryLowQual" \
-filter "QUAL > 30.0 && QUAL < 50.0 " \
--filter-name "LowQual" \
-filter "QD < 1.5 " \
--filter-name "LowQD" \
-filter "SB > -10.0 " \
--filter-name "StrandBias"
