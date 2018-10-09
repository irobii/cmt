#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# basic argument
PROJECT=$1
REF=$2
AlignPath=$3
NORMAL=$4
TUMOR=$5
INTERVAL=$6
# ID=${TUMOR%%.bam}
ID=${TUMOR%%-tumor.RGadded.marked.realigned.fixed.recal.bam}
DIR=/data/project/$PROJECT


if [ ! -d $DIR/03_somatic/mutect2 ]
then
	mkdir -p $DIR/03_somatic/mutect2
fi
AnalysisPath=$DIR/03_somatic/mutect2/
# gatk4=/opt/Yonsei/GATK/4.0.4.0/gatk

# fixed annotation data path
dbSNP=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris.vcf.gz
SV=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_structural_variations.vcf.gz # All structural variations
VAR=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_incl_consequences.vcf.gz # All consequences of the variations on the Ensembl transcriptome, as called by the variation consequence pipeline

# print start time
date

# 1. Variant call by muTect2
gatk Mutect2 \
-R $REF \
-I $AlignPath$ID"-tumor/"$TUMOR \
-I $AlignPath$ID"-normal/"$NORMAL \
-L $INTERVAL \
-tumor ${TUMOR%%.RGadded.marked.realigned.fixed.recal.bam} \
-normal ${NORMAL%%.RGadded.marked.realigned.fixed.recal.bam} \
-O $AnalysisPath$ID.mutect2.somatic.vcf

# 2. Filter mutation variant call
gatk FilterMutectCalls \
-V $AnalysisPath$ID.mutect2.somatic.vcf \
-O $AnalysisPath$ID.mutect2.somatic.filter.vcf

bgzip -c $AnalysisPath$ID.mutect2.somatic.filter.vcf >$AnalysisPath$ID.mutect2.somatic.filter.vcf.gz
tabix -p vcf $AnalysisPath$ID.mutect2.somatic.filter.vcf.gz

# 3. Select passed mutation variant call
egrep '#|PASS' $AnalysisPath$ID.mutect2.somatic.filter.vcf >$AnalysisPath$ID.mutect2.somatic.pass.vcf

# 4. Select variant types
gatk SelectVariants -select-type SNP -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.SNP.vcf
gatk SelectVariants -select-type INDEL -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.INDEL.vcf
gatk SelectVariants -select-type MNP -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.MNP.vcf
gatk SelectVariants -select-type MIXED -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.MIXED.vcf
gatk SelectVariants -select-type SYMBOLIC -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.SYMBOLIC.vcf
gatk SelectVariants -select-type NO_VARIATION -V $AnalysisPath$ID.mutect2.somatic.filter.vcf -O $AnalysisPath$ID.mutect2.somatic.NO_VARIATION.vcf

gatk SelectVariants -select-type SNP -V $AnalysisPath$ID.mutect2.somatic.pass.vcf -O $AnalysisPath$ID.mutect2.somatic.SNP.pass.vcf
gatk SelectVariants -select-type INDEL -V $AnalysisPath$ID.mutect2.somatic.pass.vcf -O $AnalysisPath$ID.mutect2.somatic.INDEL.pass.vcf
gatk SelectVariants -select-type MNP -V $AnalysisPath$ID.mutect2.somatic.pass.vcf -O $AnalysisPath$ID.mutect2.somatic.MNP.pass.vcf
gatk SelectVariants -select-type MIXED -V $AnalysisPath$ID.mutect2.somatic.pass.vcf -O $AnalysisPath$ID.mutect2.somatic.MIXED.pass.vcf

# 5. Stastics of variants
vcf-stats $AnalysisPath$ID.mutect2.somatic.SNP.vcf -p $AnalysisPath"stats/"$ID.mutect2.somatic.SNP
# vcf-stats $AnalysisPath$ID.mutect2.somatic.SNP.pass.vcf -p $AnalysisPath"stats/"$ID.mutect2.somatic.SNP.pass

# 6. bgzip & tabix
bgzip -c $AnalysisPath$ID.mutect2.somatic.vcf >$AnalysisPath$ID.mutect2.somatic.vcf.gz
tabix -p vcf $AnalysisPath$ID.mutect2.somatic.vcf.gz

bgzip -c $AnalysisPath$ID.mutect2.somatic.pass.vcf >$AnalysisPath$ID.mutect2.somatic.pass.vcf.gz
tabix -p vcf $AnalysisPath$ID.mutect2.somatic.pass.vcf.gz

# print end time
date
