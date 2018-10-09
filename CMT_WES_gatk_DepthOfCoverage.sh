#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# Basic argument should be assigned with your directory and files
REF=/Directory/Of/Ensembl/CanFam3.1/Sequence/genome.fa
AlignPath
SAMPLE
INTERVAL=$5
INPUT_BAM
SAMPLE=${INPUT_BAM%%.RGadded*}

if [ ! -d $DIR/04_CNVs/DepthOfCoverage ]
then
	mkdir -p $DIR/04_CNVs/DepthOfCoverage
fi
AnalysisPath=$DIR/04_CNVs/DepthOfCoverage/

gatk3=/Directory/Of/GATK_VERSION3/GenomeAnalysisTK.jar

# 1. Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
java -Xmx16g -Djava.io.tmpdir=$DIR/tmp -jar $gatk3 \
-T DepthOfCoverage \
-R $REF \
-I $AlignPath$SAMPLE"/"$INPUT \
-L $INTERVAL \
-dt BY_SAMPLE -dcov 5000 \
-l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 \
--minMappingQuality 20 \
--start 1 \
--stop 5000 \
--nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o $AnalysisPath$SAMPLE.DATA
