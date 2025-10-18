#!/bin/bash

# Notes: extract and analyse ligation junctions for each single cell

prefix=$1
chromsize=$2
ctcf_motifs=$3 

python ../lib/extract_ligation_junction.py --readid ${prefix}.sam.gz > ${prefix}_ligation_junction.bed
awk -v OFS="\t" '{print $1, $2, $3, ".", $5, $6}'  ${prefix}_ligation_junction.bed | sort -k1,1 -k2,2n | uniq > ${prefix}_clean_ligation_junction.bed
bedtools genomecov -i ${prefix}_clean_ligation_junction.bed -g ${chromsize} -bg > ${prefix}_clean_ligation_junction.bg 
bedGraphToBigWig ${prefix}_clean_ligation_junction.bg ${chromsize} ${prefix}_clean_ligation_junction.bw

# plot genome-average profile

computeMatrix  reference-point -R ${ctcf_motifs} \
    -S ${prefix}_clean_ligation_junction.bw  \
    -o ${prefix}_clean_ligation_junction_pileup.CTCF.matrix.gz \
    --referencePoint center -b 200 -a 200 -bs 1 -p 2

plotProfile -m ${prefix}_clean_ligation_junction_pileup.CTCF.matrix.gz \
    -out ${prefix}_clean_ligation_junction_pileup.CTCF.matrix.pdf --dpi 150 \
    --averageType sum -z CTCF  --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"

