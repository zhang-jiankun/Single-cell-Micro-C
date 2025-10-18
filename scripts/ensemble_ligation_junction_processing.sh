#!/bin/bash

chromsize=$1

# merge different single cells

cat *_clean_ligation_junction.bg  > clean_ligation_junction.bed
sort -k1,1 -k2,2n clean_ligation_junction.bed > clean_ligation_junction_sorted.bed
bedtools genomecov -i clean_ligation_junction_sorted.bed -g ${chromsize} -bg > clean_ligation_junction_sorted.bg
bedGraphToBigWig clean_ligation_junction_sorted.bg ${chromsize} clean_ligation_junction_sorted.bw 

# plot

function plot() { 
    
    bigwig=$1

    computeMatrix reference-point -R ${ctcf_motifs} -S ${bigwig} -o ${bigwig/bw/CTCF.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${rela_motifs} -S ${bigwig}  -o ${bigwig/bw/RELA.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${mef2c_motifs} -S ${bigwig}  -o ${bigwig/bw/MEF2C.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${irf4_motifs} -S ${bigwig} -o ${bigwig/bw/IRF4.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${spi1_motifs} -S ${bigwig}  -o ${bigwig/bw/SPI1.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${myc_motifs} -S ${bigwig}  -o ${bigwig/bw/MYC.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${usf1_motifs} -S ${bigwig}  -o ${bigwig/bw/USF1.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${yy1_motifs} -S ${bigwig}  -o ${bigwig/bw/YY1.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60
    computeMatrix reference-point -R ${rest_motifs} -S ${bigwig}  -o ${bigwig/bw/REST.matrix.gz} --referencePoint center -b 200 -a 200 -bs 1 -p 60

    plotProfile -m ${bigwig/bw/CTCF.matrix.gz} -out ${bigwig/bw/CTCF.matrix.pdf} --dpi 150 --averageType sum -z CTCF  --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/RELA.matrix.gz} -out ${bigwig/bw/RELA.matrix.pdf} --dpi 150 --averageType sum -z RELA --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/MEF2C.matrix.gz} -out ${bigwig/bw/MEF2C.matrix.pdf} --dpi 150 --averageType sum -z MEF2C  --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/IRF4.matrix.gz} -out ${bigwig/bw/IRF4.matrix.pdf} --dpi 150 --averageType sum -z IRF4 --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
     plotProfile -m ${bigwig/bw/SPI1.matrix.gz} -out ${bigwig/bw/SPI1.matrix.pdf} --dpi 150 --averageType sum -z SPI1  --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/MYC.matrix.gz} -out ${bigwig/bw/MYC.matrix.pdf} --dpi 150 --averageType sum -z MYC --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/USF1.matrix.gz} -out ${bigwig/bw/USF1.matrix.pdf} --dpi 150 --averageType sum -z USF1  --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
    plotProfile -m ${bigwig/bw/YY1.matrix.gz} -out ${bigwig/bw/YY1.matrix.pdf} --dpi 150 --averageType sum -z YY1 --samplesLabel ${prefix} --yAxisLabel "Number of ligation junctions"
}

plot clean_ligation_junction_sorted.bw
