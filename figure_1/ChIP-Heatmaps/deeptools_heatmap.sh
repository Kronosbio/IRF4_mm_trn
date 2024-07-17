###Steps taken

#1) install (https://deeptools.readthedocs.io/en/develop/content/installation.html):
#python setup.py install (--prefix /User/Tools/deepTools2.0)
#
#2) bam to bigwig (needs to have .BAI index in folder) - if you already have bigWigs no need.
#
#3) computeMatrix (https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html)
#
#4) plotHeatmap (https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html)

###first generating bigwigs, replicates require additional steps.

bamCoverage -b ./bams/01_0EWN_0229Kronos_MM1S-No-Treatment-1_H3K18Ac_hs-dm_i02.bam -o ./bw_reps/KB_MM1S_H3K18ac_1.bw 
bamCoverage -b ./bams/02_0EWO_0229Kronos_MM1S-No-Treatment-2_H3K18Ac_hs-dm_i04.bam -o ./bw_reps/KB_MM1S_H3K18ac_2.bw 
bamCoverage -b ./bams/11_0EWX_0229Kronos_MM1S-No-Treatment-1_H3K27Ac_hs-dm_i20.bam -o ./bw_reps/KB_MM1S_H3K27ac_1.bw  
bamCoverage -b ./bams/12_0EWY_0229Kronos_MM1S-No-Treatment-2_H3K27Ac_hs-dm_i21.bam -o ./bw_reps/KB_MM1S_H3K27ac_2.bw 
bamCoverage -b ./bams/21_0F27_0229Kronos_MM1S-No-Treatment-1_H3K27me3_hs-dm_i02.bam -o ./bw_reps/KB_MM1S_H3K27me3_1.bw 
bamCoverage -b ./bams/22_0F28_0229Kronos_MM1S-No-Treatment-2_H3K27me3_hs-dm_i04.bam -o ./bw_reps/KB_MM1S_H3K27me3_2.bw 
bamCoverage -b ./bams/31_0EX7_0229Kronos_MM1S-No-Treatment-1_IRF4_hs_i45.bam -o ./bw_reps/KB_MM1S_IRF4_1.bw 
bamCoverage -b ./bams/32_0EX8_0229Kronos_MM1S-No-Treatment-2_IRF4_hs_i48.bam -o ./bw_reps/KB_MM1S_IRF4_2.bw 
bamCoverage -b ./bams/33_0EX9_0229Kronos_MM1S-No-Treatment-1_p300_hs_i49.bam -o ./bw_reps/KB_MM1S_p300_1.bw 
bamCoverage -b ./bams/34_0EXA_0229Kronos_MM1S-No-Treatment-2_p300_hs_i50.bam -o ./bw_reps/KB_MM1S_p300_2.bw 

bamCoverage -b ./bams/barwick_MM1S_BRD4.bam -o ./bw/tfs/barwick_MM1S_BRD4.bw 
bamCoverage -b ./bams/barwick_MM1S_IKZF1.bam -o ./bw/tfs/barwick_MM1S_IKZF1.bw 
bamCoverage -b ./bams/loven_MM1S_RNA_POL2.bam -o ./bw/tfs/loven_MM1S_RNA_POL2.bw 

### merging replicates

bigWigMerge ./bw_reps/KB_MM1S_IRF4_1.bw ./bw_reps/KB_MM1S_IRF4_2.bw  ./bedgraph/KB_MM1S_IRF4.bedGraph
bigWigMerge ./bw_reps/KB_MM1S_p300_1.bw ./bw_reps/KB_MM1S_p300_2.bw  ./bedgraph/KB_MM1S_p300.bedGraph
bigWigMerge ./bw_reps/KB_MM1S_H3K27ac_1.bw ./bw_reps/KB_MM1S_H3K27ac_2.bw  ./bedgraph/KB_MM1S_H3K27ac.bedGraph
bigWigMerge ./bw_reps/KB_MM1S_H3K27me3_1.bw ./bw_reps/KB_MM1S_H3K27me3_2.bw  ./bedgraph/KB_MM1S_H3K27me3.bedGraph
bigWigMerge ./bw_reps/KB_MM1S_H3K18ac_1.bw ./bw_reps/KB_MM1S_H3K18ac_2.bw  ./bedgraph/KB_MM1S_H3K18ac.bedGraph

sort -k1,1 -k2,2n  ./bedgraph/KB_MM1S_IRF4.bedGraph >  ./bedgraph/KB_MM1S_IRF4_sort.bedGraph
bedGraphToBigWig ./bedgraph/KB_MM1S_IRF4_sort.bedGraph hg19.chrom.sizes ./bw/tfs/KB_MM1S_IRF4.bw 

sort -k1,1 -k2,2n  ./bedgraph/KB_MM1S_p300.bedGraph >  ./bedgraph/KB_MM1S_p300_sort.bedGraph
bedGraphToBigWig ./bedgraph/KB_MM1S_p300_sort.bedGraph hg19.chrom.sizes ./bw/tfs/KB_MM1S_p300.bw 

sort -k1,1 -k2,2n  ./bedgraph/KB_MM1S_H3K27ac.bedGraph >  ./bedgraph/KB_MM1S_H3K27ac_sort.bedGraph
bedGraphToBigWig ./bedgraph/KB_MM1S_H3K27ac_sort.bedGraph hg19.chrom.sizes ./bw/histone/KB_MM1S_H3K27ac.bw 

sort -k1,1 -k2,2n  ./bedgraph/KB_MM1S_H3K18ac.bedGraph >  ./bedgraph/KB_MM1S_H3K18ac_sort.bedGraph
bedGraphToBigWig ./bedgraph/KB_MM1S_H3K18ac_sort.bedGraph hg19.chrom.sizes ./bw/histone/KB_MM1S_H3K18ac.bw 

sort -k1,1 -k2,2n  ./bedgraph/KB_MM1S_H3K27me3.bedGraph >  ./bedgraph/KB_MM1S_H3K27me3_sort.bedGraph
bedGraphToBigWig ./bedgraph/KB_MM1S_H3K27me3_sort.bedGraph hg19.chrom.sizes ./bw/histone/KB_MM1S_H3K27me3.bw 

#### sorting and merge peaks from IRF4 and p300

sort -k1,1 -k2,2n ./beds/31_0EX7_0229Kronos_MM1S-No-Treatment-1_IRF4_hs_i45_peaks.bed > ./beds/IRF4_1_sort.bed
sort -k1,1 -k2,2n ./beds/32_0EX8_0229Kronos_MM1S-No-Treatment-2_IRF4_hs_i48_peaks.bed > ./beds/IRF4_2_sort.bed
sort -k1,1 -k2,2n ./beds/33_0EX9_0229Kronos_MM1S-No-Treatment-1_p300_hs_i49_peaks.bed > ./beds/p300_1_sort.bed
sort -k1,1 -k2,2n ./beds/34_0EXA_0229Kronos_MM1S-No-Treatment-2_p300_hs_i50_peaks.bed > ./beds/p300_2_sort.bed

bedtools merge -i ./beds/IRF4_1_sort.bed -i ./beds/IRF4_2_sort.bed > ./beds/IRF4.bed
bedtools merge -i ./beds/p300_1_sort.bed -i ./beds/p300_2_sort.bed > ./beds/p300.bed
cat ./beds/p300.bed ./beds/IRF4.bed > ./beds/p300_IRF4_merge.bed

##################

computeMatrix reference-point -S ./bw/tfs/* -R  ./beds/p300_IRF4_merge.bed -b 3000 -a 3000 -o computeMatrix.gz -p 7 -bs 50 --referencePoint center

plotHeatmap -m computeMatrix.gz -out figure_1_chip_tf_heatmap.pdf  --plotFileFormat 'pdf' --sortUsing mean --colorList 'white,magenta' 'white,cyan' 'white,black' 'white,black' 'white,black' 'white,black' 'white,black' --interpolationMethod 'nearest' --outFileSortedRegions ./sorted_heatmap.bed

computeMatrix reference-point -S ./bw/histone/* -R  ./sorted_heatmap.bed -b 3000 -a 3000 -o computeMatrix.gz -p 7 -bs 50 --referencePoint center

plotHeatmap -m computeMatrix.gz -out figure_1_chip_histone_heatmap.pdf  --plotFileFormat 'pdf' --colorList 'white,black' 'white,black' 'white,black'  --interpolationMethod 'nearest'
