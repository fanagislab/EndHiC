gzip -d combine_100000_iced.matrix.gz 
gzip -d combine_100000.matrix.gz

##run three rounds of EndHiC
perl ../../endhic_iterate.pl --rounds 3 --binnumstep 5  contigs_all.len combine_100000_abs.bed combine_100000.matrix combine_100000_iced.matrix 

##For the usage of other programs in the EndHiC package, please see examples in the Arabidopsis_thalina directory

##Cichorium intybus has 9 chromosomes. Using 3 rounds of EndHiC, we can get a near-chromosome-level scaffolds. To finally get all the chromosomes, manual curation with Hi-C heatmaps or guided by Juicebox are needed.

