# *EndHiC* 
EndHic is a fast and easy-to-use Hi-C scaffolding tool, using the Hi-C links from contig end regions instead of whole contig regions to assemble large contigs into chromosomal-level scaffolds. 

EndHiC takes the HiC-pro's bin matrix results as input data. After running HiC-pro, the recommended EndHiC usage for most users is to run endhic.pl or endhic_iterate.pl. When your contig assembly is quite good, then endhic.pl [one round of EndHiC] is able to finish the job; When your contig assembly is relatively fragmental, then endhic_iterate.pl [multiple rounds of EndHiC] should be used. How many rounds is needed depend on the contig assembly level, and more fragmental contigs need higher rounds of EndHiC.

Contents
========
- [Example Data](#example-data)
- [Usage instructions](#usage-instructions)
- [EndHiC programs](#endhic-programs)
	- [Basic usage](#basic-usage)
	- [Basic pipeline](#basic-pipeline)
	- [Standard pipeline](#standard-pipeline)
	- [Iterative pipeline](#iterative-pipeline)
- [EndHiC output sub-directory and files](#endhic-output-sub-directory-and-files)
- [Pre-EndHiC programs](#pre-endhic-programs)
- [Post-EndHiC programs](#post-endhic-programs)
- [Accuracy verifying programs](#accuracy-verifying-programs)
- [Reference](#reference)
- [Contact](#contact)

## Example Data:

Two example data are included in EndHiC package:
```
  git clone git@github.com:fanagislab/EndHiC.git

  cd EndHiC/z.testing_data/Arabidopsis_thalina
  sh work.sh 

  cd EndHiC/z.testing_data/Cichorium_intybus
  sh work.sh
``` 
>**Note**
> The Arabidopsis_thalina testing data shows the usage of endhic.pl on long-continuous contig assembly, while the Cichorium_intybus testing data shows the usage of endhic_iterate.pl with relatively shorter contig assembly.

An example data for detecting the assembly errors in the contigs is also included in EndHiC package:
```
  cd EndHiC/z.testing_data/Detect_errors_in_contigs
  sh work.sh
```

## Usage instructions:

1. Input files: (example, human hifiasm + hic-pro)
	- `hifiasm.fa.len`: Includes two column: contig_id contig_length	

	- `humanHiC_100000_abs.bed`: Generated by Hic-pro, 100-kb bins, bed format file

	- `humanHiC_100000.matrix`: Generated by Hic-pro, 100-kb bins, raw matrix file

	- `humanHiC_100000_iced.matrix`: Generated by Hic-pro, 100-kb bins, normalized matrix file


## EndHiC programs: 
the ranks shows invoking relationship
```
endhic_iterate.pl

	endhic.pl							

		endhic_ctgEnd_pipeline.pl

			ctgContact_from_ctgEndContacts.pl
			
			turningpoint_by_lineartransform.pl
			
			scaffold_by_trueCtgContact.pl
			
			cluster_and_classify_GFA.pl
			
			order_and_orient_GFA.pl
```

### Basic usage: 

run endhic with specified contig end size and specified contact cutoff
- **Step 1:**  calculate the HiC contact values among contigs, using Hi-C links data from fixed-size contig ends
	```
	ctgContact_from_ctgEndContacts.pl \
		--binsize 100000 \
		--binnum 10 hifiasm.fa.len \
		humanHiC_100000_abs.bed \
		humanHiC_100000.matrix \
		> humanHiC_100000.matrix.100000.10.CtgContact
	```

- **Step 2:** adjust the contig contacts, and perform linear transformation, to find the turning point
	```
	turningpoint_by_lineartransform.pl \
		humanHiC_100000.matrix.100000.10.CtgContact \
		> humanHiC_100000.matrix.100000.10.CtgContact.adjustTransform \
		2> humanHiC_100000.matrix.100000.10.CtgContact.turningPoint
	```
- **Step 3:** build contig graph by assigning links to contigs whose contact is larger than a given cutoff, and also satisfy reciprocal best requirement
	```
	scaffold_by_trueCtgContact.pl \
		--contacts 147.07 \
		--reciprocalmax  hifiasm.fa.len  \
		humanHiC_100000_iced.matrix.100000.10.CtgContact \
		> humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa 
	```

- **Step 5:** Identify linear and circular topology in the contig graph
	```
	cluster_and_classify_GFA.pl \
		humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa \
		> humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa.topology
	```

- **Step 6:** Output cluster results with order and orientation information
	```
	order_and_orient_GFA.pl \
		--size 2000000 \
		humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa \
		humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa.topology \
		> humanHiC_100000_iced.matrix.100000.10.CtgContact.overCutoff.1.0.reciprocalMax.gfa.cluster
	```

### Basic pipeline:

- **Step 1:** run endhic with specified contig end size, in various automatically determined contact cutoff, using Hic-pro **raw** matrix data
	```
	endhic_ctgEnd_pipeline.pl \
		--binsize 100000 \
		--binnum 10 \
		hifiasm.fa.len \
		humanHiC_100000_abs.bed \
		humanHiC_100000.matrix
	```

- **Step 2:** run endhic with specified contig end size, in various automatically determined contact cutoff, using Hic-pro **normalized** matrix data
	```
	endhic_ctgEnd_pipeline.pl \
		--binsize 100000 \
		--binnum 10 \
		hifiasm.fa.len \
		humanHiC_100000_abs.bed \
		humanHiC_100000_iced.matrix
	```

### Standard pipeline: 

>**Note**
> run only one round of endhic.pl, when contig assembly is quite good

run endhic with various contig end size, in various automatically determined contact cutoff, using Hic-pro raw and normalized matrix data. At most cases, this can generate chromosome-level scaffolds

```
endhic.pl  \
	hifiasm.fa.len \
	humanHiC_100000_abs.bed \
	humanHiC_100000.matrix \
	humanHiC_100000_iced.matrix
```

### Iterative pipeline: 

>**Note**
> run multiple rounds of endhic.pl, when contig assembly is not so good

If a single run of endhic.pl can't finish the scaffolding task, i.e. the number of resulting clusters is more than that of chromosomes, iterative running of endhic.pl is recommended. In each loop, the contig end size is increasing. In this way, the problems caused by the repeat sequences on the contig ends will be overcomed.

Using default parameters of endhic_iterate.pl
```
endhic_iterate.pl  \
	--rounds 3  \
	--binnumstep 5 \
	hifiasm.fa.len \
	humanHiC_100000_abs.bed \
	humanHiC_100000.matrix \
	humanHiC_100000_iced.matrix
```

For more shorter contigs, try to run more rounds with smaller increasing of contig end sizes
```
endhic_iterate.pl  \
	--rounds 15 \
	--binnumstep 1 \
	hifiasm.fa.len \
	humanHiC_100000_abs.bed \
	humanHiC_100000.matrix humanHiC_100000_iced.matrix
```
## EndHiC output sub-directory and files


In **01.contig_end_contact_results/**

- `humanHiC_100000.matrix.*.CtgContact`
    Contig contact file, with 7 columns (#CtgId1 CtgId2  EndContact Ctg1Pos Ctg2Pos UsedBinNum1     UsedBinNum2)

- `humanHiC_100000.matrix.*.CtgContact.adjustTransform`
    Contig contact, adjusted, and linear transformed, to find the turning point

- `humanHiC_100000.matrix.*.CtgContact.turningPoint`
    Automatically inferred turning point, which will be used as the basic value for the contig contact cutoff

In **02.GFA_contig_graph_results/**

- `humanHiC_100000.matrix.*.CtgContact.overCutoff.1.0.gfa`
    Contig graph in GFA format, contact value over cutoff, can be viewed in Bandage software

- `humanHiC_100000.matrix.*.CtgContact.overCutoff.1.0.reciprocalMax.gfa`
    Contig graph in GFA format,  contact value over cutoff, and satisfy reciprocal best, can be viewed in Bandage software

In **03.cluster_order_orient_results/**

- `humanHiC_100000.matrix.*.CtgContact.overCutoff.1.0.reciprocalMax.gfa.topology`
    Topology of the contig graph, identify linear or circular groups

- `humanHiC_100000.matrix.*.CtgContact.overCutoff.1.0.reciprocalMax.gfa.cluster`
    Scaffold results, including cluster, order, and orientation information


In **04.summary_and_merging_results/**

- `z.EndHiC.A.results.summary`
    Summary and analysis results for the first loop, merging all the raw and iced results

- `z.EndHiC.A.results.summary.cluster`  [Final EndHiC Result]
    Final scaffold results with high robustness, merging all the raw and iced results
    This is recommeded to be the final endhic result.

- `z.EndHiC.A.results.summary.cluster.gfa` 
    GFA format of the final scaffold results, which can be graphically viewed in Bandage software


Instruction of ***.summary** file:
  - **Part 1**: Number of clusters under each condition
  - **Part 2**: Statistics of all Cluster units
  - **Part 3**: Statistics of merged Cluster units
  - **Part 4**: Statistics of stable (high frequency) cluster units
  - **Part 5**: Statistics of stable cluster units (redundant short contigs removed)
  - **Part 6**: Included contigs, total number, total length 


Format of ***.cluster file**: 
  - **column 1**: Cluster id, sorted by cluster length   
  - **column 2**: Number of contigs included in this cluster 
  - **column 3**: Cluster length, total length of contigs in this cluster 
  - **column 4**: robustness, i.e. appearance times in the results from various contig end sizes and contact cutoffs
  - **column 5**: Included contigs with order and orientation, separated by ";", and "+-" means strands
             e.g. ptg000046l-;ptg000079l+;ptg000058l-;ptg000047l+ (equivalent to ptg000047l-;ptg000058l+;ptg000079l-;ptg000046l+)



## Pre-EndHiC programs

- **Step 1:** draw HiC heatmap for contigs, helpful to find assembly errors in contigs, each point stands for window size 10*100000 = 1 M
	```
	matrix2heatmap.py \
		humanHiC_100000_abs.bed \
		humanHiC_100000.matrix \
		10
	```

- **Step 2:** mapping the unitigs (p_utg) to contigs (p_ctg), only tested for Hifiasm contigs:
	```
	perl ./map_utg_to_ctg.pl \
		human.p_ctg.noseq.gfa \
		human.p_utg.noseq.gfa \
		> human.utg_to_ctg.map \
		2> human.utg_to_ctg.map.gfa
	```

- **Step 3:** Identify the assembly errors in contigs by Hi-C heatmaps and unitig breaks:

	use both the HiC data and the utg_to_ctg.map data
	```
	perl asm_error_check.pl \
		human_100000_abs.bed \
		human_100000.matrix \
		human.utg_to_ctg.map \
		> human.assmebly_errors.position
	```

	>**Note**
	>only use the HiC data, when the utg_to_ctg.map is not available

	```
	perl asm_error_check.pl \
		human_100000_abs.bed \
		human_100000.matrix  \
		> human.assmebly_errors.position
	```

	The **result** file contains 10 columns:
	- **Ctg_id**: contig id
	- **Error_loc**: assembly error position inferred from HiC data
	- **Inter_contact**: the count of HiC links between the two bins (500kb-apart) crossing the assembly error position
	- **Inter_cutoff**: the cutoff of Inter_contact, which is 1/10 of the median values from all the inter-bins (500kb-apart)
	- **Intra_contact1**: the count of HiC links within the left bin crossing the assembly error position
	- **Intra_contact2**: the count of HiC links within the right bin crossing the assembly error position
	- **Intra_cutoff**: the cutoff of Intra_contact1 and Intra_contact2, which is 1/10 of the median values of all the intra-bins
	- **Break_pos**: assembly error position inferred from unitig break point
	- **Utg1_id**: the left unitig crossing the assembly error position 
	- **Utg2_id**: the right unitig crossing the assembly error position 



## Post-EndHiC programs


- **Step 1:** Convert to AGP and Fasta format files

	convert cluster format file to AGP format file
	```
	cluster2agp.pl \
		z.EndHiC.A.results.summary.cluster  \
		hifiasm.fa.len  \
		> scaffolds.agp
	```
	convert AGP format file to Fasta format file
	```
	agp2fasta.pl \
		scaffolds.agp  \
		hifiasm.fa \
		> scaffolds.fa
	```
- **Step 2:** Draw Hi-C heatmaps for the EndHiC scaffolds
	convert the contig bed file to cluster bed file
	```
	cluster2bed.pl \
		humanHiC_100000_abs.bed  \
		z.EndHiC.B.results.summary.cluster \
		> clusterB_100000_abs.bed \
		> 2> clusterB.id.len
	```
	draw HiC heatmap for endhic scaffolds, each point represents window size 20*100000 = 2 M  
	```
	matrix2heatmap.py \
		clusterB_100000_abs.bed \
		humanHiC_100000.matrix  \
		20
	```

- **Step 3:** Mapping unclustered short contigs to each cluster

	calculate the HiC contact values among contigs, using Hi-C links data from half contig(i.e. max contig end size)
	```
	perl ../../../ctgContact_from_ctgEndContacts.pl  \
		--binsize 100000  \
		--binnum -1  \
		hifiasm.fa.len  \
		humanHiC_100000_abs.bed  \
		humanHiC_100000.matrix \
		> humanHiC_100000.matrix.halfContig.ctgContact
	```

	normalize the contig contact by the used bin numbers, only keep the max contact from head vs head, head vs tail, tail vs head, tail vs tail comparisons
	```
	ctgContact_normalize_distance.pl  \
		--normalize  humanHiC_100000.matrix.halfContig.ctgContact \
		> humanHiC_100000.matrix.halfContig.ctgContact.normalize
	```

	mapping the unclustered short contigs into Endhic clusters, with specified cutoff
	```
	shortCtgs_to_cluster.pl \
		--contact 1 \
		--times 2   \
		z.EndHiC.A.results.summary.cluster  \
		hifiasm.fa.len  \
		humanHiC_100000.matrix.halfContig.ctgContact.normalize \
		> shortCtgs.mapped.to.clusters.list
	```

- **Step 4:** Convert EndHiC result to juicebox compatible file formats, which can be viewed in Juicebox

	convert EndHiC .cluster file into juicebox .assembly file
	```
	cluster_to_juciebox_assembly.pl \
		contigs.fa.len \
		z.EndHiC.A.results.summary.cluster \
		> draft.assembly
	```
	
	index the contig sequence file
	```
	bwa index \
		draft.fa
	```
	
	generate the enzyme cutting sites file draft_MboI.txt
	```
	juicer/misc/generate_site_positions.py \
		MboI \
		draft draft.fa
	```
	
	generate the Hi-C reads alignment file aligned/merged_nodups.txt
	prepare data: put the HiC reads under ./fastq/; put the contig sequence file and index files under ./reference/;
	```
	juicer/CPU/juicer.sh \
		-S early \
		-g draft \
		-s MboI \
		-z ./references/draft.fa \
		-y ./draft_MboI.txt \
		-p ./references/draft.fa.size \
		-t 50 \
		-D juicer/CPU
	```
	
	generate the hic input file draft.hic for viewing in juicebox
	```
	3d-dna/visualize/run-assembly-visualizer.sh \
		draft.assembly \
		merged_nodups.txt  
	```
	
	For more instructions, please refer to the help pages of juicer and juicebox.


## Accuracy verifying programs

- **Step 1:** Run endhic using the max contact values from bin pairs of two compared contigs
	```
	endhic_maxBin_pipeline.pl \
		--binsize 1000000  \
		hifiasm.fa.len \
		humanHiC_100000_abs.bed \
		humanHiC_100000.matrix
	```
	
- **Step 2:** Apply Hierarchical clustering algorithm with contig distance converted from Hi-C contact values derived from contig end regions

	normalize the contig contact by the used bin numbers, and converted to distance values ranging from 0 to 1
	```
	ctgContact_normalize_distance.pl \
		humanHiC_100000.matrix.100000.10.CtgContact \
		> humanHiC_100000.matrix.100000.10.CtgContact.distance
	```

	generate all the middle procedure results of the Hierarchical clustering algorithm
	```
	hcluster_contigs.pl \
		--verbose \
		-type min \
		humanHiC_100000.matrix.100000.10.CtgContact.distance \
		hifiasm.fa.len \
		> humanHiC_100000.matrix.100000.10.CtgContact.distance.hcluster.one \
		2> humanHiC_100000.matrix.100000.10.CtgContact.distance.hcluster
	```

	Find the suitable stop loop, which represents correct chromosomes, by giving expected crhomosome number and minimum chromosome length cutoff
	```
	hcluster_suitable_stop.pl \
		--chr_num  23 \
		--chr_len  20000000  \
		humanHiC_100000.matrix.100000.10.CtgContact.distance.hcluster  \
		>  humanHiC_100000.matrix.100000.10.CtgContact.distance.hcluster.need
	```


- **Step 3:** Compare two clusters

	only compare the clustering information, not consider order and orientation information
	```
	cluster_compare.pl  \
		human.contigs.minimap2.cluster  \
		z.EndHiC.A.results.summary.cluster  \
		> z.EndHiC.A.results.summary.cluster.vs.ref
	```

## Reference

Sen Wang, Hengchao Wang, Fan Jiang, Anqi Wang, Hangwei Liu, Hanbo Zhao, Boyuan Yang, Dong Xu, Yan Zhang, Wei Fan. EndHiC: assemble large contigs into chromosomal-level scaffolds using the Hi-C links from contig ends. (2021)  https://arxiv.org/abs/2111.15411v1


## Contact 
Contact any of these authors for help:
**Wei Fan**, 0000-0001-5036-8733  fanwei@cass.cn  or fanweiagis@126.com   
**Sen Wang**,  0000-0001-9793-4472  wangsen1993@163.com   
**Hengchao Wang**, 0000-0002-8754-4195  wanghengchao000@qq.com  
**Fan Jiang**, 0000-0003-1359-0970	 greatjf@163.com  
**Yan Zhang**,	0000-0003-2281-7807 milrazhang@163.com  
