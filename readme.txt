partition.py 

	it takes haplotype or genotype data (see SampleData.txt) as input, and find blocks according to Fout Gamete Test.
	
	outputs include: 	BLOCKS	|	CORES	|	SUSPECT SNPS	|	MAXIMAL-K-COVER
	
	our blocks satisfy the following condition:
		A.	no recombination occurs within a block
		B.	every block is a maximal length local Perfect Phelogeny Tree
		C.	allowing overlap between blocks
		Biological meaning: each block contains several SNPs and inherits as a whole unit, which could be directly used to association study, infering genome history etc.
		
	our cores satisfy the following condition:
		A.	in dependent of the scanning order (Left-to-Right or Right-to-Left)
		B.	the number of cores is the minimal number of blocks to cover the entire dataset
		C.	no recombination occurs within a core
		D.	no overlap between cores
		Biological meaning: need to be elucidated
		
	our suspect SNPs satisfy the following condition:
		A.	once removed, the minimal number of blocks to cover the entire dataset will reduce
		Biological meaning: they are indication of genotyping error, gene conversion and homoplasy.
		
	our maximal-k-cover satisfies the following condition:
		A.	it is a set of blocks that are maximal in their overlap and minimal in their number to cover the entire dataset
		Biological meaning: need to be elucidated


to test use:
	python genotypeScan_sample_usage.py genSampleData.txt
	python haplotypeScan_sample_usage.py hapSampleData.txt
 
