# 

Output of `retro_hunter` :

One row per retrogene called with:
1. **gene name**
2. **number of junctions** where split reads supporting a retrogene have been identified. This should match the number of junctions listed in 6. Junctions must have at least one split read and >= 1 more split read or pair to be counted.
3. **number of split reads** number of "spliced" reads in total across all junctions identified for the gene.
4. **number of split pairs** number of split read paris in total across all junctions identified for the gene. Note: if this is zero for all genes, it is likely that the wrong reference was used (eg. hg19 instead of hg38). Split reads would still be reported even with the wrong reference.
5. **library size** total number of reads sequenced for the sample
6. **tab delimated list of junctions** Each junction information will be comma delimited with the following information: 
  * **junction position (chrom:intron_start-intron_end)** Note: positions might be out by a few bases. 
  * **number of spliced reads** across the junction
  * **number of split pairs** across the junction. Note: this is only an estimate. Split pairs are assigned to junctions based on read distance to the closest exon boundary, which may not be correct. Exons which are close together in the genome (with small intron) could lead to split pair counts without a retrogene.
  * **left intron coverage/left exon coverage** where intron coverage is taken 4bp into the intron and exon coverage is taken 1 bp into the exon. We expect to see a large drop in coverage across exon boundaries for retrogenes 
  * **right intron coverage/left intron coverage**
