# 

Output of `retro_hunter` :

One row per retrogene called with:
1. gene anme
2. number of junctions
3. number of spliced reads
4. number of split pairs
5. sample library size
6. a tab delimated list of junctions. Each junction information will be comma delimited with the following information: 
  * junction position (chrom:intron_start-intron_end)
  * number of spliced reads across the junction
  * number of split pairs across the junction (*)
  * left intron coverage/left exon coverage. - we expect to see a large drop in coverage across exon boundaries for retrogenes 
  * right intron coverage/left intron coverage.
