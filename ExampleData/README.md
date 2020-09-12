## Description of example data 
  

1. "Ex_genotype.txt" 
   - This file consists of 12,301 SNP markers (row) and 200 samples (column).
   - SNP markers belong to chromosomes 1-12 and are named SNP1, SNP2, SNP 3, ... and SNP12301.
   - Samples are named Sample1, Sample2, Sample3 ..., and Sample200.
   - Markers are coded as -1, 0, 1, and 2 along missing, homozygous reference, heterozygous, and homozygous alternative genotypes. 



2. "Ex_phenotype.txt" 
   - This file is the phenotype data of the 200 samples, and consistis of two columns, the Sample ID and the phenotype value.
   - The phenotype name in the example file are named "Phenotype".



3. "Ex_gwas.txt"
   - This file is the result of a genome-wide association study between the 200 samples and their phenotypes.
   - This file consists of "SNP marker", "Chromosome", "Position", and "P.value" columns in order.
   


4. "Ex_marker_information.txt"
   - This file is the information of 12,301 SNP markers, and consists of "SNP marker", "Chromosome", and "Position" columns in order.
   - If a GWAS reuslt file is not provided ("Ex_gwas.txt"), GMStool internally estimates marker effects and performs marker selection and prediction based on the priority of the marker effects. At this time, this file is used to obtain the information of the markers.
   
   
