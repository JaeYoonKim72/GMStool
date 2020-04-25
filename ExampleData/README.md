## Description of example data 
  

Before using the example file, extract the "Ex_genotype.zip" file.


1. "Ex_genotype.txt" 
   - This file consists of 36,091 SNP markers (row) and 413 samples (column).
   - SNP markers belongs to chromsomes 1-12, and samples were named as Sample1, Sample2, Sample3 ..., and Sample413.
   - Markers are coded as -1, 0, 1, and 2 along missing, homozygous reference, heterozygous, and homozygous alternative genotypes. 



2. "Ex_phenotype.txt" 
   - This file is the phenotype data of 413 samples, and consistis of two columns, the Sample ID and the phenotype value.
   - The phenotype column in the example file was named "Phenotype".



3. "Ex_gwas.txt"
   - This file is the result of a genome-wide association study between the 413 sample samples and their phenotypes.
   - This file has 4 columns, and consists of "SNP marker", "Chromosome", "Position", and "P.value" columns in order.
   

4. "Ex_marker_information.txt"
   - This file is the information of 36,091 SNP markers, and consists of "SNP marker", "Chromosome", and "Position" columns in order.
   - If a GWAS reuslt file is not provided ("Ex_gwas.txt"), GMStool internally estimates marker effects and performs marker selection and prediction based on the priority of the marker effects. At this time, this file is used to obtain the information of the markers.
   
   
