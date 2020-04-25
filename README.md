# GMStool
  - GWAS-based Marker Selection Tool for Genomic Prediction from Genomic Data


## 1. Introduction

  - GMStool is a tool of selecting an optimal marker-set for predicting a phenotype. 
  - This tool is based on genome-wide association study, and heuristically searches optimal markers for the phenotype prediction through statistic and machine-learning methods (RR-BLUP and Random forest). Then, it performs genomic predictions using several statistical and machine/deep learning models (RR-BLUP, Random forest, DNN, and CNN) and finally presents the best prediction model with the optimal marker-set. 
  - GMStool consists of three phases, preparation, marker selection, and final modeling. Three input files, genotype, phenotype, and GWAS result files, are prepared in the preparation phase, and optimal markers are then selected in the marker selectinon phase ("GMStools.MS.v1.R" or "GMStools.MS.MultiThreading.v1.R" scripts). Last, final prediction using the optimal marker-set is conducted in the final modeling phase ("GMStools.FM.v1.R" script).
  
  
## 2. Flow-chart
  - Flow-chart of GMStool is as follows:

![그림11111](https://user-images.githubusercontent.com/49300659/80271666-40c93f00-86fd-11ea-81f9-08dc33b51163.jpg)


## 3. Run
  #### 3-1. Preparation phase
 - GMStool essentially requires three input files, genotype, phenotype, and GWAS result files. Preset file is optional.
 - "Genotype file" consists of markers (rows) and samples (columns), and genotypes are coded as -1, 0, 1, and 2 along missing, homozygous reference, heterozygous, and homozygous alternative genotypes. 
 - "Phenotype file" consists of samples (rows) and phenotype values (a column). Only one phenotype column is acceptable to GMStool, and an output directory is created based on the phenotype column name (see https://github.com/JaeYoonKim72/GMStool/tree/master/Results).
 - "GWAS result file" consists of SNPIDs (marker names), chromosome number, physical position, and p-value columns, in order. Additional columns may be present in the GWAS results file, but these four columns must be organized in order. 
 - "Preset file" means a list of markers that must be selected, and consists of a column with marker names (optional). 

![그림111111111](https://user-images.githubusercontent.com/49300659/80273485-57779200-870d-11ea-8cd0-1297dd98b052.jpg)


  #### 3-2. Marker selection phase

  - Marker selection phase is executed by either "GMStools.MS.v1.R" or "GMStools.MS.MultiThreading.v1.R" scripts. The only difference between the two scripts is whether multithreading is performed. All other options are identical.
  
    
    Usage: 
    GMStools.MS.v1.R -m [METHOD] -g [GENO] -p [PHENO] -gw [GWAS] -i [INFO] -pre [PRESET] 
                     -cv [CV] -a [ACC] -d [INCREMENT] -is [INTIAL SNPs] -ss [SNPS_SELECTED] 
                     -gpu [GPU_USAGE] -all [ALL_SNPs] -t [TIME]
                             
        Description of arguments:
             -m METHOD,         Selection method (RRblup, RF, or RRblup_RF).
             -g GENO,           Genotype file (Essential).
             -p PHENO,          Phenotype file (Essential).
             -gw GWAS,          GWAS result file. If GWAS file is not provided, GMStool calculates marker effects internally (Essential or optional).
             -i INFO,           Marker information file. Required if GWAS file is not provided (Optional; Default NULL).
             -pre PRESET,       Marker list to be selected in advance (Optional; Default NULL).
             -cv CV,            The number of cross validation (Default 3).
             -a ACC,            Goal of correlation rate (Default 0.9).
             -d INCREMENT,      Increament of correlation rate in marker selection (Default 0.001).
             -is INITIAL_SNPS,  The number of initial markers to be selected (>=2) (Default 5).
             -ss SNPS_SELECTED, The number of markers to be selected at one time (Default 1).
             -gpu GPU_USAGE,    If TRUE, RR-BLUP is calculated using GPU (Default FALSE).
             -all ALL_SNPs,     If TRUE, correlation rates of all markers for validation sets are calculated, but it takes a lot of time (Default FALSE).
             -t TIME,           Runtime cut-off per each CV (Default 1 hour).
  
  
  - -m option specifies the selection method to be used, and RRblup and RF can be selected. If you want to use both RRblup and RF, put "_" between the two methods and specify the -m option to "RRblup_RF".
  - -g and -p options are mandatory, and specify the genotype and phenotype files prepared in the previous phase.
  - -gw option specifies the gwas result file obtained in the previous phase. If this option is not given by the user, GMStool internally estimates the effects of markers for conducting the marker selection. When RR-BLUP is selected as the selection method, marker effects are derived from the coefficients of genotype variables of the RR-BLUP model. Also, when RF is selected, variable importance values of the RF model are estimated as marker effects. Although GMStool has the functions to estimate marker effects internally, it is recommended to use a separate GWAS result file with -gw option.
  - -i option means the marker information file, and is used when the -gw option is not given. 
  - -pre option specifies markers that must be selected.
  - -cv option means k value in k-fold cross validation, and indicates the number of cross validation.
  - -a option specifies the target accuracy of the markers to be selected.
  - -d option is an increment value of accuracy, and a marker to be selected must be higher than the accuracy of the previous marker plus the increment value.
  - -is option means the number of top markers to select initially from the priority of GWAS markers. If the preset option is defined (-pre), the -is option is ignored and the preset markers are considered initial markers.
  - -ss option indicates the number of markers to select at one time in the marker selection algorithm. It is recommended to select one marker at one time.
  - -gpu option determines whether to use the GPU when calculating the RR-BLUP method. Depending on the GPU and system settings, it may not be possible in some computation environments.
  - -all option determines whether to calculate the accuracy of all markers for the validation-set in each CV. This accuracy can be used as a reference for the minimum accuracy that the finally selected markers should have.
  - -t option 
  
  





        Example of marker selection not using multithreading
    
          Rscript GMStools.MS.v1.R \   
    
                     -m RRblup_RF \                       # Chose the selection methods (RRblup, RF, or RRblup_RF)
                         
                     -g ExampleData/Ex_genotype.txt \     # Genotype file
                         
                     -p ExampleData/Ex_phenotype.txt \    # Phenotype file
                         
                     -gw ExampleData/Ex_gwas.txt \        # GWAS result file
                         
                     -cv 3 \                              # Cross validation value
                         
                     -a 0.9 \                             # Target accuracy 
                          
                     -is 5 \                              # The number of initial SNPs
                                                 
                     -t 4                                 # Computational time for each CV



        Example of marker selection using multithreading
    
          Rscript GMStools.MS.MultiThreading.v1.R \   
    
                     -m RRblup_RF \                       # Chose the selection methods (RRblup, RF, or RRblup_RF)
                         
                     -g ExampleData/Ex_genotype.txt \     # Genotype file
                         
                     -p ExampleData/Ex_phenotype.txt \    # Phenotype file
                         
                     -gw ExampleData/Ex_gwas.txt \        # GWAS result file
                         
                     -cv 3 \                              # Cross validation value
                         
                     -a 0.9 \                             # Target accuracy 
                          
                     -is 5 \                              # The number of initial SNPs
                                                 
                     -t 4                                 # Computational time for each CV
                     
                     


  #### 3-3. Final modeling phase
  
    Usage: 
         GMStools.FM.v1.R -m [MODEL] -d [DIR] -gw [GWAS] -i [INFO] -pe [PERMUTATION] -gpu [GPU_USAGE] -t [TIME]

    Description of arguments:
         -m MODEL,        Prediction model (RRblup, RF, DNN, or CNN).
         -d DIR,          Result directory of marker selection (Essential).
         -gw GWAS,        GWAS result file. If GWAS file is not provided, marker information file should be provided (Essential or optional).
         -i INFO,         Marker information file. Required if GWAS file is not provided (Optional; Default NULL).
         -pe PERMUTATION, The number of permutations per each modeling (Default 50).
         -gpu GPU_USAGE,  If TRUE, DNN and CNN are calculated using GPU (Default FALSE).
         -t TIME,         Runtime cut-off for permutatios of each modeling (Default 1 hour).



    Example of final modeling
    
        Rscript GMStools.FM.v1.R \   
    
                     -m RRblup_RF_DNN_CNN \               # Chose the prediction models (RRblup, RF, DNN, or CNN)
                         
                     -d Results/Phenotype_RRblup_RF_PN_CV3_Ini5_Sel1_with_gwas/ \ # The path of result directory of marker selection
                         
                     -gw ExampleData/Ex_gwas.txt \        # GWAS result file
                         
                     -pe 50 \                             # The number of permutations for each modeling
                         
                     -gpu TRUE \                          # Whether to use the GPU when modeling DNN or CNN 
                                                                           
                     -t 1                                 # Runtime cut-off for permutatios of each modeling
                     
                     
   #### 3-4. Running screen            
  
  
![그림5](https://user-images.githubusercontent.com/49300659/80275712-0754fb80-871e-11ea-8df3-01e375c511fd.jpg)

  

## 4. Results

  #### 4-1. Marker selection phase
  
  ![re1](https://user-images.githubusercontent.com/49300659/80275083-fa360d80-8719-11ea-960c-bb9c7ef82dc1.jpg)


  #### 4-2. Final modeling phase
  
  ![re12](https://user-images.githubusercontent.com/49300659/80275247-02db1380-871b-11ea-86f0-dc1c0b6b696a.jpg)
  ![re123](https://user-images.githubusercontent.com/49300659/80275249-04a4d700-871b-11ea-898d-0db4e7d1d0e5.jpg)

## 5. Requirement

GMStool basically requires R version 3.6.1 or higher, and needs 13 additional libraries.

- data.table, tidyverse, ggplot2, ggpmisc, caret, randomForest, rrBLUP, tensorflow, keras, gpuR, doMC, foreach, and iterators


## 6. Contact

jaeyoonkim72@gmail.com

likemun@gmail.com


## 7. Citation

- Paper is under review.
