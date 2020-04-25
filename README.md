# GMStool
  - GWAS-based Marker Selection Tool for Genomic Prediction from Genomic Data

## 1. Introduction

GWAS GAPIT pipeline was based on GAPIT libaraies developed by Lipka and Tang, and constructed by Jae-Yoon Kim.

This analysis pipeline uses a VCF file as input file and performs a genome-wide association study.

Source code was written in Python and R languages and supported on windows and linux platforms.


## 2. Flow-chart

The flow-chart is as follows:

![그림11111](https://user-images.githubusercontent.com/49300659/80271666-40c93f00-86fd-11ea-81f9-08dc33b51163.jpg)


## 3. Run
  #### 3-1. Preparation phase





  #### 3-2. Marker selection phase

    Usage: GMStools.MS.v1.R -m [METHOD] -g [GENO] -p [PHENO] -gw [GWAS] -i [INFO] -pre [PRESET] 
                            -cv [CV] -a [ACC] -d [DELTA] -is [INTIAL SNPs] -ss [SEL_SNPS] 
                            -gpu [GPU_USE] -all [ALL_SNPs] -t [TIME]


    Description of arguments:
         -m METHOD,        Selection method (RRblup, RF, or RRblup_RF).
         -g GENO,          Genotype file (Essential)
         -p PHENO,         Phenotype file (Essential)
         -gw GWAS,         GWAS result file. If GWAS file is not provided, GMStool calculates marker effects internally.
         -i INFO,          Marker information file. Required if GWAS file is not provided (Optional; Default NULL)
         -pre PRESET,      Marker list to be selected in advance (Optional; Default NULL)  
         -cv CV,           The number of cross validation (Default 3)
         -a ACC,           Goal of correlation rate (Default 0.9)
         -d DELTA,         Increament of correlation rate for selection (Default 0.001)
         -is INITIAL_SNPS, The number of initial markers to be selected (>=2) (Default 5)
         -ss SEL_SNPS,     The number of markers to be selected at one time (Default 1)
         -gpu GPU_USE,     If TRUE, RR-BLUP is calculated using GPU (Default FALSE)
         -all ALL_SNPs,  If TRUE, correlation rates of all markers for validation sets is calculated, but it takes a lot of time (Default FALSE)
         -t TIME,          Runtime cut-off for each CV (Default 1)



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
  

    Example of final modeling
    
        Rscript GMStools.FM.v1.R \   
    
                     -m RRblup_RF_DNN_CNN \               # Chose the prediction models (RRblup, RF, DNN, or CNN)
                         
                     -d Results/Phenotype_RRblup_RF_PN_CV3_Ini5_Sel1_with_gwas/ \ # The path of result directory of marker selection
                         
                     -gw ExampleData/Ex_gwas.txt \        # GWAS result file
                         
                     -pe 50 \                             # The number of permutations for each modeling
                         
                     -gpu TRUE \                          # Whether to use the GPU when modeling DNN or CNN 
                                                                           
                     -t 1                                 # Computational time for each each modeling
                     
  
  
  

## 4. Results

Result files are provided with a total of 25 files including a result table and the following 3 images.



## 5. Requirement

GMStool basically requires R version 3.6.1 or higher, and needs 13 additional libraries.

- data.table, tidyverse, ggplot2, ggpmisc, caret, randomForest, rrBLUP, tensorflow, keras, gpuR, doMC, foreach, and iterators


## 6. Contact

jaeyoonkim72@gmail.com

likemun@gmail.com


## 7. Citation

Pulication is under review.
