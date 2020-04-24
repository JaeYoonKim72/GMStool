## Description of modules of GMStool


The directory named "src" and the R function files should always be with the files "GMStools.MS.v1.R" and "GMStools.MS.MultiThreading.v1.R".


  1. "load_files.R" is used in the preparation phase of "GMStools.MS.v1.R" and "GMStools.MS.MultiThreading.v1.R", and loads genotype, phenotype, and gwas result files 


  2. "GMS_main.R" is only used in the "GMStools.MS.MultiThreading.v1.R", and utilized as the main module for multithreading. 


  3. "RRb_func.R" and "RF_func.R" are the functions used in the marker selection and final modeling phases of "GMStools.MS.v1.R" and "GMStools.MS.MultiThreading.v1.R".


  4. "DNN_func.R" and "CNN_func.R" are the functions only used in the final modeling phase of "GMStools.MS.v1.R" and "GMStools.MS.MultiThreading.v1.R". These functions support GPU high-speed computation.


  5. "gpu.mixed.solve.R" is used with "RRb_func.R", and 

