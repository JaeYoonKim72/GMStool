##########
#GMStool v1 for Marker selecion, made by Seongmun Jeong and JaeYoon Kim in 2020-09-09.
##########

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-m", "--method", default = "RRB", help = "Selection method (RRB, BTS, RRB_BTS)")
parser$add_argument("-g", "--geno", default = NULL, help = "Genotype file")
parser$add_argument("-p", "--pheno", default = NULL, help = "Phenotype file")
parser$add_argument("-gw", "--gwas", default = NULL, help = "GWAS result file")
parser$add_argument("-i", "--info", default = NULL, help = "Marker information file. Required if GWAS file is not provided")
parser$add_argument("-t", "--test", default = NULL, help = "Test sample list file. This file contains the sample names of the test set")
parser$add_argument("-pre", "--preset", default = NULL, help = "Marker list to be selected in advance (Default NULL)")
parser$add_argument("-cv", "--cv", type = "integer", default = 5, help = "# of cross validation (>=3) (Default 3)")
parser$add_argument("-c", "--corr", type = "double", default = 0.95, help = "Traget correlation rate for the validation set (Defualt 0.95)")
parser$add_argument("-d", "--increment", type = "double", default = 0.00005, help = "Increament of correlation rate in marker selection (Default 0.00005)")
parser$add_argument("-is", "--initial_snps", type = "integer", default = 5, help = "# of initial SNPs to be selected (>=2) (Default 5)")
parser$add_argument("-ss", "--snps_selected", type = "integer", default = 1, help = "# of SNPs to be selected at one time (Default 1)")
parser$add_argument("-gpu", "--GPU_usage", default = FALSE, help = "If TRUE, RRB is calculated using GPU (Only available in Linux) (Default FALSE)")
parser$add_argument("-all", "--all_snps", default = FALSE, help = "If TRUE, correlation rate of all markers for the validation set in each CV is calculated, but it takes a lot of time (Default FALSE)")

#(0) Prepareing analysis
args <- parser$parse_args()
genofile <- args$geno
phenofile <- args$pheno
gwasfile <- args$gwas
prefile <- args$preset
testsfile <- args$test
infofile <- args$info
cv <- as.numeric(args$cv)
delta <- as.numeric(args$increment)
sel_snps <- as.numeric(args$snps_selected)
acc1 <- as.numeric(args$corr)
mmm <- args$method
gpu_use <- as.logical(args$GPU_usage)
ini_snps <- args$initial_snps
allm <- as.logical(args$all_snps)

if (is.null(infofile)){
    infofile1 = "not provided"
} else {
    infofile1 = infofile
}
if (is.null(prefile)){
    prefile1 = "not provided"
} else {
    prefile1 = prefile
}
if (is.null(gwasfile)){
    gwasfile1 = "not provided"
} else {
    gwasfile1 = gwasfile
}
if (is.null(genofile) | is.null(phenofile)){
  cat("Error1: invalide genotype and phenotype files.", '\n')
  quit(save = "no")
}
if (mmm %in% c("RRB", "BTS", "RRB_BTS", "BTS_RRB")==FALSE){
  cat("Error2: invalide selection method.", '\n')
  quit(save = "no")
}
cat("\n\n")
cat("Reading files ==================================================", "\n")
cat("1)  Genotype file is", genofile, "\n")
cat("2)  Phenotype file is", phenofile, "\n")
cat("3)  GWAS result file is", gwasfile1, "\n")
cat("4)  Marker information file is", infofile1, "\n")
cat("5)  Pre-selected SNP file is", prefile1, "\n")
cat("6)  Test sample list file is", testsfile, "\n\n")

cat("Reading selection arguments ====================================", "\n")
cat("7)  Prediction method is", mmm, "\n")
cat("8)  Usage of GPU is", gpu_use, "\n")
cat("9)  Target accuracy is", acc1, "\n")
cat("10) Cross validation k is", cv, "\n")
cat("11) Increament cut-off of correlation rate (delta) is", delta, "\n")
cat("12) # of initial SNPs to be selected is", ini_snps, "\n")
cat("13) # of SNPs to be selected is at one time", sel_snps, "\n")
Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02, gc.profiling = TRUE)
cat("================================================================", "\n\n")


#(0) Initial check
if (cv < 3) {
  cat("Error3: Cross validation argument is much small. The argument is recomended to more than 3.", '\n')
  quit(save = "no")
}

if (is.null(testsfile)) {
  cat("Error4: invalid test sample list file.", '\n')
  quit(save = "no")
}

if (mmm != "RRB" & ini_snps <= 4){
  cat(paste0("Error5: ", mm, " method needs the number of initial SNPs more than 5."), '\n')
  quit(save = "no")
}
if ((is.null(gwasfile) & mmm == "RRB") | (is.null(gwasfile) & mmm == "RRB_BTS")){
  cat("Warnning1: GWAS result file was not defined. Marker effects will be derived using a linear-mixed model.", '\n')
}
if (is.null(gwasfile) & mmm == "BTS"){
  cat("Warnning2: GWAS result file was not defined. 
                Marker effects will be derived using the Random forest model, but it takes a lot of time.", '\n')
}
if (is.null(gwasfile) & is.null(infofile)){
  cat("Warnning3: SNP information file was not defined. 
                Plots for SNP distribution across chromosomes is not generated.", '\n')
}
if (allm == TRUE){
  cat("Warnning4: Calculation of correlation rate for all markers takes a lot of time.", '\n')
}


#(1) Load assicated packages, and set gpu/cpu and threads.
##(1-1) Load associated packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
library(ggplot2)
library(ggpmisc)
library(data.table)
library(tidyverse)
library(caret)
library(randomForest)
options(warn=-1)
source("./src/RRb_func.R")
source("./src/RF_func.R")
source("./src/load_files.R")
##(1-2) Set gpu/cpu
if (mmm == "RRB" | mmm == "RRB_BTS") {
  library(rrBLUP)
  if (gpu_use) {
    library(gpuR)
    source("./src/gpu.mixed.solve.R")
    mixed.solve <- gpu.mixed.solve
  }
} 

#(2)Load and preprocess input datasets
load_data = load_files(genofile, phenofile, testsfile, gwasfile, infofile)
geno1 = load_data$genotype
phenotype = load_data$phenotype
ix = load_data$ix
ix1 = ix
pheno1 = load_data$pheno
trains = load_data$train
tests = load_data$test

#(3) Load pre-selected markers
if (!is.null(prefile)){
  preset_data <- fread(prefile, header=T, check.names=F, data.table=FALSE)
  preset_name <- as.vector(unlist(preset_data))
  PN = length(preset_name)
  INIbk = ini_snps
  cat(paste0("Preset: # of preset markers: ", length(preset_name), sep= ""), '\n')
  preset_fname <- intersect(preset_name, colnames(geno1))
  cat(paste0("Preset: # of preset markers filtered: ", length(preset_fname), sep= ""), '\n')
  write.table(preset_fname, file = "Preset.filtered_markers.txt", quote=F, row.names = FALSE, col.names = "Pre_markers_filtered")
  PC = paste("PY", PN, "-", length(preset_fname), sep="")
  if (length(preset_fname) >= ini_snps){
    ini_snps = 0
    cat("- Preset markers more than initial markers, so initial marker argument is ignored.", "\n\n")
  } else if ( (0 < length(preset_fname)) &  (length(preset_fname) < ini_snps)) {
    ini_snps = 0
    cat("- Preset markers less than inital markers, # of inital marker is regarded to # of preset markers.", "\n\n")
  } else {
    ini_snps = ini_snps
    preset_fname = NULL
    cat(paste0("- No preset markers due to filtering, so inital makrer argument is selected as like defiend ", ini_snps, "\n\n"), sep="")
  }
} else {
  ini_snps = ini_snps
  INIbk = ini_snps
  preset_fname = NULL
  cat(paste0("Initial marker: GWAS top ", ini_snps, "\n\n"), sep="")
  PC = "PN"
}
ini_snps_bk = ini_snps


#(4) Preparing train/validation/test sets from whole datasets
set.seed(200)
Ftest_samples = tests
Ftest_geno = geno1[Ftest_samples,]
Ftest_pheno = phenotype[Ftest_samples]
Train_samples = setdiff(names(phenotype), Ftest_samples)
geno2 <- geno1[Train_samples,]
phenotype1 <- phenotype[Train_samples]
names(phenotype1) <- rownames(pheno1)[Train_samples]
cv_samples <- sample(1:cv, nrow(geno2), replace = TRUE) 
cv_val_mean = mean(table(cv_samples))
cv_tra_mean = 0
for ( i in 1:cv){
    cv_mean = sum(table(cv_samples[which(cv_samples != i)]))
    cv_tra_mean = c(cv_tra_mean, cv_mean)
}
cv_tra_mean = sum(cv_tra_mean) / cv

#(5) Making output directory and Saving temporary files
CVN = paste("CV", cv, sep="")
INI = paste("Ini", INIbk, sep="")
SEL = paste("Sel", sel_snps, sep="")
TEN = paste("Tes", length(Ftest_samples), sep="")
dir.create("Results")
setwd("Results")
if (is.null(gwasfile)){
  w_dir <- paste0(colnames(pheno1)[2],"_", mmm, "_", PC, "_", CVN, "_", INI, "_", SEL, "_", TEN, "_wo_gwas/")
} else {
  w_dir <- paste0(colnames(pheno1)[2],"_", mmm, "_", PC, "_", CVN, "_", INI, "_", SEL, "_", TEN, "_with_gwas/")
}
if(!dir.exists(w_dir)){dir.create(w_dir)}
setwd(w_dir)
saveRDS(Ftest_geno, file = "Final_test_genotype.rds", compress = "gzip")
saveRDS(Ftest_pheno, file = "Final_test_phenotype.rds", compress = "gzip")
rm(list=c("Ftest_geno", "Ftest_pheno", "geno1"))
invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))


#(6) Staring Analysis.
N <- ncol(geno2)
cat(paste0("Train set: Analysis start with ", nrow(geno2), " samples and ", ncol(geno2), " SNPs", '\n\n'), sep="")
init_selsnp = sel_snps
all_train_acc = NULL
all_val_acc = NULL
selected_train_acc = NULL
selected_val_acc = NULL
CV_results = list()
nSamGenSum = list()
tsum = list()
inisum = list()
if (mmm == "RRB_BTS"){
  mmm1 = c("RRB", "BTS")
} else {
  mmm1 = mmm
}
t_time = Sys.time()
for (mm in mmm1){
   for (j in 1:cv){
     ini_snps = ini_snps_bk
     o_delta = 1
     sel_snps = init_selsnp
     cat(paste0("CV Number ", j, ": Start","\n"), sep="")
     f_time = Sys.time()

     #(6-1)Cheking prediction (correlation) rate for all markers
     cat(paste0("CV Number ", j, ": Check prediction (correlation) rate for all markers using ", mm, " model","\n"), sep="")
     train_samples = which(cv_samples != j)
     train_geno = geno2[train_samples,]
     val_geno = geno2[-train_samples,]
     train_pheno = phenotype1[train_samples] 
     val_pheno = phenotype1[-train_samples]
  
     nTrainSam = length(train_samples)
     nValSam = length(val_pheno)
     nGeno = dim(train_geno)[2]
     nSamGen = c(nTrainSam, nValSam, nGeno)
     nSamGenSum[[paste0("CV-", j, "_", mm)]] = nSamGen
  
     if (allm == TRUE){
       if (mm == "RRB"){
         fit = RRb_func(train_pheno, train_geno, val_geno)
         val_pred = fit$val_predicted
         train_pred = fit$train_predicted
         fit_model = fit$model
       }
       if (mm == "BTS"){
         fit = RF_func(train_pheno, train_geno, val_geno, "MS")
         val_pred = fit$val_predicted
         train_pred = fit$train_predicted
         fit_model = fit$model
       }

       train_acc = cor(train_pred, train_pheno, use = "complete")
       val_acc = cor(val_pred, val_pheno, use = "complete")
       all_train_acc = c(all_train_acc, as.vector(train_acc)) 
       all_val_acc = c(all_val_acc, as.vector(val_acc))
     } else {
       train_acc = "-"
       val_acc = "-"
       all_train_acc = c(all_train_acc, as.vector(train_acc))
       all_val_acc = c(all_val_acc, as.vector(val_acc))
     }
     cat(paste0("CV Number ", j, ": Prediction (correlation) rate for all markers using ", mm, " is ", train_acc, " train acc ", val_acc, " val acc ", "\n"), sep="")
     cat(paste0("CV Number ", j, ": Prediction running time for all markers is ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"), sep="")
     inisum[[paste0("CV-", j, "_", mm)]] <- format(difftime(Sys.time(), f_time), usetz = TRUE)

     #(6-2)Selecting initial markers 
     ix1 = ix
     if(length(preset_fname) != 0){
       preset_idx = match(preset_fname, ix1)
       ix1 = ix1[-preset_idx]
       pp_ix = preset_fname
       if (ini_snps > 0) {
         if (mm == "BTS"){
           init_set = ix1[1:ini_snps]
           train_ix = c(pp_ix, init_set)
         } else {
           init_set_t = ix1[1:ini_snps]
           init_set = c(pp_ix, init_set_t)
           zv = apply(train_geno[,init_set], 2, function(x) length(unique(x)) == 1)
           if (sum(cv) == 1){
              highCorr = NULL
           } else {
             cor_mat = cor(train_geno[,init_set][, !zv])
             highCorr = findCorrelation(cor_mat, 0.9)  
           }
           if (length(highCorr) == 0){
             train_ix = init_set
           } else {
             init_set_tx = match(init_set_t, init_set)
             init_set_tx_int = intersect(init_set_tx, highCorr)
             if (length(init_set_tx_int) > 0 ){
               train_ix = init_set[-init_set_tx_int]
             } else {
               train_ix = init_set
             }
           }
         }
       } else {
         train_ix = pp_ix
       }
       rm(pp_ix)  
     } else {
       if (mm == "BTS"){
         init_set = ix1[1:ini_snps]
         train_ix = init_set
       } else {
         init_set = ix1[1:ini_snps]
         zv = apply(train_geno[,init_set], 2, function(x) length(unique(x)) == 1)
         if (sum(cv) == 1){
           highCorr = NULL
         } else {
           cor_mat = cor(train_geno[,init_set][, !zv])
           highCorr = findCorrelation(cor_mat, 0.9)  
         }
         if (length(highCorr) == 0){
           train_ix = ix1[1:ini_snps]
         } else {
           train_ix = init_set[-highCorr]
         }
       }
     }
  
     r=0
     except_ix = NULL
     except_tix = NULL
     except_ixc = 0
     add_ix = train_ix
     result_train_acc = 0        
     result_val_acc = 0
     ms_out = list()
     t_add = train_ix
     train_itr_geno = train_geno[,train_ix]
     val_itr_geno = val_geno[,train_ix]
     invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))
     #(6-3)Cheking prediction (correlation) rate for initial markers
     if (mm == "RRB"){
       fit_itr = RRb_func(train_pheno, train_itr_geno, val_itr_geno)
       val_pred_itr_pheno = fit_itr$val_predicted
       train_pred_itr_pheno = fit_itr$train_predicted
       fit_itr_model = fit_itr$model
     }
     if (mm == "BTS"){
       fit_itr = RF_func(train_pheno, train_itr_geno, val_itr_geno, "MS")
       val_pred_itr_pheno = fit_itr$val_predicted
       train_pred_itr_pheno = fit_itr$train_predicted
       fit_itr_model = fit_itr$model
     }

     train_itr_acc = cor(train_pred_itr_pheno, train_pheno, use = "complete")
     val_itr_acc = cor(val_pred_itr_pheno, val_pheno, use = "complete")
  
  #(6-4)Checking initial result
     ms_out[[1]] = c(train_itr_acc, train_ix)
     k=2 
     cat(paste0("CV Number ", j, ": Initial marker set's correlation rate is ", round(train_itr_acc,6), " in train, ", round(val_itr_acc,6), " in validation ", "\n"), sep="")
     result_train_acc = c(result_train_acc, train_itr_acc)
     result_val_acc = c(result_val_acc, val_itr_acc)
  
  #(6-5)Starting marker selection
     i=1
     max_noid = 0
     max_noid2 = 0
     s_time = Sys.time()
     cat(paste0("CV Number ", j, ": Start marker selection at ", format(s_time, usetz = T), "\n"), sep="")
     readc = length(train_ix)
     readcc = length(train_ix)
     break_n = round(length(ix1) * 0.2)
     cat(paste0(" -- Stop condition: if the correlation rate does not increase ", break_n, " times consecutively (20% of the total number of markers), the selection process is stopped.", "\n"), sep="")
     while (val_itr_acc < acc1){
       readcc = readcc + 1
       if (max_noid > 0 & max_noid == break_n){
         cat(paste0("STOP condition1: correlation rate did not increase ", break_n, " times in succession, so stop marker selection", "\n\n"), sep="")
         break
       }
       c5 = i %% 500
       if (c5 == c(0)){
         cat(paste0(i, " markers were read...", "\n"), sep="")
       }
    #(6-6) Selecting SNPs
       if (mm == "BTS"){
         add_snps = ix1[(ini_snps + (sel_snps*(i-1)+1)):(ini_snps + (sel_snps*(i)))]
         new_train_ix <- c(train_ix, add_snps)
         if (sum(is.na(new_train_ix)) >= 1){
           cat("STOP condition2: All markeres were read, so stop this CV", '\n\n')
           break
         }
       } else {
         add_snps_t = ix1[(ini_snps + (sel_snps*(i-1) +1)):(ini_snps + (sel_snps*(i)))]
         init_set <- c(train_ix, add_snps_t)
         if (sum(is.na(init_set)) >= 1){
           cat("STOP condition2: All markeres were read, so stop this CV", '\n\n')
           break
         }
         add_snps_t_ix = match(add_snps_t, init_set)
         zv <- apply(train_geno[,init_set], 2, function(x) length(unique(x)) == 1)
         cor_mat <- cor(train_geno[,init_set][, !zv])
         highCorr <- findCorrelation(cor_mat, 0.9)
         if (length(highCorr) == 0){
           new_train_ix <- init_set
           add_snps <- add_snps_t
         } else {
           if (sum(highCorr %in% add_snps_t_ix) >=1) {
             add_snps_t_ix2 = add_snps_t_ix[-match(highCorr, add_snps_t_ix)[!is.na(match(highCorr, add_snps_t_ix))]]
             new_train_ix = c(train_ix, init_set[add_snps_t_ix2])
             add_snps= init_set[add_snps_t_ix2]
           } else {
             new_train_ix <- init_set
             add_snps <- add_snps_t
           }
         }
       }
    
    #(6-7) Cheking prediction (correlation) rate for selected markers 
       train_itr_geno <- train_geno[,new_train_ix]
       val_itr_geno <- val_geno[,new_train_ix]
       if (mm == "RRB"){
         fit_itr = RRb_func(train_pheno, train_itr_geno, val_itr_geno)
         train_pred_itr_pheno = fit_itr$train_predicted
         val_pred_itr_pheno = fit_itr$val_predicted
         fit_itr_model = fit_itr$model
       }
       if (mm == "BTS"){
         fit_itr = RF_func(train_pheno, train_itr_geno, val_itr_geno, "MS")
         train_pred_itr_pheno = fit_itr$train_predicted
         val_pred_itr_pheno = fit_itr$val_predicted
         fit_itr_model = fit_itr$model
       }

       acc_train = cor(train_pred_itr_pheno, train_pheno, use = "complete")
       acc_val = cor(val_pred_itr_pheno, val_pheno, use = "complete")

       check = ((acc_train - train_itr_acc >= delta) & (acc_val - val_itr_acc >= delta))

       if (check == FALSE) {
         except_ix = c(except_ix, add_snps)
         max_noid = max_noid + 1
         max_noid2 = max_noid2 + 1
       } else {
         readc = readc + 1
         ms_out[[k]] = c(acc_train, add_snps)
         train_ix = c(train_ix, add_snps)
         t_add = c(t_add, add_snps)
         train_itr_acc = acc_train
         val_itr_acc = acc_val
         add_ix = c(add_ix, add_snps)
         result_train_acc = c(result_train_acc, train_itr_acc)
         result_val_acc = c(result_val_acc, val_itr_acc)
         t_delta = result_train_acc[length(result_train_acc)] - result_train_acc[length(result_train_acc)-1]
         v_delta = result_val_acc[length(result_val_acc)] - result_val_acc[length(result_val_acc)-1] 
         cat(paste0("Add: ", round(train_itr_acc,6)," (delta: ", round(t_delta, 6), ") train corr., ", round(val_itr_acc,6)," (delta: ", round(v_delta, 6), ") val corr., ", readc, " markers (", readcc, " reading), time: ", format(difftime(Sys.time(), s_time), usetz = TRUE), "\n", sep=""))
         max_noid = 0
         max_noid2 = 0
         r=0
         sel_snps = init_selsnp
       }
    
    #(6-8) Setting nex iteration
       if (mm == "BTS"){
         if (sel_snps != 1){
           if ((length(add_ix) + length(except_ix)) >= (sel_snps * 15)){
             except_ixc = except_ixc + 1
             except_tix = c(except_tix, except_ix)
             if (r>=2 & max_noid2 >=28  ){
               except_tix = c(except_tix, except_ix)
               t_tix = c(except_tix, add_ix)
               search_ix = unique(t_tix)
               ix1 <- ix1[!(ix1 %in% search_ix)]
               max_noid2 = 0
               except_ix <- NULL
               add_ix <- NULL
               ini_snps = 0
               i <- 1
               r=0
               except_ixc <- 0
               except_tix <- NULL
             }
             if (sum(table(except_tix) == 1) >= (sel_snps * 15)) {
               t_tix = c(except_tix, add_ix)
               search_ix = unique(t_tix)
               ix1 <- ix1[!(ix1 %in% search_ix)]
               except_ix <- NULL
               add_ix <- NULL
               ini_snps = 0
               i <- 1
               except_ixc <- 0
               except_tix <- NULL
               r= r+1
             }
             if (except_ixc == 1 ){
               search_ix = c(add_ix, except_tix)
               ix1_t1 = setdiff(ix1, search_ix)
               set.seed(200)
               except_ix2 = sample(except_tix)
               ix1 = c(except_ix2, ix1_t1)
               except_ix <- NULL
               add_ix <- NULL
               ini_snps = 0
               i <- 1
               r= r+1
             } else if (except_ixc > 1) {
               before_except_add = intersect(except_tix, add_ix)
               except_tix = except_tix[!(except_tix %in% before_except_add)]
               except_tix_table = table(except_tix)
               except_once = names(except_tix_table[except_tix_table <= 1])
               except_more = names(except_tix_table[except_tix_table > 1])
               except_tix = except_once
               search_ix = c(add_ix, except_more, except_once)
               ix1_t1 = setdiff(ix1, search_ix)
               set.seed(200)
               except_ix2 = sample(except_once)
               ix1 = c(except_ix2, ix1_t1)
               except_ix <- NULL
               add_ix <- NULL
               ini_snps = 0
               i <- 1
               r=r+1
             } else {
               a=0
             }
           } else {
             i <- i + 1
           }
         } else {
           i <- i + 1
         }
       }
    
       if (mm == "RRB" ) {
         if (i == length(ix1)/sel_snps){
           ix_sum <- c(add_ix, except_ix)
           ix1 <- ix1[-which(ix1 %in% ix_sum)]
           add_ix <- NULL
           except_ix <- NULL
           ini_snps = 0
           i <- 1    
         } else {
           i <- i+1
         }
       }
       k = k+1
       invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))
     }
   
  
  #(6-9) Saving and iteration result
     selected_train_acc = c(selected_train_acc, train_itr_acc)
     selected_val_acc = c(selected_val_acc, val_itr_acc)
     invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))
     cat(paste0("CV Number ", j, ": # of selected markers is ", length(train_ix), ", and correlation rates for train and validation are ", round(train_itr_acc,6), ', ', round(val_itr_acc,6), "\n"), sep="")
     cat(paste0("CV Number ", j, ": Running time for marker selection is ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"),sep="")
     cat(paste0("CV Number ", j, ": End ", mm, " model."  ,"\n\n\n"), sep="")
     tsum[[paste0("CV-", j, "_", mm)]] = format(difftime(Sys.time(), f_time), usetz = TRUE)

     CV_name = paste0("CV-", j, "_", mm)
     CV_results[[CV_name]] = train_ix
     
     a <- plyr::ldply(ms_out, rbind)   
     count <- length(which(!is.na(as.vector(a[1,-1])))) 
     index <- count
     for (i in 2:nrow(a)){  
       count <- count + length(which(!is.na(as.vector(a[i,-1]))))
       index <- c(index, count)
     }
     a <- data.frame(Index = index, a)   
     colnames(a)[c(2,3)] <- c("Corr", "Marker")
     a$Corr <- as.numeric(as.vector(a$Corr))
     marker_save = paste0("CV-", j, "_", mm, ".MarkerList.log", sep="")
     write.table(a, file = marker_save, row.names = F, sep='\t', quote=F)
  }
}

#(7) Load and save all result files
cat(paste0("END: Total running time is ", format(difftime(Sys.time(), t_time), usetz = TRUE), "\n"), sep="")
CV_results = CV_results[order(names(CV_results), decreasing=FALSE)]
nSamGenSum = nSamGenSum[order(names(nSamGenSum), decreasing=FALSE)]
inisum = inisum[order(names(inisum), decreasing=FALSE)]
tsum = tsum[order(names(tsum), decreasing=FALSE)]
marker_vec <- unlist(CV_results)
marker_table <- table(marker_vec)[order(table(marker_vec), decreasing = TRUE)]
marker_select <- names(marker_table)
msummary_list <- list()
name_tmp = NULL
for (m1 in mmm1){
  tt = paste("CV-", seq(cv), "_", m1, sep="")
  name_tmp = c(name_tmp, tt)
}
for (iter in names(CV_results)){
  msummary_list[[iter]] = sapply(marker_select, function(x) as.integer(x %in% CV_results[[iter]]))
}
marker_summary = data.frame(msummary_list, check.names = FALSE)
marker_summary = marker_summary[,name_tmp]
marker_summary = data.frame("Num"=seq(1,dim(marker_summary)[1]), "Marker"=rownames(marker_summary), marker_summary, "Total_score"=apply(marker_summary,1,sum), check.names = FALSE)
write.table(marker_summary, file=paste("CV_", mmm, "_Marker_scores.txt", sep=""), quote=FALSE, sep="\t", row.names = FALSE)
write(marker_select, file = paste("CV_Marker_lists.txt", sep=""))

if (allm == TRUE){
  m_all_train_acc = round(mean(all_train_acc),6)
  sd_all_train_acc = round(sd(all_train_acc),6)
  m_all_val_acc = round(mean(all_val_acc),6)
  sd_all_val_acc = round(sd(all_val_acc),6)
  add_row = c("Total", mmm, cv_tra_mean, cv_val_mean, dim(geno2)[2], paste(round(mean(all_train_acc), 6), " (",round(sd(all_train_acc),6),")", sep=""), paste(round(mean(all_val_acc), 6), " (",round(sd(all_val_acc),6),")", sep=""), paste(round(mean(selected_train_acc), 6), " (",round(sd(selected_train_acc),6),")", sep=""), paste(round(mean(selected_val_acc), 6), " (",round(sd(selected_val_acc),6),")", sep=""), length(preset_fname), paste(length(unique(marker_vec)), " (", sum(table(marker_vec) == (cv)), ")", sep=""), "-", format(difftime(Sys.time(), t_time), usetz = TRUE)) 
  cat(paste0("END: Train correlation rate of whole markers is ", m_all_train_acc, " ", paste("\U00B1"), " ", sd_all_train_acc, "\n"), sep="")
  cat(paste0("END: Validation correlation rate of whole markers is ", m_all_val_acc, " ", paste("\U00B1"), " ", sd_all_val_acc, "\n"), sep="")
} else {
  add_row = c("Total", mmm, cv_tra_mean, cv_val_mean, dim(geno2)[2], "-", "-", paste(round(mean(selected_train_acc), 6), " (",round(sd(selected_train_acc),6),")", sep=""), paste(round(mean(selected_val_acc), 6), " (",round(sd(selected_val_acc),6),")", sep=""), length(preset_fname), paste(length(unique(marker_vec)), " (", sum(table(marker_vec) == (cv)), ")", sep=""), "-", format(difftime(Sys.time(), t_time), usetz = TRUE)) 
}
cat(paste0("END: Train correlation rate of selected markers is ", round(mean(selected_train_acc),6), " ", paste("\U00B1"), " ", round(sd(selected_train_acc),6), "\n"),sep="")
cat(paste0("END: Validation correlation rate of selected markers is ", round(mean(selected_val_acc),6), " ", paste("\U00B1"), " ", round(sd(selected_val_acc),6), "\n"),sep="")

names(all_train_acc) = name_tmp
names(selected_train_acc) = name_tmp
names(all_val_acc) = name_tmp
names(selected_val_acc) = name_tmp

tsummary = data.frame("Train_Corr_for_all_markers"= all_train_acc[names(CV_results)], "Val_Corr_for_all_markers"= all_val_acc[names(CV_results)], "Train_Corr_for_selected_markers"=selected_train_acc[names(CV_results)], "Val_Corr_for_selected_markers"=selected_val_acc[names(CV_results)])
tsummary = data.frame("CV"=unlist(lapply(strsplit(names(CV_results), "_"), function(x) x[1])), "Model"=unlist(lapply(strsplit(names(CV_results), "_"), function(x) x[2])), t(data.frame(nSamGenSum, check.names = FALSE))[names(CV_results),], tsummary, "Pre-Selected_markers"=length(preset_fname), "Total-Selected_markers"=unlist(lapply(CV_results, length)), t(data.frame(inisum, check.names = FALSE))[names(CV_results),], t(data.frame(tsum, check.names = FALSE))[names(CV_results),])
colnames(tsummary) = c("CV", "Model", "Train_Samples", "Val_Samples", "Used_Markers", "Train_Corr_for_all_markers", "Val_Corr_for_all_markers", "Train_Corr_for_selected_markers", "Val_Corr_for_selected_markers",  "Pre_Selected_markers", "Total_Selected_markers", "Initial_runtime", "Total_runtime")
tsummary = tsummary[order(tsummary$Model),]
tsummary = rbind(as.matrix(tsummary), add_row)
write.table(tsummary, file=paste("CV_", mmm, "_Selection_summary.txt", sep=""), quote=FALSE, sep="\t", row.names = FALSE)

FtrainGeno = geno2[,marker_select]
FtrainPheno = phenotype1
saveRDS(FtrainGeno, file = "Final_train_genotype.rds", compress = "gzip")
saveRDS(FtrainPheno, file = "Final_train_phenotype.rds", compress = "gzip")
cat("\n\n")

#############################################################################

