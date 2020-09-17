##########
#GMStool v1 for Marker selecion, made by Seongmun Jeong and JaeYoon Kim in 2020-09-17.
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
parser$add_argument("-c", "--corr", type = "double", default = 1.0, help = "Target correlation rate for the validation set. If the correlation rate for the validation set reaches the target correlation rate, the marker selection is stopped and the results up to that point is returned (Defualt 1.0)")
parser$add_argument("-d", "--increment", type = "double", default = 0.00005, help = "Increament of correlation rate. In order for a marker to be selected, the correlation rate of the marker must be greater than or equal to the value of the correlation rate of the immediately previously selected marker plus this increment value. This condition is applied to the train validation set, so the marker is selected only when satisfying this condition for the validation set (Default 0.00005)")
parser$add_argument("-x", "--stopx", type = "double", default = 20, help = "Stop condition: this condition uses a number corresponding to x percent of the total number of input markers. If the correlation rates of validation sets do not improved as many as this number in a row, the marker selection for the corresponding CV is stopped and the results up to that point are returned (Default 20)")
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
breakx <- as.numeric(args$stopx) / 100
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
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(caret))
library(ggplot2)
library(ggpmisc)
library(data.table)
library(tidyverse)
library(caret)
library(foreach)
library(iterators)
library(doParallel)
options(warn=-1)
source("./src/GMS_main.R")
source("./src/load_files.R")
source("./src/RRb_func.R")
source("./src/RF_func.R")
##(1-2) Set gpu/cpu
if (mmm == "RRB" | mmm == "RRB_BTS") {
  library(rrBLUP)
  if (gpu_use) {
    library(gpuR)
    source("./gpu.mixed.solve.R")
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


#(4) Preparing train/test/validation sets from whole datasets
Ftest_samples = tests
Ftest_geno = geno1[Ftest_samples,]
Ftest_pheno = phenotype[Ftest_samples]
Train_samples = setdiff(names(phenotype), Ftest_samples)
geno2 <- geno1[Train_samples,]
phenotype1 <- phenotype[Train_samples]
names(phenotype1) <- rownames(pheno1)[Train_samples]
cv_samples = rep(1:cv, ceiling(nrow(geno2)/cv))
cv_samples = cv_samples[1:nrow(geno2)]
set.seed(200)
cv_samples = sample(cv_samples)
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
  w_dir <- paste0(colnames(pheno1)[2],"_", mmm, "_", PC, "_", CVN, "_", INI, "_", SEL, "_", TEN, "_MT_wo_gwas/")
} else {
  w_dir <- paste0(colnames(pheno1)[2],"_", mmm, "_", PC, "_", CVN, "_", INI, "_", SEL, "_", TEN, "_MT_with_gwas/")
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
init_selsnp <- sel_snps
if (mmm == "RRB_BTS"){
  mmm1 = c("RRB", "BTS")
} else {
  mmm1 = mmm
}
t_time <- Sys.time()

registerDoParallel(detectCores()/(cv))

all_train_acc = NULL
all_val_acc = NULL
selected_train_acc = NULL
selected_val_acc = NULL
CV_results <- list()
nSamGenSum <- list()
tsum <- list()
inisum <- list()

for (mm in mmm1){
  cat(paste0("Marker selection: ", mm, " (multithreading) ...", "\n"), sep="")
	results <- foreach(j=1:cv) %dopar%
	  GMS_main(ini_snps_bk, init_selsnp, j, mm, cv_samples, geno2, phenotype1, preset_fname, ix1, breakx)
	assign(paste0("results_", mm), value = results)
}


for (mm in mmm1){
	results = get(x = paste0("results_", mm))
	for (k in 1:cv){
		all_train_acc <- c(all_train_acc, results[[k]]$all_train_acc)
                all_val_acc <- c(all_val_acc, results[[k]]$all_val_acc)
		selected_train_acc <- c(selected_train_acc, results[[k]]$selected_train_acc)
                selected_val_acc <- c(selected_val_acc, results[[k]]$selected_val_acc)
		CV_results[[paste0("CV-", k, "_", mm)]] <- results[[k]]$CV_results
		nSamGenSum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$nSamGenSum
		tsum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$tsum
		inisum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$inisum
	}
}


#(7) Load and save all result files
cat("\n")
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

