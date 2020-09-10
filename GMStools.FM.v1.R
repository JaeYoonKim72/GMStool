##########
#GMStool v1 for Final modeling, made by Jae-Yoon Kim and Seongmun Jeong in 2020-09-09.
##########

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-m", "--model", default = "RRB", help = "Prediction model (RRB_RF_DNN_CNN)")
parser$add_argument("-d", "--dir", default = NULL, help = "Result directory of marker selection")
parser$add_argument("-gw", "--gwas", default = NULL, help = "GWAS result file. If GWAS file is not provided, marker information file should be provided.")
parser$add_argument("-i", "--info", default = NULL, help = "Marker information file.  Required if GWAS file is not provided")
parser$add_argument("-pe", "--permutation", default = 50, help = "The number of permutations per each modeling")
parser$add_argument("-gpu", "--GPU_usage", default = FALSE, help = "If TRUE, DNN and CNN are calculated using GPU")
parser$add_argument("-t", "--time", type = "double", default = 1, help = "Runtime cut-off for permutatios of each modeling")

#(0) Prepareing analysis
prd = getwd()
args <- parser$parse_args()
dr <- paste(prd, args$dir, sep="/")
mm <- args$model
gwasfile <- paste(prd, args$gwas, sep="/")
infofile <- paste(prd, args$info, sep="/")
gpu_use <- as.logical(args$GPU_usage)
time_cutoff <- args$time
per <- as.integer(args$permutation)

if (is.null(args$info)){
    infofile1 = "not provided"
} else {
    infofile1 = infofile
}
if (is.null(args$gwas)){
    gwasfile1 = "not provided"
} else {
    gwasfile1 = gwasfile
}
if (is.null(args$info) & is.null(args$gwas)){
  cat("Error1: provide gwas result file or marker information file", '\n')
  quit(save = "no")
}
mms = unlist(strsplit(mm, "_"))
for (ms in mms){
    if (ms %in% c("RRB", "RF", "DNN", "CNN")==FALSE){
    cat("Error2: invalide selection method", '\n')
    quit(save = "no")
    }
}

cat("\n\n")
cat("Reading arguments ==================================================", "\n")
cat("1)  Directory which markers selected is", dr, "\n")
cat("2)  Prediction model is", mm, "\n")
cat("3)  GWAS result file is", gwasfile1, "\n")
cat("4)  Marker information file is", infofile1, "\n")
cat("5)  Usage of GPU is", gpu_use, "\n")
cat("6)  # of permutation is", per, "\n")
cat("7)  Limit of permutation runtime is", time_cutoff, "hours", "\n")
cat("================================================================", "\n\n")


#(1) Load assicated packages, and set gpu/cpu and threads.
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(ggplot2)
library(ggpmisc)
library(data.table)
library(tidyverse)
library(caret)
options(warn=-1)
source("./src/RRb_func.R")
source("./src/RF_func.R")
source("./src/DNN_func.R")
source("./src/CNN_func.R")
setwd(dr)

##(1-2) Set gpu/cpu
if ("DNN" %in% mms | "CNN" %in% mms) {
  suppressPackageStartupMessages(library(tensorflow))
  suppressPackageStartupMessages(library(keras))
  library(tensorflow)
  library(keras)
  tensorflow::tf$random$set_seed(200)
  tf$compat$v1$set_random_seed(200)
  if (gpu_use == FALSE) {
    cpus <- tf$config$experimental$list_physical_devices("CPU")
    tf$config$experimental$set_visible_devices(cpus[1], "CPU")
  }
  metric_cor <- custom_metric("metric_cor", function(y_true, y_pred){
    fsp = y_pred - k_mean(y_pred)
    fst = y_true - k_mean(y_true)
    devP = k_std(y_pred)
    devT = k_std(y_true)
    return(k_mean(fsp*fst)/(devP*devT))
  })
  metric_r2_score <- custom_metric("r2_score", function(y_true, y_pred){
    fsp = y_pred - k_mean(y_pred)
    fst = y_true - k_mean(y_true)
    devP = k_std(y_pred)
    devT = k_std(y_true)
    return(k_square((k_mean(fsp*fst)/(devP*devT))))
  })
} 
if ("RRB" %in% mms | "RF" %in% mms) {
  library(rrBLUP)
  suppressPackageStartupMessages(library(randomForest))
  library(randomForest)
}


#(2) Prepare data set
cat(paste0("Prepareing: Reading test and train datasets.", "\n\n"), sep="")
FtrainGenoOTH = readRDS("Final_train_genotype.rds")
FtrainPheno = readRDS("Final_train_phenotype.rds")
FtestGenoOTH = readRDS("Final_test_genotype.rds")
FtestPheno = readRDS("Final_test_phenotype.rds")
Markerlist = as.character(read.table("CV_Marker_lists.txt", header=F)[,1])
FtrainGenoOTH = FtrainGenoOTH[,Markerlist]
FtestGenoOTH = FtestGenoOTH[,Markerlist]
cat(paste0("Marker: ", length(Markerlist), " input markers", "\n\n"), sep="")


#(3) Physical positional sort for CNN modeling
if ("CNN" %in% mms){
    cat("\n\n")
    cat("CNN: input markers are sorted in ascending order according to chromosome number and physical position", '\n')
    train_ms = colnames(FtrainGenoOTH)
    test_ms = colnames(FtestGenoOTH)
    if (!is.null(args$gwas)){
        cat(" - using GWAS result file", "\n\n")
        gwas_results <- data.frame(fread(gwasfile, header=T, check.names=FALSE, data.table=FALSE), row.names=1)
        train_ms_info = gwas_results[train_ms,]
        train_ms_sort_df = train_ms_info[order(train_ms_info[,1], train_ms_info[,2]),]
        train_ms_sort = rownames(train_ms_sort_df)
        FtrainGenoCNN = FtrainGenoOTH[,train_ms_sort]
        FtestGenoCNN = FtestGenoOTH[,train_ms_sort]
    }
    if (!is.null(args$info)){
        cat(" - using marker-information file", "\n\n")
        infotable = data.frame(fread(infofile, header=T, check.names=FALSE, data.table=FALSE), row.names=1)
        train_ms_info = infotable[train_ms,]
        train_ms_sort_df = train_ms_info[order(train_ms_info[,1], train_ms_info[,2]),]
        train_ms_sort = rownames(train_ms_sort_df)
        FtrainGenoCNN = FtrainGenoOTH[,train_ms_sort]
        FtestGenoCNN = FtestGenoOTH[,train_ms_sort]
    }
    if (!is.null(args$gwas) & !is.null(args$info)){
        cat(" - Warning: since neither GWAS file or Marker-information file was defined by the user, the input markers are sorted in the most selected order in all CVs (original sorting state).", "\n\n")
        FtrainGenoCNN = FtrainGenoOTH
        FtestGenoCNN = FtestGenoOTH
    }
}


#(4) Final Modelling
Total_list = list()
for (mm in mms){
  FinalMarker = Markerlist
  Total_list[["mm"]] = list("Model"=NULL, "TrainCorr"=NULL, "TestCorr"=NULL, "Bat"=NULL)
  FmodelList = list()
  FcorrList = list("Corr"=NULL, "Model"=NULL, "Bat"=NULL, "Predict"=NULL)
  if (mm == "CNN" | mm == "DNN" | mm == "RF"){
    cat(paste0("Start ", mm, " permutation. ", "\n"), sep="")
    if (mm == "CNN"){
        FtrainGeno = FtrainGenoCNN
        FtestGeno = FtestGenoCNN
    } else {
        FtrainGeno = FtrainGenoOTH
        FtestGeno = FtestGenoOTH
    }

    f_time <- Sys.time()
    for (ite in seq(1,per)){
      if (difftime(Sys.time(), f_time, units = "hours") > time_cutoff){
        cat("STOP: Computational time is over, so stop permutation", '\n\n')
        break
      }
      set.seed(ite)
      train_samples <- sample(1:nrow(FtrainGeno), size = nrow(FtrainGeno)*0.8, replace = FALSE)
      train_geno <- FtrainGeno[train_samples,]
      val_geno <- FtrainGeno[-train_samples,]
      train_pheno <- FtrainPheno[train_samples]  
      val_pheno <- FtrainPheno[-train_samples]
      cat(paste0(" - ", mm, ": permutation ", ite," ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"), sep="")
      if (mm == "DNN"){
        tensorflow::tf$random$set_seed(200)
        fit = DNN_func(train_pheno, train_geno, val_geno)  
        val_pred = fit$val_predicted
        train_pred = fit$train_predicted
        train_corr = cor(train_pheno, train_pred)
        val_corr = cor(val_pheno, val_pred)
        fit_model = fit$model
        fit_batch = fit$batch
      }
      if (mm == "CNN"){
        tensorflow::tf$random$set_seed(200)
        fit = CNN_func(train_pheno, train_geno, val_geno)
        val_pred = fit$val_predicted
        train_pred = fit$train_predicted
        train_corr = cor(train_pheno, train_pred)
        val_corr = cor(val_pheno, val_pred)
        fit_model = fit$model
        fit_batch = fit$batch
      }
      if (mm == "RF"){
        fit = RF_func(train_pheno, train_geno, val_geno, "FM")
        val_pred = fit$val_predicted
        train_pred = fit$train_predicted
        train_corr = cor(train_pheno, train_pred)
        val_corr = cor(val_pheno, val_pred)
        fit_model = fit$model
        fit_batch = 1
        FmodelList[[ite]] = fit_model
      }
      FcorrList$TCorr = c(FcorrList$TCorr, train_corr)
      FcorrList$VCorr = c(FcorrList$VCorr, val_corr)
      FcorrList$Model = c(FcorrList$Model, fit_model)
      FcorrList$Bat = c(FcorrList$Bat, fit_batch)
    }
    MedCorrMADT = mad(FcorrList$TCorr, na.rm=TRUE)
    MedCorrMADV = mad(FcorrList$VCorr, na.rm=TRUE)
    FcorrList$TCorr[is.na(FcorrList$TCorr)] = 0
    FcorrList$VCorr[is.na(FcorrList$VCorr)] = 0
    MaxCorrV = max(FcorrList$VCorr)
    MaxPhenoV = FcorrList$Pred[which(FcorrList$VCorr == MaxCorrV)]
    MaxCorrT = FcorrList$TCorr[which(FcorrList$VCorr == MaxCorrV)]

    if (mm == "CNN" | mm == "DNN"){
      MaxModel = FcorrList$Model[which(FcorrList$VCorr == MaxCorrV)]
      MaxBatch = FcorrList$Bat[which(FcorrList$VCorr == MaxCorrV)]
      if (mm == "CNN"){
        mind = dim(FtestGeno)[1]
        mnum = dim(FtestGeno)[2]
        FtestGeno2 <- FtestGeno + 1
        dim(FtestGeno2) <- c(mind, 1, mnum)
        FtestGeno2 <- tf$cast(FtestGeno2, tf$float32) /3
        FtestGeno2 <- tf$expand_dims(FtestGeno2, axis=-1L)
        Fpredicted_pheno = predict(MaxModel[[1]], FtestGeno2, batch_size=MaxBatch)
        Fcorr = cor(Fpredicted_pheno, FtestPheno)
      } else {
        FtestGeno2 = FtestGeno + 1
        FtestGeno2 <- tf$cast(FtestGeno2, tf$float32) /3
        Fpredicted_pheno = predict(MaxModel[[1]], FtestGeno2, batch_size=MaxBatch)
        Fcorr = cor(Fpredicted_pheno, FtestPheno)
      }
    } else{
      MaxModel = FmodelList[[which(FcorrList$VCorr == MaxCorrV)]]
      FtestGeno2 = (FtestGeno + 1) / 3
      MaxBatch = 1
      Fpredicted_pheno = predict(MaxModel, FtestGeno2)
      Fcorr = cor(Fpredicted_pheno, FtestPheno)
    }
    
    Total_list[[mm]]$Model = MaxModel
    Total_list[[mm]]$TrainCorr = MaxCorrT
    Total_list[[mm]]$ValCorr = MaxCorrV
    Total_list[[mm]]$Bat = MaxBatch
    Total_list[[mm]]$Predict = Fpredicted_pheno
    Total_list[[mm]]$TestCorr = Fcorr
    Total_list[[mm]]$TrainCorrlist = FcorrList$TCorr
    Total_list[[mm]]$ValCorrlist = FcorrList$VCorr
    Total_list[[mm]]$MedCorrMADV = MedCorrMADV
    Total_list[[mm]]$MedCorrMADT = MedCorrMADT
    
    cat(paste0(" - ", mm, ": Final correlations for training and validation dataset are ", round(MaxCorrT,6), ", ", round(MaxCorrV,6), "\n"), sep="")
    cat(paste0(" - ", mm, ": Median absoulte deviations for ", per, " replications are ", round(MedCorrMADT, 6), ", ", round(MedCorrMADV, 6), "\n"), sep="")
    cat(paste0(" - ", mm, ": Final correlation for test dataset is ", round(Fcorr, 6), "\n\n"), sep="")

  } else {
    FtrainGeno = FtrainGenoOTH
    FtestGeno = FtestGenoOTH
    all_train_acc <- NULL
    all_val_acc <- NULL
    e_mat <- NULL
    beta <- NULL
    f_time <- Sys.time()
    cat(paste0("Start ", mm, " permutation. ", "\n"), sep="")
    for (ite in seq(1,per)){
      cat(paste0(" - ", mm, ": Permutation ", ite," ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"), sep="")
      set.seed(ite)
      train_samples <- sample(1:nrow(FtrainGeno), size = nrow(FtrainGeno)*0.8, replace = FALSE)
      train_geno <- FtrainGeno[train_samples,]
      val_geno <- FtrainGeno[-train_samples,]
      train_pheno <- FtrainPheno[train_samples]
      val_pheno <- FtrainPheno[-train_samples]
    
      fit = RRb_func(train_pheno, train_geno, val_geno)
      e_mat <- cbind(e_mat, fit$u)
      beta <- c(beta, fit$beta)
      train_e = rowMeans(e_mat)
      train_beta <- mean(beta)
    
      train_valid = data.matrix(train_geno) %*% as.matrix(train_e)
      train_pred <- train_valid + c(train_beta)
      train_acc <- cor(train_pred, train_pheno, use = "complete")
      val_valid = data.matrix(val_geno) %*% as.matrix(train_e)
      val_pred <- val_valid + c(train_beta)
      val_acc <- cor(val_pred ,val_pheno, use = "complete")
      all_train_acc <- c(all_train_acc, as.vector(train_acc))
      all_val_acc <- c(all_val_acc, as.vector(val_acc))
    }
    FtestGeno2 = (FtestGeno + 1)/3
    test_valid = data.matrix(FtestGeno2) %*% as.matrix(train_e)
    test_pred <- test_valid + c(train_beta)
    FinalCorr = cor(test_pred, FtestPheno)
    df1 <- data.frame(Phenotype = FtestPheno, GEBV = test_pred)
    df1 = data.frame(Sample=rownames(df1), Obs.Phenotype=df1[,1], Pred.Phenotype=df1[,2])
    df2 <- data.frame(Num=seq(1,length(names(train_e))), Marker=names(train_e), e=train_e, beta=train_beta)
    df3 <- data.frame(Repeat=seq(1,per), Train.Correlation=all_train_acc, Val.Correlation=all_val_acc)
    add_row1 = c("Train-Validation (Average)", paste(round(mean(all_train_acc), 6), " (",round(sd(all_train_acc),6),")", sep=""), paste(round(mean(all_val_acc), 6), " (",round(sd(all_val_acc),6),")", sep=""))
    add_row2 = c("Test", "-", paste(round(FinalCorr, 6), " (",round(sd(all_train_acc),6),")", sep=""))
    df3 = rbind(as.matrix(df3), add_row1)
    df3 = rbind(as.matrix(df3), add_row2)
    Total_list[[mm]]$Model1 = df1
    Total_list[[mm]]$Model2 = df2
    Total_list[[mm]]$Model3 = df3
    Total_list[[mm]]$TrainCorr = round(mean(all_train_acc),6)
    Total_list[[mm]]$ValCorr = round(mean(all_val_acc),6)
    Total_list[[mm]]$Bat = 1
    Total_list[[mm]]$Predict = test_pred
    Total_list[[mm]]$TestCorr = FinalCorr
    
    cat(paste0(" - ", mm, ": Final correlations for training and validation dataset are ", round(mean(all_train_acc),6), ", ", round(mean(all_val_acc),6),"\n"), sep="")
    cat(paste0(" - ", mm, ": Standard deviation for ", per, " replications is ", round(sd(all_train_acc),6), ", ", round(sd(all_val_acc),6), "\n"), sep="")
    cat(paste0(" - ", mm, ": Final correlation for test dataset is ", round(FinalCorr,6), "\n\n"), sep="")
  }
}


#(5) Cheking Best Model
cat(paste0("Final: Checking Best model.", "\n", sep=""))
TestCorr_vec = NULL
TrainCorr_vec = NULL
ValCorr_vec = NULL
MM_name = NULL
Batch_vec = NULL
for (mm in mms){
  Trainmm = Total_list[[mm]][["TrainCorr"]]
  Valnmm = Total_list[[mm]][["ValCorr"]]
  Testmm = Total_list[[mm]][["TestCorr"]]
  Batchmm = Total_list[[mm]][["Bat"]]
  TrainCorr_vec = c(TrainCorr_vec, Trainmm)
  ValCorr_vec = c(ValCorr_vec, Valnmm)
  TestCorr_vec = c(TestCorr_vec, Testmm)
  MM_name = c(MM_name, mm)
  Batch_vec = c(Batch_vec, Batchmm)
}
Max_idx = which(ValCorr_vec == max(ValCorr_vec))
Max_name = MM_name[Max_idx]
Max_TestCorr = TestCorr_vec[Max_idx]
Max_TrainCorr = TrainCorr_vec[Max_idx]
Max_ValCorr = ValCorr_vec[Max_idx]
cat(paste0(" - Best Model: ", Max_name, ", Train, validation and test correlation: ", round(Max_TrainCorr,6), ", ", round(Max_ValCorr,6), ", ", round(Max_TestCorr,6),"\n", sep=""))
for (idx in seq(1,length(MM_name))[-Max_idx]){
  cat(paste0(" - Other Models: ", MM_name[idx], ", Train, validation and test correlation: ", round(TrainCorr_vec[idx],6), ", ", round(ValCorr_vec[idx],6), ", ", round(TestCorr_vec[idx],6),"\n", sep=""))
}
Tdf = data.frame("Num"=seq(1:length(MM_name)), "Model"=MM_name, "Batch"=Batch_vec, "Train_Corr."=round(TrainCorr_vec,6), "Val_Corr."=round(ValCorr_vec, 6), "Test_Corr."=round(TestCorr_vec,6), check.names=FALSE)
Tdf = Tdf[order(Tdf[,5],decreasing = TRUE),]
Tdf[,1] = seq(1:length(MM_name))
SaveTdf = paste0("Final_modeling_results_for_all_models.txt")
write.table(Tdf, file = SaveTdf, row.names=F, quote=F, sep="\t")

cat(paste0(" - Saving all prediction models as h5 or rds.", "\n", sep=""))
cat(paste0(" - Saving all plot and correlation files.", "\n\n", sep=""))
for (mm in mms){
  if (mm == Max_name){
    Prefix = "Best_final_"
  } else {
    Prefix = "Final_"
  }
  if (mm == "CNN" | mm == "DNN" | mm == "RF"){
    if (mm == "CNN" | mm == "DNN"){
      MaxModel = Total_list[[mm]]$Model 
      MaxBatch = Total_list[[mm]]$Bat
      MaxTrainCor = Total_list[[mm]]$TrainCorr
      MaxTestCor = Total_list[[mm]]$TestCorr
      MaxValCor = Total_list[[mm]]$ValCorr
      SaveNameMax = paste0(Prefix, mm, "_prediction_max_model_", length(FinalMarker), "m_", round(Total_list[[mm]]$TestCorr * 100,0), "c.h5")
      MaxModel[[1]] %>% save_model_hdf5(SaveNameMax)
    }
    if (mm == "RF") {
      MaxModel = Total_list[[mm]]$Model
      MaxBatch = Total_list[[mm]]$Bat
      MaxTrainCor = Total_list[[mm]]$TrainCorr
      MaxTestCor = Total_list[[mm]]$TestCorr
      MaxValCor = Total_list[[mm]]$ValCorr
      SaveNameMax = paste0(Prefix, mm, "_prediction_max_model_", length(FinalMarker), "m_", round(Total_list[[mm]]$TestCorr * 100,0), "c.rds")
      saveRDS(MaxModel, SaveNameMax)
    }

    SaveGEBV = paste0(Prefix, mm, "_predPhenotype_list_", length(FinalMarker), "m_", round(Total_list[[mm]]$TestCorr * 100,0), "c.txt")
    SaveCor = paste0(Prefix, mm, "_permutation_results_", length(FinalMarker), "m_", round(Total_list[[mm]]$TestCorr * 100,0), "c.txt")
    df1 <- data.frame(Obs.Phenotype = FtestPheno, Pred.Phenotype = Total_list[[mm]]$Predict)
    df1 = data.frame(Sample=rownames(df1), Obs.Phenotype=df1[,1], Pred.Phenotype=df1[,2])
    write.table(df1, file = SaveGEBV, row.names=F, quote=F, sep="\t")
  
    df2 <- data.frame(Repeat=seq(1,length(Total_list[[mm]]$TrainCorrlist)), Train.Correlation=Total_list[[mm]]$TrainCorrlist, Val.Correlation=Total_list[[mm]]$ValCorrlist)
    add_row1 = c("Validation (Max)", paste(round(MaxTrainCor, 6), " (",round(Total_list[[mm]]$MedCorrMADT,6),")", sep=""),  paste(round(MaxValCor, 6), " (",round(Total_list[[mm]]$MedCorrMADV,6),")", sep=""))
    add_row2 = c("Test (of Max)", "-", round(MaxTestCor, 6))
    df2 = rbind(df2, add_row1, add_row2)
    write.table(df2, file = SaveCor, row.names=F, quote=F, sep="\t")
    
    my.formula <- y~x
    g0 <- ggplot(df1, aes(x=Obs.Phenotype, y = Pred.Phenotype)) + stat_smooth(method = "lm", se = FALSE, color = "black", formula = my.formula) + ylab("Pred. phenotype") + xlab("Obs. phenotype") +
      stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), size=6, parse = TRUE) + geom_point() +
      theme_bw() +
      theme(axis.text = element_text(size = 18),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=22))
    ggsave(filename = paste0(Prefix, mm, "_predPhenotype_plot.png"), plot = g0, dpi = 600)
    
  } else {
    SavePar = paste0(Prefix, mm, "_prediction_model_", length(FinalMarker), "m_", round(FinalCorr * 100,0), "c.txt")
    SaveCor = paste0(Prefix, mm, "_permutation_results_", length(FinalMarker), "m_", round(FinalCorr * 100,0), "c.txt")
    SaveGEBV = paste0(Prefix, mm, "_predPhenotype_list_", length(FinalMarker), "m_", round(FinalCorr * 100,0), "c.txt")
    write.table(Total_list[[mm]]$Model1, file = SaveGEBV, row.names=F, quote=F, sep="\t")
    write.table(Total_list[[mm]]$Model2, file = SavePar, row.names=F, quote=F, sep="\t")
    write.table(Total_list[[mm]]$Model3, file = SaveCor, row.names=F, quote=F, sep="\t")

    my.formula <- y~x
    g0 <- ggplot(Total_list[[mm]]$Model1, aes(x=Obs.Phenotype, y = Pred.Phenotype)) + stat_smooth(method = "lm", se = FALSE, color = "black", formula = my.formula) + ylab("Pred. phenotype") + xlab("Obs. phenotype") +
      stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), size=6, parse = TRUE) + geom_point() +
      theme_bw() +
      theme(axis.text = element_text(size = 18),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=22))
    ggsave(filename = paste0(Prefix, mm, "_predPhenotype_plot.png"), plot = g0, dpi = 600)
  }
}


#(6) Final Marker distribution plot and marker information
if (!is.null(args$gwas)){
  snp_cn=1
  chr_cn = 2
  pos_cn = 3
  pv_cn = 4
  gwas_results <- fread(gwasfile, header=T, check.names=FALSE, data.table=FALSE)
  sel_ix <- match(FinalMarker, gwas_results[,snp_cn])
  df2 <- data.frame(Marker = gwas_results[,snp_cn][sel_ix], Chromosome = gwas_results[sel_ix, chr_cn], Position = gwas_results[sel_ix, pos_cn], "-log(P-value)" = log10(gwas_results[sel_ix, pv_cn]), check.names=F)
  df3 <- as.data.frame(table(factor(df2$Chromosome, levels = c(1:max(gwas_results[,chr_cn])))))
  colnames(df3) <- c("Chromosome", "Count")
  g2 <- ggplot(df3, aes(Chromosome, Count)) +
    geom_bar(stat = "identity", width = .6, fill = "tomato3") +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), plot.background = element_blank()) +
    geom_text(aes(label = Count), vjust = -0.3, size = 5) + theme_classic(base_size=22)  ## origianl size 3.5
  g2
  marker_select_save = paste0("Best_final_", Max_name, "_selected_markers_Chromosomal_distribution.png", sep="")
  ggsave(filename = marker_select_save, plot = g2, dpi = 600)
  write.table(df2, file = paste0("Best_final_", Max_name, "_selected_markers_Info.txt"), row.names=F, quote=F, sep="\t")
}

if (!is.null(args$info)){  
  infotable = fread(infofile, header=T, check.names=FALSE, data.table=FALSE)
  efftable = fread("..\\GWAS_marker_effect_ordered.csv", header=T, check.names=FALSE, data.table=FALSE)
  df2 <- cbind(FinalMarker, infotable[match(FinalMarker, infotable[,1]),c(2,3)], efftable[match(FinalMarker, efftable[,1]),2])
  colnames(df2) <- c("Marker", "Chromosome", "Position", "Effect")
  df3 <- data.frame(table(factor(df2$Chromosome, levels = c(1:max(df2$Chromosome)))))
  colnames(df3) <- c("Chromosome", "Count")
  g2 <- ggplot(df3, aes(Chromosome, Count)) +
    geom_bar(stat = "identity", width = .6, fill = "tomato3") +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), plot.background = element_blank()) +
    geom_text(aes(label = Count), vjust = -0.3, size = 5) + theme_classic(base_size=22)  ## origianl size 3.5
  marker_select_save =  paste0("Best_final_",  Max_name, "_selected_markers_Chromosomal_distribution.png", sep="")
  ggsave(filename = marker_select_save, plot = g2, dpi = 600)
  write.table(df2, file = paste0("Best_final_", Max_name, "_selected_markers_Info.txt"), row.names=F, quote=F, sep="\t")
}

if (!is.null(args$gwas) & !is.null(args$info)){
  efftable = fread("..\\GWAS_marker_effect_ordered.csv", header=T, check.names=FALSE, data.table=FALSE)
  df2 <- data.frame(Num=seq(1:length(FinalMarker)), Marker=selected_markers, Effect=efftable[match(selected_markers, efftable[,1]),2])
  write.table(df2, file = paste0("Best_final_", Max_name, "_selected_markers_Info.txt"), row.names=F, quote=F, sep="\t")
}

cat(paste0("END: Modeling and prediction.","\n"), sep="")
cat("\n\n")

#################################################################
