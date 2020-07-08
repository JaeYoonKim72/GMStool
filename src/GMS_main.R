GMS_main <- function(ini_snps_bk, init_selsnp, j, mm, cv_samples, geno2, phenotype1, preset_fname, ix){
  library(caret)
  library(rrBLUP)
  sink(file = paste0("CV-", j, "_", mm, ".log"))
  ini_snps = ini_snps_bk
  o_delta <- 1
  sel_snps <- init_selsnp
  cat(paste0("CV Number ", j, ": Start","\n"), sep="")
  f_time <- Sys.time()
  #(8-1)Cheking prediction (correlation) rate for all markers

  cat(paste0("CV Number ", j, ": Check prediction (correlation) rate for all markers using ", mm, " model.","\n"), sep="")
  train_samples <- which(cv_samples != j)
  train_geno <- geno2[train_samples,]
  test_geno <- geno2[-train_samples,]
  train_pheno <- phenotype1[train_samples]
  test_pheno <- phenotype1[-train_samples]

  nTrainSam = length(train_samples)
  nTestSam = length(test_pheno)
  nGeno = dim(train_geno)[2]
  nSamGen = c(nTrainSam, nTestSam, nGeno)
  #nSamGenSum[[paste0("CV-", j, "_", mm)]] = nSamGen


  if (allm == TRUE){
    if (mm == "RRblup"){
      fit = RRb_func(train_pheno, train_geno, test_geno)
      pred_pheno = fit$predicted
      fit_model = fit$model
    }
    if (mm == "RF"){
      fit = RF_func(train_pheno, train_geno, test_geno)
      pred_pheno = fit$predicted
      fit_model = fit$model
    }
    train_acc <- cor(pred_pheno, test_pheno, use = "complete")
    #all_train_acc <- c(all_train_acc, as.vector(train_acc))
  } else {
    train_acc <- "-"
    #all_train_acc <-c(all_train_acc, as.vector(train_acc))
  }
  cat(paste0("CV Number ", j, ": Prediction (correlation) rate for all markers using ", mm, " is ", train_acc, "\n"), sep="")
  cat(paste0("CV Number ", j, ": Prediction running time for all markers is ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"), sep="")
  #inisum[[paste0("CV-", j, "_", mm)]] <- format(difftime(Sys.time(), f_time), usetz = TRUE)
  inisum1 <- format(difftime(Sys.time(), f_time), usetz = TRUE)

  #(8-2)Selecting initial markers
  ix1 <- ix
  if(length(preset_fname) != 0){
    preset_idx = match(preset_fname, ix1)
    ix1 <- ix1[-preset_idx]
    pp_ix = preset_fname
    if (ini_snps > 0) {
      if (mm == "RF"){
        init_set <- ix1[1:ini_snps]
        train_ix <- c(pp_ix, init_set)
      } else {
        init_set_t <- ix1[1:ini_snps]
        init_set <- c(pp_ix, init_set_t)
        zv <- apply(train_geno[,init_set], 2, function(x) length(unique(x)) == 1)
        if (sum(cv) == 1){
          highCorr = NULL
        } else {
          cor_mat <- cor(train_geno[,init_set][, !zv])
          highCorr <- findCorrelation(cor_mat, 0.9)
        }
        if (length(highCorr) == 0){
          train_ix <- init_set
        } else {
          init_set_tx = match(init_set_t, init_set)
          init_set_tx_int = intersect(init_set_tx, highCorr)
          if (length(init_set_tx_int) > 0 ){
            train_ix <- init_set[-init_set_tx_int]
          } else {
            train_ix <- init_set
          }
        }
      }
    } else {
      train_ix <- pp_ix
    }
    rm(pp_ix)
  } else {
    if (mm == "RF"){
      init_set <- ix1[1:ini_snps]
      train_ix = init_set
    } else {
      init_set <- ix1[1:ini_snps]
      zv <- apply(train_geno[,init_set], 2, function(x) length(unique(x)) == 1)
      if (sum(cv) == 1){
        highCorr = NULL
      } else {
        cor_mat <- cor(train_geno[,init_set][, !zv])
        highCorr <- findCorrelation(cor_mat, 0.9)
      }
      if (length(highCorr) == 0){
        train_ix <- ix1[1:ini_snps]
      } else {
        train_ix <- init_set[-highCorr]
      }
    }
  }

  r=0
  except_ix <- NULL
  except_tix <- NULL
  except_ixc <- 0
  add_ix <- train_ix
  result_acc <- 0
  ms_out <- list()
  t_add <- train_ix
  train_itr_geno <- train_geno[,train_ix]
  test_itr_geno <- test_geno[,train_ix]
  invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))

  #(8-3)Cheking prediction (correlation) rate for initial markers
  if (mm == "RRblup"){
    fit_itr = RRb_func(train_pheno, train_itr_geno, test_itr_geno)
    pred_itr_pheno = fit_itr$predicted
    fit_itr_model = fit_itr$model
  }
  if (mm == "RF"){
    fit_itr = RF_func(train_pheno, train_itr_geno, test_itr_geno)
    pred_itr_pheno = fit_itr$predicted
    fit_itr_model = fit_itr$model
  }

  train_itr_acc <- cor(pred_itr_pheno, test_pheno, use = "complete")

  #(8-4)Checking initial result
  ms_out[[1]] <- c(train_itr_acc, train_ix)
  k=2
  cat(paste0("CV Number ", j, ": Initial marker set's prediction (correlation) rate is ", round(train_itr_acc,6), "\n"), sep="")
  result_acc <- c(result_acc, train_itr_acc)

  #(8-5)Starting marker selection
  i=1
  max_noid <- 0
  max_noid2 <- 0
  s_time <- Sys.time()
  cat(paste0("CV Number ", j, ": Start marker selection at ", format(s_time, usetz = T), "\n"), sep="")


  while (train_itr_acc < acc1){
    if (difftime(Sys.time(), s_time, units = "hours") > time_cutoff){
      cat("STOP: Computational time is over, so stop this CV", '\n\n')
      break
    }
    if (max_noid > 0 & max_noid %% 5000 == 0){
      cat("STOP: Low delta counts reached 5000 counts, so stop this CV", '\n\n')
      break
    }

    #(8-6) Selecting SNPs
    if (mm == "RF"){
      add_snps = ix1[(ini_snps + (sel_snps*(i-1)+1)):(ini_snps + (sel_snps*(i)))]
      new_train_ix <- c(train_ix, add_snps)
    } else {
      add_snps_t = ix1[(ini_snps + (sel_snps*(i-1) +1)):(ini_snps + (sel_snps*(i)))]
      init_set <- c(train_ix, add_snps_t)
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

    #(8-7) Cheking prediction (correlation) rate for selected markers
    train_itr_geno <- train_geno[,new_train_ix]
    test_itr_geno <- test_geno[,new_train_ix]
    if (mm == "RRblup"){
      fit_itr = RRb_func(train_pheno, train_itr_geno, test_itr_geno)
      pred_itr_pheno = fit_itr$predicted
      fit_itr_model = fit_itr$model
    }
    if (mm == "RF"){
      fit_itr = RF_func(train_pheno, train_itr_geno, test_itr_geno)
      pred_itr_pheno = fit_itr$predicted
      fit_itr_model = fit_itr$model
    }

    acc_test <- cor(pred_itr_pheno, test_pheno, use = "complete")

    if (acc_test - train_itr_acc <= delta){
      except_ix <- c(except_ix, add_snps)
      max_noid <- max_noid + 1
      max_noid2 <- max_noid2 + 1
    } else {
      ms_out[[k]] <- c(acc_test, add_snps)
      train_ix <- c(train_ix, add_snps)
      t_add <- c(t_add, add_snps)
      train_itr_acc <- acc_test
      add_ix <- c(add_ix, add_snps)
      result_acc <- c(result_acc, train_itr_acc)
      o_delta <- result_acc[length(result_acc)] - result_acc[length(result_acc)-1]
      cat(paste0("Add: ", round(train_itr_acc,6)," (delta: ", round(o_delta, 6), "), time: ", format(difftime(Sys.time(), s_time), usetz = TRUE), "\n", sep=""))
      max_noid <- 0
      max_noid2 <- 0
      r=0
      sel_snps <- init_selsnp
    }

    #(8-8) Setting nex iteration
    if (mm == "RF"){
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
            set.seed(10)
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
            set.seed(10)
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

    if (mm == "RRblup" ) {
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


  #(8-9) Saving and iteration result
  #selected_train_acc <- c(selected_train_acc, train_itr_acc)
  invisible(gc(verbose=FALSE, reset=TRUE, full=TRUE))
  cat(paste0("CV Number ", j, ": # of selected markers is ", length(train_ix), ", and correlation rate is ", train_itr_acc, "\n"), sep="")
  cat(paste0("CV Number ", j, ": Running time for marker selection is ", format(difftime(Sys.time(), f_time), usetz = TRUE), "\n"),sep="")
  cat(paste0("CV Number ", j, ": End ", mm, " model."  ,"\n\n\n"), sep="")
  #tsum[[paste0("CV-", j, "_", mm)]] = format(difftime(Sys.time(), f_time), usetz = TRUE)
  tsum1 = format(difftime(Sys.time(), f_time), usetz = TRUE)

  CV_name = paste0("CV-", j, "_", mm)
  #CV_results[[CV_name]] = train_ix

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
  sink(file = NULL)
  return(list(all_train_acc = train_acc,
              selected_train_acc = train_itr_acc,
              CV_results = train_ix,
              nSamGenSum = nSamGen,
              tsum = tsum1,
              inisum = inisum1))
}
