RF_func = function (train_pheno, train_geno, test_geno, check){
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  set.seed(200)
  if (check == "FM"){
    train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]/3), ntree=1000, metric="RMSE", maxnodes=NULL)
  } else {
    train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=round(dim(train_geno)[1]), ntree=100, metric="RMSE", maxnodes=NULL)
  }
  train_pred <- predict(train_RF, newdata = train_geno)
  val_pred <- predict(train_RF, newdata = test_geno)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=train_RF)
  return(return_value)
}
