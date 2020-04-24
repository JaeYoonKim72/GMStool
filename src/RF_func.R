RF_func = function (train_pheno, train_geno, test_geno){
  set.seed(10)
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=dim(train_geno)[1], ntrees=100, metric="RMSE")
  pred_pheno <- predict(train_RF, newdata = test_geno)
  return_value = list("predicted"=pred_pheno, "model"=train_RF)
  return(return_value)
}
