RRb_func = function (train_pheno, train_geno, test_geno) {
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_BLUP <- mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
  train_e = as.matrix(train_BLUP$u)
  train_valid = train_geno %*% train_e
  train_pred = train_valid + c(train_BLUP$beta)
  val_valid = test_geno %*% train_e
  val_pred <- val_valid + c(train_BLUP$beta)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=train_BLUP, "u"=train_BLUP$u, "beta"=train_BLUP$beta)
  return(return_value)
}
