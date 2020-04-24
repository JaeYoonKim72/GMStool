RRb_func = function (train_pheno, train_geno, test_geno) {
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_BLUP <- mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
  train_e = as.matrix(train_BLUP$u)
  pheno_valid = test_geno %*% train_e
  pred_pheno <- pheno_valid + c(train_BLUP$beta)
  return_value = list("predicted"=pred_pheno, "model"=train_BLUP, "u"=train_BLUP$u, "beta"=train_BLUP$beta)
  return(return_value)
}
