p_impute <- function(genotype) {

  #calculate probabilities of the three genotypes for each column (SNP)
  n <- nrow(genotype)
  p_0 <- apply(genotype, 2, function(x){sum(x == 0, na.rm = T)/(n - sum(is.na(x)))})
  p_1 <- apply(genotype, 2, function(x){sum(x == 1, na.rm = T)/(n - sum(is.na(x)))})
  p_2 <- apply(genotype, 2, function(x){sum(x == 2, na.rm = T)/(n - sum(is.na(x)))})
  p <- data.frame(p_0, p_1, p_2)
  
  #make sure the probalilities add up to one so that it won't be a problem for the sampling function
  p <- t(apply(p, 1, function(x){x/sum(x)}))
  
  #make a table for indices of missing genotypes
  NA_indices <- which(is.na(genotype), arr.ind = T)
  
  #replace missing genotypes by sampling from (0,1,2) based on probabilities given in table p
  for (i in 1:nrow(NA_indices)) {
    x <- NA_indices[i, ]
    genotype[x[1], x[2]] <- sample(c(0:2), 1, replace = T, prob = p[x[2], ])
  }
  
  genotype
}