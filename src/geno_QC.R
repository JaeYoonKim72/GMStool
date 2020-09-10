geno_QC <- function(genotype, phenotype, maf_cutoff = 0.01, miss_cutoff = 0.01){
  sample_miss <- apply(genotype, 1, function(x) length(which(is.na(x))))/ncol(genotype)
  snp_miss <- apply(genotype, 2, function(x) length(which(is.na(x))))/nrow(genotype)
  n0 <- apply(genotype==0,1,sum,na.rm=T)
  n1 <- apply(genotype==1,1,sum,na.rm=T)
  n2 <- apply(genotype==2,1,sum,na.rm=T)
  n <- n0 + n1 + n2
  p <- ((2*n0)+n1)/(2*n)
  q <- 1 - p
  maf <- pmin(p, q)
  maf_ix <- which(maf < maf_cutoff)
  miss_sample_ix <- which(sample_miss >= miss_cutoff)
  snp_rm_ix <- union(maf_ix, which(snp_miss >= miss_cutoff))
  if (length(miss_sample_ix) > 0){
    genotype <- genotype[-miss_sample_ix,]
    phenotype <- phenotype[-miss_sample_ix]
    #pheno1 <- pheno1[-miss_sample_ix,]
    cat(paste0("- ", length(miss_sample_ix), " samples were filtered", '\n'), sep="")
  } else {
    cat(paste0("- ", "No samples were filtered.", "\n"), sep="")
  }
  if (length(snp_rm_ix) > 0){
    genotype <- genotype[,-snp_rm_ix]
    ix <- ix[which(ix %in% colnames(genotype))]
    cat(paste0("- ", length(snp_rm_ix), " SNPs were filtered", '\n\n'), sep="")
  } else{
    cat(paste0("- No SNPs were filtered.", "\n\n"), sep="")
  }
  return(list(genotype = genotype, phenotype = phenotype))
}
