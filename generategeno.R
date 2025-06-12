
# Function for testing individual PGS variability with different sample sizes, heritability, # of SNPs

PGS.uncertainty <- function(N, M, h2, S = -1){
  # N: # of individuals
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var.neg = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas.neg = rnorm(M, mean = 0, sd = sqrt(var.neg)) # draw betas with variance var.neg
  
  beta.neg.est = betas.neg + 1/N * rnorm(M, 0, 1) # add noise proportional to sample size
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  
  pgs = X %*% beta.neg.est # compute PGS accounting for negative selection
  
  ## WRONG
  pgs.var = 1/N * rowSums(X^2) # variance in individual PGS estimates due to noise in the effect size estimates

  
  return(list(pgs = pgs, variance = pgs.var))
}



PGS.sample.uncertainty <- function(N, M, h2, S = -1, n.samples = 1000){
  # N: # of individuals
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  # n.samples: # of samples where PGS and its variance are computed from
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var.neg = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas.neg = rnorm(M, mean = 0, sd = sqrt(var.neg)) # draw true betas with variance var.neg
  
  beta.samples = matrix(nrow = n.samples, ncol = M) # collect n.samples of beta estimates
  
  for (i in 1:n.samples){
    beta.samples[i,] = betas.neg + 1/N * rnorm(M, 0, 1) # add noise proportional to sample size
  }
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  beta.mean = colMeans(beta.samples) # mean of the sampled beta estimates, acts as the point estimate
  pgs.samples = matrix(nrow = n.samples, ncol = N) # each column is a sample of PGSs for an individual
  
  pgs.samples = beta.samples %*% t(X)

  
  pgs = colMeans(pgs.samples) # point estimates for PGSs
  pgs.var = apply(pgs.samples, 2, var)
  
  return(list(pgs = pgs, variance = pgs.var))
}

