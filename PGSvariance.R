
# Function for testing individual PGS variability with different sample sizes, heritability, # of SNPs

PGS.uncertainty <- function(N, M, h2, S = -1){
  # N: # of individuals in GWAS used to estimate effect sizes
  # n: # of individuals whose PGS is to be computed
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  
  ps = runif(M, 0.01, 0.99) #  allele 1 frequencies for M SNPs
  var = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas = rnorm(M, mean = 0, sd = sqrt(var)) # draw betas with variance var.neg
  
  sigma.est = 1/N *(1 - var) # effects of SNPs can be estimated with variance proportional to their allele frequency (rare alleles have larger effect sizes, and larger effect sizes can be estimated with smaller uncertainty)
  
  beta.est = rnorm(M, betas, sqrt(sigma.est)) # add noise proportional to sample size, formula from paper 3
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes for n individuals assuming allele frequencies follow HWE
  idx.filter = apply(X, 2, sd) > 0 # filter out SNPs with no variation
  X = X[, idx.filter]
  X = scale(X, center = TRUE, scale = TRUE)
  beta.est = beta.est[idx.filter]
  sigma.est = sigma.est[idx.filter]
  
  pgs = X %*% beta.est # compute PGS accounting for negative selection
  #return(dim(sigma.est))
  pgs.var = X^2 %*% sigma.est # variance in individual PGS estimates due to noise in the effect size estimates
  
  return(list(pgs = pgs, variance = pgs.var))
}

PGS.uncertainty.liability <- function(N, M, h2, S = -1, K){
  # N: # of individuals
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  # K: disease prevalence
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var.neg = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas.neg = rnorm(M, mean = 0, sd = sqrt(var.neg)) # draw betas with variance var.neg
  
  t  <- qnorm(1-K) # threshold liability for disease based on prevalence K. People with liability > t have the disease
  z = dnorm(t) # density of standard normal distribution at t
  betas.scaled = betas.neg*K*(1-K)/z # transform betas to liability scale
  sigma.est = 1/N *(1 - var.neg)
  beta.neg.est = rnorm(M, betas.scaled, sqrt(sigma.est)) # add noise proportional to sample size, formula from paper 3
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  idx.filter = apply(X, 2, sd) > 0
  X = X[, idx.filter]
  X = scale(X, center = TRUE, scale = TRUE)
  beta.neg.est = beta.neg.est[idx.filter]
  
  pgs = X %*% beta.neg.est # compute PGS accounting for negative selection
  
  #sigma.est = 1/N *(1 - h2/M) 
  pgs.var = sigma.est * rowSums(X^2) # variance in individual PGS estimates due to noise in the effect size estimates
  
  return(list(pgs = pgs, variance = pgs.var))
}

PGS.sample.uncertainty.liability <- function(N, M, h2, S = -1, K, n.samples = 1000){
  # N: # of individuals
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  # K: disease prevalence
  # n.samples: # of samples where PGS and its variance are computed from
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas = rnorm(M, mean = 0, sd = sqrt(var.neg)) # draw true betas with variance var.neg
  
  t  <- qnorm(1-K) # threshold liability for disease based on prevalence K. People with liability > t have the disease
  z = dnorm(t) # density of standard normal distribution at t
  betas.scaled = betas*K*(1-K)/z # transform betas to liability scale
  beta.samples = matrix(nrow = n.samples, ncol = M) # collect n.samples of beta estimates
  
  for (i in 1:n.samples){
    beta.samples[i,] = rnorm(M, betas.scaled, sqrt(1/N *(1 - var))) # add noise proportional to sample size, formula from paper 3
  }
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  
  idx.filter = apply(X, 2, sd) > 0 # filter out variants with no variation across individuals
  X.filtered = X[, idx.filter]
  beta.samples = beta.samples[, idx.filter] 
  X.scaled = scale(X.filtered) # standardize genotypes to have mean 0 and variance 1
  beta.mean = colMeans(beta.samples) # mean of the sampled beta estimates, acts as the point estimate
  pgs.samples = matrix(nrow = n.samples, ncol = N) # each column is a sample of PGSs for an individual
  
  pgs.samples = beta.samples %*% t(X.scaled)
  
  
  pgs = colMeans(pgs.samples) # point estimates for PGSs
  pgs.var = apply(pgs.samples, 2, var) # variances of individual PGS estimates
  
  return(list(pgs.samples = pgs.samples, pgs = pgs, variance = pgs.var))
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
    beta.samples[i,] = rnorm(M, betas.neg, sqrt(1/N *(1 - h2/M))) # add noise proportional to sample size, formula from paper 3
  }
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  
  idx.filter = apply(X, 2, sd) > 0 # filter out variants with no variation across individuals
  X.filtered = X[, idx.filter]
  beta.samples = beta.samples[, idx.filter] 
  X.scaled = scale(X.filtered) # standardize genotypes to have mean 0 and variance 1
  beta.mean = colMeans(beta.samples) # mean of the sampled beta estimates, acts as the point estimate
  pgs.samples = matrix(nrow = n.samples, ncol = N) # each column is a sample of PGSs for an individual
  
  pgs.samples = beta.samples %*% t(X.scaled)

  
  pgs = colMeans(pgs.samples) # point estimates for PGSs
  pgs.var = apply(pgs.samples, 2, var) # variances of individual PGS estimates
  
  return(list(pgs = pgs, variance = pgs.var))
}

PGS.sample.rankings <- function(N, M, n, h2, n.samples, S = -1, t = 0.9){
  # N: # of individuals in GWAS used to estimate effect sizes
  # M: # of SNPs
  # n: # of individuals whose PGS is to be computed
  # h2: heritability
  # n.samples: # of PGS estimates to be sampled
  # S: relative contributions of common vs rare variants
  # q: PGS quantile to be considered
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas = rnorm(M, mean = 0, sd = sqrt(var)) # draw betas with variance var.neg
  
  X <- matrix(rbinom(n * M, size = 2, prob = rep(ps, each = n)), nrow = n, ncol = M) # generate genotypes for n individuals assuming allele frequencies follow HWE (n x M)
  
  # filter out SNPs with no variation
  idx.filter = apply(X, 2, sd) > 0
  X = X[, idx.filter]
  X = scale(X, center = TRUE, scale = TRUE)
  betas = betas[idx.filter] 
  
  M.true = length(betas)
  
  sigma.est = 1/N *(1 - var) # effects of SNPs can be estimated with variance proportional to their allele frequency (rare alleles have larger effect sizes, and larger effect sizes can be estimated with smaller uncertainty)
  sigma.est = sigma.est[idx.filter]
  
  # sample effect estimates
  beta.samples = matrix(nrow = n.samples, ncol = M.true) # collect n.samples of beta estimates (n.samples x M (or less if filtered SNPs))
  for (i in 1:n.samples){ # sample beta estimates and collect into beta.samples matrix
    beta.samples[i,] = rnorm(M.true, betas, sqrt(sigma.est)) # add noise proportional to sample size, formula from paper 3
  }
  
  # compute PGS estimates
  pgs.est = X %*% t(beta.samples) # PGS samples (n x n.samples)
  pgs.mean = colMeans(pgs.est) # PGS point estimates for each individual
  
  # compute rankings of PGS estimates
  rankings = matrix(nrow = n, ncol = n.samples) # a (n x n.samples) matrix of rankings
  
  for (i in 1:n.samples) {
    ind = 1:n # indices of PGS point estimates
    rankings[, i] = as.integer(rank(pgs.est[, i], ties.method = "first")) # order indices from smallest PGS to largest
      #order(pgs.est[,i], decreasing = FALSE) 
    
  }
  t.ind = ceiling(t*n) # cutoff index
  rank.prop = apply(rankings, 1, function(x){mean(x > t.ind)x
    #return(mean(above))
  })
  
  return(list(
    geno = X, # a (n x M) matrix of genotypes of n individuals on M SNPs
    beta = betas, # a vector of the true effects of M SNPs
    p = ps, # a vector of allele frequencies at M SNPs
    pgs = pgs.mean, # a vector of PGS point estimates for n individuals
    pgs.est = pgs.est, # a (n x n.samples) matrix of PGS samples
    pgs.rankings = rankings, # a (n x n.samples) matrix of rankings of n.samples PGS estimates for n individuals
    prop.over.t = rank.prop # a vector of proportions of PGS values over t:th quantile for each individual
    )) 
}
  
ordered.estimates <- function(pgs.est, pgs.var, t, rho){ # uses confidence intervals
  # pgs.est: list of PGS point estimates for each individual
  # pgs.var: list of variance of the PGS point estimate for each individual
  # t: above PGS quantile to be considered
  # rho: level of confidence intervals
  
  pgs.est = as.numeric(pgs.est)
  pgs.var = as.numeric(pgs.var)
  N = length(pgs.est)
  ind = 1:N # indices of PGS point estimates
  ind = order(pgs.est, decreasing = FALSE) # order indices from smallest PGS to largest
  q.ind = ceiling(t * N) # cutoff index
  
  above.ind = ind[(q.ind + 1):N] # top (1-t)*100 % individuals
  below.ind = ind[1:q.ind] # individuals below q.ind
  
  above.pgs = pgs.est[above.ind] # highest PGS estimates
  above.pgs.var = pgs.var[above.ind] # variances of highest PGS estimates
  below.pgs = pgs.est[below.ind] # lowest t*100 % PGS estimates
  below.pgs.var = pgs.var[below.ind] # variances of lowest PGS estimates
  
  # confidence intervals:
  z = qnorm((1 + rho)/2) # z = 1.96 for 95 % CI
  above.pgs.lower = above.pgs - z*sqrt(above.pgs.var) # lower end
  above.pgs.upper = above.pgs + z*sqrt(above.pgs.var) # upper end
  
  below.pgs.lower = below.pgs - z*sqrt(below.pgs.var) # lower end
  below.pgs.upper = below.pgs + z*sqrt(below.pgs.var) # upper end
  
  cutoff.val = quantile(pgs.est, probs = t) # value at t:th quantile
  
  # proportion of individuals having CIs entirely above the t-quantile cutoff
  ci.above.t = above.pgs.lower > cutoff.val # number of lower CIs above t:th quantile, returns a logical vector
  proportion.certain.above = mean(ci.above.t) # proportion of pgs CIs entirely above t among those pgs estimates above t
  if (is.null(proportion.certain.above)) {
    proportion.certain.above = 0.0
  }
  # proportion of individuals having CIs entirely below the t-quantile cutoff
  ci.below.t = below.pgs.upper < cutoff.val # number of upper CIs below t:th quantile, returns a logical vector
  proportion.certain.below = mean(ci.below.t) # proportion of pgs CIs entirely below t among those pgs estimates below t
  if (is.null(proportion.certain.below)) {
    proportion.certain.below = 0.0
  }
 
  return(list(
    above.pgs = above.pgs, # above t % PGS estimates ...
    above.pgs.var = above.pgs.var, # ... and their corresponding variances
    above.pgs.lower = above.pgs.lower, # lower end of rho level CIs
    above.pgs.upper = above.pgs.upper, # upper end of rho level CIs
    cutoff.val = cutoff.val, # PGS value at t:th quantile
    prop.certain.above = proportion.certain.above, # proportion of PGSs certainly above t ...
    certain.above.indices = above.ind[ci.above.t], # ... and their indices
    below.pgs = below.pgs, # below t % PGS estimates ...
    below.pgs.var = below.pgs.var, # ... and their corresponding variances
    below.pgs.lower = below.pgs.lower, # lower end of rho level CIs
    below.pgs.upper = below.pgs.upper, # upper end of rho level CIs
    prop.certain.below = proportion.certain.below, # proportion of PGSs certainly below t ...
    certain.below.indices = below.ind[ci.below.t] # ... and their indices
  
  ))
  
}

compute.rank.correlations <- function(pgs.ranks, n.pairs = 1000, seed = 42) {
  set.seed(seed)
  #pgs.ranks = t(pgs.ranks) # n.sample rankings for the n indviduals, t(n x n.samples) matrix = (n.samples x n) matrix
  n.samples = ncol(pgs.ranks) # of samples (ranks) for each individual
  
  all.pairs = combn(n.samples, 2) # all unique combinations of indices of samples (2 x jotain) matrix
  nof.pairs = ncol(all.pairs) # of unique pairs
  sample.ind = sample(nof.pairs, min(n.pairs, length(all.pairs))) # samples the n.pairs pairs from the possible pairs nof.pairs
  
  sampled.pairs = as.matrix(all.pairs[,sample.ind, drop = FALSE]) # n.pairs based on sample.ind (2 x n.pairs) matrix
  
  corrs = apply(sampled.pairs, 2, function(pair){cor(pgs.ranks[, pair[1]], pgs.ranks[, pair[2]], method = "spearman")}) # compute Spearman correlations for each pairs
  return(corrs)
}

generate.geno <- function(M, n, h2, S = -1) {
  # M: # of SNPs
  # n: # of individuals whose PGS is to be computed
  # h2: heritability
  # S: relative contributions of common vs rare variants
  
ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
var = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
betas = rnorm(M, mean = 0, sd = sqrt(var)) # draw betas with variance var

X <- matrix(rbinom(n * M, size = 2, prob = rep(ps, each = n)), nrow = n, ncol = M) # generate genotypes for n individuals assuming allele frequencies follow HWE (n x M)

# filter out SNPs with no variation
idx.filter = apply(X, 2, sd) > 0
X = X[, idx.filter]
X = scale(X, center = TRUE, scale = TRUE)
betas = betas[idx.filter] 
var = var[idx.filter]

return(list(
  X = X, # a (n x M) matrix of genotypes
  var = var, # a M-vector of variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  beta = betas # a M-vector of true effect sizes at M SNPs
))
}

PGS.sample.rankings.2 <- function(N, n.samples, X, var, beta, t = 0.9){
  # N: # of individuals in GWAS used to estimate effect sizes
  # n.samples: # of PGS estimates to be sampled
  # X: # a (n x M) matrix of genotypes
  # var: a M-vector of variance proportional to negative selection coefficient for each M SNPs
  # beta: a M-vector of true effect sizes at M SNPs
  # t: PGS quantile to be considered
  
  M = ncol(X)
  n = nrow(X)
  sigma.est = 1/N *(1 - var) # effects of SNPs can be estimated with variance proportional to their allele frequency (rare alleles have larger effect sizes, and larger effect sizes can be estimated with smaller uncertainty)
  
  # sample effect estimates
  beta.samples = matrix(nrow = n.samples, ncol = M) # collect n.samples of beta estimates (n.samples x M (or less if filtered SNPs))
  for (i in 1:n.samples){ # sample beta estimates and collect into beta.samples matrix
    beta.samples[i,] = rnorm(M, beta, sqrt(sigma.est)) # add noise proportional to sample size, formula from paper 3
  }
  
  # compute PGS estimates
  pgs.est = X %*% t(beta.samples) # PGS samples (n x n.samples)
  pgs.mean = colMeans(pgs.est) # PGS point estimates for each individual
  
  # compute rankings of PGS estimates
  rankings = matrix(nrow = n, ncol = n.samples) # a (n x n.samples) matrix of rankings
  
  for (i in 1:n.samples) {
    ind = 1:n # indices of PGS point estimates
    rankings[, i] = as.integer(rank(pgs.est[, i], ties.method = "first")) # order indices from smallest PGS to largest
    #order(pgs.est[,i], decreasing = FALSE) 
    
  }
  t.ind = ceiling(t*n) # cutoff index
  rank.prop = apply(rankings, 1, function(x){mean(x > t.ind)
    #return(mean(above))
  })
  
  return(list(
    p = ps, # a vector of allele frequencies at M SNPs
    pgs = pgs.mean, # a vector of PGS point estimates for n individuals
    pgs.est = pgs.est, # a (n x n.samples) matrix of PGS samples
    pgs.rankings = rankings, # a (n x n.samples) matrix of rankings of n.samples PGS estimates for n individuals
    prop.over.t = rank.prop # a vector of proportions of PGS values over t:th quantile for each individual
  )) 
}


