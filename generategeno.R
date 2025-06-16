
# Function for testing individual PGS variability with different sample sizes, heritability, # of SNPs

PGS.uncertainty <- function(N, M, h2, S = -1){
  # N: # of individuals
  # M: # of SNPs
  # h2: heritability
  # S: relative contributions of common vs rare variants
  
  ps = runif(M, 0, 1) #  allele 1 frequencies for M SNPs
  var.neg = h2/((2*ps*(1-ps))^S*M) # variance proportional to negative selection coefficient, scaled to match heritability h2 and M SNPs
  betas.neg = rnorm(M, mean = 0, sd = sqrt(var.neg)) # draw betas with variance var.neg
  
  beta.neg.est = rnorm(M, betas.neg, sqrt(1/N *(1 - h2/M))) # add noise proportional to sample size, formula from paper 3
  
  X <- matrix(rbinom(N * M, size = 2, prob = rep(ps, each = N)), nrow = N, ncol = M) # generate genotypes assuming allele frequencies follow HWE
  idx.filter = apply(X, 2, sd) > 0
  X = X[, idx.filter]
  X = scale(X, center = TRUE, scale = TRUE)
  beta.neg.est = beta.neg.est[idx.filter]
  
  pgs = X %*% beta.neg.est # compute PGS accounting for negative selection
  
  sigma.est = 1/N *(1 - h2/M) 
  pgs.var = sigma.est * rowSums(X^2) # variance in individual PGS estimates due to noise in the effect size estimates
  
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
  pgs.var = apply(pgs.samples, 2, var)
  
  return(list(pgs = pgs, variance = pgs.var))
}

ordered.estimates <- function(pgs.est, pgs.var, t, rho){
# pgs.est: list of PGS point estimates for each individual
# pgs.var: list of variance of the PGS point estimate for each individual
# t: above PGS quantile to be considered
# rho: level of confidence intervals
  
  N = length(pgs.est)
  ind = 1:N # indices of PGS point estimates
  ind = ind[order(pgs.est[ind], decreasing = FALSE)] # order indices from smallest PGS to largest
  q.ind = ceiling(t * N) # cutoff index
  
  above.ind = ind[(q.ind + 1):N] # top (1-t)*100 % individuals
  below.ind = ind[1:q.ind] # individuals below q.ind
  
  above.pgs = pgs.est[above.ind] # highest PGS estimates
  above.pgs.var = pgs.var[above.ind] # variances of highest PGS estimates
  below.pgs = pgs.est[below.ind] # lowest t*100 % PGS estimates
  below.pgs.var = pgs.var[below.ind] # variances of lowest PGS estimates
  
  # confidence intervals
  z = qnorm((1 + rho) / 2) # z = 1.96 for 95 % CI
  above.pgs.lower = above.pgs - z*sqrt(above.pgs.var) # lower end
  above.pgs.upper = above.pgs + z*sqrt(above.pgs.var) # upper end
  
  below.pgs.lower = below.pgs - z*sqrt(below.pgs.var) # lower end
  below.pgs.upper = below.pgs + z*sqrt(below.pgs.var) # upper end
  
  
  # proportion of individuals having CIs entirely above the t-quantile cutoff
  cutoff.val = quantile(pgs.est, probs = t) # value at t:th quantile
  ci.above.cutoff = above.pgs.lower > cutoff.val # number of lower CIs above t:th quantile, returns a logical vector
  proportion.certain.above = mean(ci.above.cutoff) # proportion of pgs CIs entirely above t among those pgs estimates above t
  
  # proportion of individuals having CIs entirely below the t-quantile cutoff
  ci.below.cutoff = below.pgs.upper < cutoff.val # number of upper CIs below t:th quantile, returns a logical vector
  proportion.certain.below = mean(ci.below.cutoff) # proportion of pgs CIs entirely below t among those pgs estimates below t
  
  return(list(
    above.pgs = above.pgs, # above t % PGS estimates ...
    above.pgs.var = above.pgs.var, # ... and their corresponding variances
    above.pgs.lower = above.pgs.lower, # lower end of rho level CIs
    above.pgs.upper = above.pgs.upper, # upper end of rho level CIs
    cutoff.val = cutoff.val, # PGS value at t:th quantile
    prop.certain.above = proportion.certain.above, # proportion of PGSs certainly above t ...
    confident.above.indices = above.ind[ci.above.cutoff], # ... and their indices
    below.pgs = below.pgs, # below t % PGS estimates ...
    below.pgs.var = below.pgs.var, # ... and their corresponding variances
    below.pgs.lower = below.pgs.lower, # lower end of rho level CIs
    below.pgs.upper = below.pgs.upper, # upper end of rho level CIs
    prop.certain.below = proportion.certain.below, # proportion of PGSs certainly below t ...
    confident.below.indices = below.ind[ci.below.cutoff] # ... and their indices
  
  ))
  
}

