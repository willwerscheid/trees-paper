sim_4pops <- function(pop_sizes,
                      branch_sds,
                      indiv_sd,
                      n_genes = 1000,
                      admix_prop = 0.5,
                      constrain_F = FALSE,
                      seed = 666) {
  set.seed(seed)

  n <- sum(pop_sizes)
  p <- n_genes

  FF <- matrix(rnorm(7 * p, sd = rep(branch_sds, each = p)), ncol = 7)
  if (constrain_F) {
    FF_svd <- svd(FF)
    FF <- FF_svd$u
    FF <- t(t(FF) * branch_sds * sqrt(p))
  }

  LL <- matrix(0, nrow = n, ncol = 7)
  LL[, 1] <- 1
  LL[, 2] <- rep(c(1, 1, 0, 0, admix_prop), times = pop_sizes)
  LL[, 3] <- rep(c(0, 0, 1, 1, 1 - admix_prop), times = pop_sizes)
  LL[, 4] <- rep(c(1, 0, 0, 0, 0), times = pop_sizes)
  LL[, 5] <- rep(c(0, 1, 0, 0, admix_prop), times = pop_sizes)
  LL[, 6] <- rep(c(0, 0, 1, 0, 1 - admix_prop), times = pop_sizes)
  LL[, 7] <- rep(c(0, 0, 0, 1, 0), times = pop_sizes)

  # Only true for trees with no admixture:
  divmat <- matrix(nrow = n, ncol = 4)
  divmat[, 1] <- LL[, 1]
  divmat[, 2] <- LL[, 2] - LL[, 3]
  divmat[, 3] <- LL[, 4] - LL[, 5]
  divmat[, 4] <- LL[, 6] - LL[, 7]

  if(!isTRUE(all.equal(branch_sds[2], branch_sds[-1])) ||
     !isTRUE(all.equal(pop_sizes[1], pop_sizes[-5]))) {
    if (pop_sizes[5] > 0) {
      divmat[(n - pop_sizes[5] + 1):n, 2] <- NA
    }
  }

  E <- matrix(rnorm(n * p, sd = indiv_sd), nrow = n)

  pops <- rep(LETTERS[1:length(pop_sizes)], times = pop_sizes)

  return(list(Y = LL %*% t(FF) + E, LL = LL, FF = FF, divmat = divmat, pops = pops))
}

sim_6pops <- function(pop_sizes,
                      branch_sds,
                      indiv_sd,
                      n_genes = 1000,
                      admix_prop = 0.5,
                      seed = 666) {
  set.seed(seed)

  n <- sum(pop_sizes)
  p <- n_genes

  FF <- matrix(rnorm(11 * p, sd = rep(branch_sds, each = p)), ncol = 11)
  LL <- matrix(0, nrow = n, ncol = 11)
  LL[, 1] <- 1
  LL[, 2] <- rep(c(1, 1, 1, 0, 0, 0, admix_prop), times = pop_sizes)
  LL[, 3] <- rep(c(0, 0, 0, 1, 1, 1, 1 - admix_prop), times = pop_sizes)
  LL[, 4] <- rep(c(1, 0, 0, 0, 0, 0, 0), times = pop_sizes)
  LL[, 5] <- rep(c(0, 0, 0, 1, 0, 0, 1 - admix_prop), times = pop_sizes)
  LL[, 6] <- rep(c(0, 1, 1, 0, 0, 0, admix_prop), times = pop_sizes)
  LL[, 7] <- rep(c(0, 0, 0, 0, 1, 1, 0), times = pop_sizes)
  LL[, 8] <- rep(c(0, 1, 0, 0, 0, 0, 0), times = pop_sizes)
  LL[, 9] <- rep(c(0, 0, 1, 0, 0, 0, admix_prop), times = pop_sizes)
  LL[, 10] <- rep(c(0, 0, 0, 0, 1, 0, 0), times = pop_sizes)
  LL[, 11] <- rep(c(0, 0, 0, 0, 0, 1, 0), times = pop_sizes)


  divmat <- matrix(nrow = n, ncol = 6)
  divmat[, 1] <- LL[, 1]
  divmat[, 2] <- LL[, 2] - LL[, 3]
  divmat[, 3] <- LL[, 4] - LL[, 6]
  divmat[, 4] <- LL[, 5] - LL[, 7]
  divmat[, 5] <- LL[, 8] - LL[, 9]
  divmat[, 6] <- LL[, 10] - LL[, 11]

  E <- matrix(rnorm(n * p, sd = indiv_sd), nrow = n)

  pops <- rep(LETTERS[1:length(pop_sizes)], times = pop_sizes)

  return(list(Y = LL %*% t(FF) + E, LL = LL, FF = FF, divmat = divmat, pops = pops))
}

sim_star <- function(pop_sizes,
                     branch_sds,
                     indiv_sd,
                     n_genes = 1000,
                     admix_prop = 0.5,
                     seed = 666) {
  set.seed(seed)

  n <- sum(pop_sizes)
  p <- n_genes
  K <- length(pop_sizes)

  FF <- matrix(rnorm(p * K, sd = rep(branch_sds, each = p)), ncol = K)
  LL <- matrix(0, nrow = n, ncol = K)

  for (k in 1:K) {
    vec <- rep(0, K)
    vec[k] <- 1
    LL[, k] <- rep(vec, times = pop_sizes)
  }

  divmat <- LL

  E <- matrix(rnorm(n * p, sd = indiv_sd), nrow = n)

  pops <- rep(LETTERS[1:length(pop_sizes)], times = pop_sizes)

  return(list(Y = LL %*% t(FF) + E, LL = LL, FF = FF, divmat = divmat, pops = pops))
}
