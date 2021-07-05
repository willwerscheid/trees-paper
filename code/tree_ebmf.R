library(flashier)
library(magrittr)

# Initialize an all-ones factor.
ones <- matrix(1, nrow = nrow(dat), ncol = 1)
ls.soln <- t(solve(crossprod(ones), crossprod(ones, dat)))

# Add the all-ones factor and a second factor.
fl <- flash.init(dat) %>%
  flash.init.factors(list(ones, ls.soln)) %>%
  flash.fix.loadings(kset = 1, mode = 1) %>%
  flash.backfit() %>%
  flash.add.greedy(
    Kmax = 1, prior.family = c(prior.point.laplace(), prior.normal())
  )

# Start with the second factor.
current_k <- 2
K <- fl$n.factors

while(current_k <= K && K < Kmax) {
  # Partition individuals according to the sign of the current factor's
  #   loadings.
  splus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] > 0), ncol = 1)
  sminus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] < 0), ncol = 1)

  # Add two new factors, one for the positive loadings and one for the
  #   negative. Fix everything else at zero.
  if (sum(splus) > 0 && sum(sminus) > 0) {
    ls.soln.plus <- t(solve(crossprod(splus),
                            crossprod(splus, dat - fitted(fl))))
    ls.soln.minus <- t(solve(crossprod(sminus),
                             crossprod(sminus, dat - fitted(fl))))
    EF <- list(cbind(splus, sminus),
               cbind(ls.soln.plus, ls.soln.minus))

    fl <- fl %>%
      flash.init.factors(EF) %>%
      flash.fix.loadings(
        kset = K + 1:2, mode = 1L, is.fixed = (EF[[1]] == 0)
      ) %>%
      flash.backfit(kset = K + 1:2) %>%
      flash.nullcheck(remove = TRUE)
  }

  # Move on to the next factor.
  current_k <- current_k + 1
  K <- fl$n.factors
}

# Unfix the loadings and backfit.
fl <- fl %>%
  flash.fix.loadings(kset = 1:K, mode = 1L, is.fixed = FALSE) %>%
  flash.backfit()

# EBcovMF:
if (is_ebcovmf) {
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf

  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    fl <- flash.init(dat_minuss2) %>%
      flash.init.factors(
        EF = fl$flash.fit$EF, EF2 = fl$flash.fit$EF2
      ) %>%
      flash.backfit()

    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
}
