if (exists("test") && test) {
  ntrials <- 2
} else {
  ntrials <- 20
}
cat("Generating output for EBMF methods comparisons:", ntrials, "trials...\n\n")


tree_fit <- function(dat, prior = prior.point.laplace(), Kmax = Inf, backfit = FALSE) {
  ones <- matrix(1, nrow = nrow(dat), ncol = 1)
  ls.soln <- t(solve(crossprod(ones), crossprod(ones, dat)))

  fl <- flash.init(dat) %>%
    flash.set.verbose(0) %>%
    flash.init.factors(list(ones, ls.soln)) %>%
    flash.fix.loadings(kset = 1, mode = 1) %>%
    flash.backfit() %>%
    flash.add.greedy(
      Kmax = 1,
      prior.family = c(prior, prior.normal())
    )

  current_k <- 2
  K <- 2

  while(current_k <= K && K < Kmax) {
    splus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] > 0), ncol = 1)
    sminus <- matrix(1L * (fl$loadings.pm[[1]][, current_k] < 0), ncol = 1)

    if (sum(splus) > 0 && sum(sminus) > 0) {
      ls.soln.plus <- t(solve(crossprod(splus), crossprod(splus, dat - fitted(fl))))
      ls.soln.minus <- t(solve(crossprod(sminus), crossprod(sminus, dat - fitted(fl))))
      EF <- list(cbind(splus, sminus), cbind(ls.soln.plus, ls.soln.minus))

      next_fl <- fl %>%
        flash.init.factors(EF) %>%
        flash.fix.loadings(kset = K + 1:2, mode = 1L, is.fixed = (EF[[1]] == 0)) %>%
        flash.backfit(kset = K + 1:2)
      if (backfit) {
        next_fl <- next_fl %>%
          flash.backfit()
      }

      if (any(next_fl$pve[K + 1:2] > 1e-6)) {
        fl <- next_fl
      }
    }

    current_k <- current_k + 1
    K <- fl$n.factors
  }

  fl$loadings.lfsr[[1]][, 1] <- 0
  fl$loadings.lfsr[[1]][is.na(fl$loadings.lfsr[[1]])] <- 1

  return(fl)
}

cov_fit <- function(covmat, prior = prior.point.laplace(), Kmax = 1000) {
  fl <- flash.init(covmat) %>%
    flash.set.verbose(0) %>%
    flash.add.greedy(Kmax = Kmax, prior.family = prior)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- flash.init(covmat_minuss2) %>%
      flash.set.verbose(0) %>%
      flash.add.greedy(Kmax = Kmax, prior.family = prior)
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  return(fl)
}

cov_tree_fit <- function(covmat, prior = prior.point.laplace(), Kmax = 1000) {
  fl <- tree_fit(covmat, prior, Kmax)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- tree_fit(covmat_minuss2, prior, Kmax)
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  return(fl)
}

ebnm_div <- function(x, s, g_init, fix_g, output, admix = FALSE) {
  if (!fix_g) {
    opt_fn <- function(par) {
      lambda <- exp(par[1])
      nu <- exp(par[2])
      if (admix) {
        g <- ashr::unimix(rep(1/4, 4), c(0, -nu, -nu, lambda), c(0, -nu, lambda, lambda))
      } else {
        g <- ashr::unimix(rep(1/3, 3), c(0, -nu, lambda), c(0, -nu, lambda))
      }

      ebnm_res <- ebnm::ebnm_ash(
        x,
        s,
        g_init = g,
        fix_g = FALSE,
        output = "log_likelihood"
      )
      return(-ebnm_res$log_likelihood)
    }
    opt_res <- optim(
      par = c(log(max(c(1, x))), log(max(c(1, -x)))),
      fn = opt_fn,
      method = "L-BFGS-B"
    )

    lambda <- exp(opt_res$par[1])
    nu <- exp(opt_res$par[2])
    if (admix) {
      g_init <- ashr::unimix(rep(1/4, 4), c(0, -nu, -nu, lambda), c(0, -nu, lambda, lambda))
    } else {
      g_init <- ashr::unimix(rep(1/3, 3), c(0, -nu, lambda), c(0, -nu, lambda))
    }
  }

  return(ebnm::ebnm_ash(x, s, g_init = g_init, fix_g = fix_g, output = output))
}


est_loadings <- function(sim_data, K = 4) {
  all_est <- do_fits(sim_data, K = K)
  all_perms <- mapply(
    function(LL, lfsr) perm_est_to_true(LL, lfsr, sim_data$divmat),
    all_est$LL,
    all_est$lfsr,
    SIMPLIFY = FALSE
  )

  return(all_perms)
}

do_fits <- function(sim_data, K) {
  all_LL <- list()
  all_lfsr <- list()

  covmat <- sim_data$Y %*% t(sim_data$Y) / ncol(sim_data$Y)

  # Backfit from SVD.
  ebmf_bf_res <- flash.init(sim_data$Y) %>%
    flash.set.verbose(0) %>%
    flash.init.factors(
      svd(sim_data$Y, nu = K, nv = K),
      prior.family = c(prior.point.laplace(), prior.normal())
    ) %>%
    flash.backfit()
  all_LL$ebmf_bf <- ebmf_bf_res$loadings.pm[[1]]
  all_lfsr$ebmf_bf <- ebmf_bf_res$loadings.lfsr[[1]]

  # Greedy (no backfit).
  ebmf_greedy_res <- flash.init(sim_data$Y) %>%
    flash.set.verbose(0) %>%
    flash.add.greedy(Kmax = K)
  all_LL$ebmf_greedy <- ebmf_greedy_res$loadings.pm[[1]]
  all_lfsr$ebmf_greedy <- ebmf_greedy_res$loadings.lfsr[[1]]

  # Admix priors (no backfit).
  ebmf_admix_res <- flash.init(sim_data$Y) %>%
    flash.set.verbose(0) %>%
    flash.add.greedy(Kmax = K, prior = c(as.prior(ebnm_div, admix = TRUE), prior.normal()))
  all_LL$ebmf_admix <- ebmf_admix_res$loadings.pm[[1]]
  all_lfsr$ebmf_admix <- ebmf_admix_res$loadings.lfsr[[1]]

  # Div priors (no backfit).
  ebmf_div_res <- flash.init(sim_data$Y) %>%
    flash.set.verbose(0) %>%
    flash.add.greedy(Kmax = K, prior = c(as.prior(ebnm_div), prior.normal()))
  all_LL$ebmf_div <- ebmf_div_res$loadings.pm[[1]]
  all_lfsr$ebmf_div <- ebmf_div_res$loadings.lfsr[[1]]

  # EBcovMF
  ebcovmf_res <- cov_fit(covmat, Kmax = K)
  all_LL$ebcovmf <- ebcovmf_res$loadings.pm[[1]]
  all_lfsr$ebcovmf <- ebcovmf_res$loadings.lfsr[[1]]

  # EBcovMF with div priors.
  ebcovmf_div_res <- cov_fit(covmat, prior = as.prior(ebnm_div), Kmax = K)
  all_LL$ebcovmf_div <- ebcovmf_div_res$loadings.pm[[1]]
  all_lfsr$ebcovmf_div <- ebcovmf_div_res$loadings.lfsr[[1]]

  # Tree (no backfit).
  ebmf_tree_res <- tree_fit(sim_data$Y, Kmax = K)
  all_LL$ebmf_tree <- ebmf_tree_res$loadings.pm[[1]]
  all_lfsr$ebmf_tree <- ebmf_tree_res$loadings.lfsr[[1]]

  # EBcovMF with tree constraints
  ebcovmf_tree_res <- cov_tree_fit(covmat, Kmax = K)
  all_LL$ebcovmf_tree <- ebcovmf_tree_res$loadings.pm[[1]]
  all_lfsr$ebcovmf_tree <- ebcovmf_tree_res$loadings.lfsr[[1]]
  return(list(LL = all_LL, lfsr = all_lfsr))
}

perm_est_to_true <- function(est_LL, est_lfsr, true_LL) {
  K <- ncol(true_LL)
  all_perms <- permutations(K, K)
  best_acc <- rep(0, K)

  for (i in 1:nrow(all_perms)) {
    perm <- all_perms[i, ]
    next_LL <- est_LL[, perm]
    next_lfsr <- est_lfsr[, perm]

    acc_zero <- colSums(next_lfsr * (true_LL == 0))
    acc_lambda <- colSums((1 - next_lfsr) * (next_LL > 0) * (true_LL > 0))
    acc_nu <- colSums((1 - next_lfsr) * (next_LL < 0) * (true_LL < 0))
    acc_lambda_flip <- colSums((1 - next_lfsr) * (next_LL < 0) * (true_LL > 0))
    acc_nu_flip <- colSums((1 - next_lfsr) * (next_LL > 0) * (true_LL < 0))

    col_accs <- acc_zero + acc_lambda + acc_nu
    col_accs_flip <- acc_zero + acc_lambda_flip + acc_nu_flip

    flip <- (col_accs_flip > col_accs)
    next_acc <- pmax(col_accs, col_accs_flip)

    if (sum(next_acc) > sum(best_acc)) {
      best_acc <- next_acc
      perm_LL <- t(t(next_LL) * (1 - 2L * flip))
    }
  }

  # Normalize loadings.
  perm_LL <- t(t(perm_LL) / apply(perm_LL, 2, function(x) max(abs(x))))

  return(list(perm_LL = perm_LL, acc = best_acc))
}

accuracy_experiment <- function(ntrials, sim_fn, sim_param_fn, K = 4) {
  tib <- tibble()
  for (i in 1:ntrials) {
    cat("  Trial", i, "\n")
    sim_par <- sim_param_fn(i)
    sim_par$seed <- i
    sim_data <- do.call(sim_fn, sim_par)
    res <- est_loadings(sim_data, K = K)
    acc <- sapply(res, function(x) sum(x$acc) / length(res[[1]]$perm_LL))
    tib <- tib %>%
      bind_rows(acc)
  }
  cat("\n")

  tib <- tib %>%
    mutate(Trial = row_number()) %>%
    pivot_longer(-Trial, names_to = "Method", values_to = "Accuracy") %>%
    mutate(Method = factor(Method, levels = names(res)))

  return(tib)
}


make_loadings_tib <- function(sim_data, res) {
  all_loadings <- do.call(rbind, (lapply(res, `[[`, "perm_LL")))

  colnames(all_loadings) <- paste0("k=", 1:ncol(all_loadings))

  acc_df <- tibble(Method = names(res),
                   Accuracy = sapply(res, function(x) sum(x$acc) / length(res[[1]]$perm_LL)))

  n <- nrow(sim_data$Y)
  tib <- as_tibble(all_loadings) %>%
    add_column(Idx = rep(1:n, length(res)),
               Pop = rep(sim_data$pops, length(res)),
               Method = rep(names(res), each = n)) %>%
    pivot_longer(-(Idx:Method), names_to = "Factor", values_to = "Loading") %>%
    left_join(acc_df, by = "Method") %>%
    mutate(Method = factor(Method, levels = names(res)),
           Accuracy = round(Accuracy, 3))

  return(tib)
}


# Unbalanced tree simulation:

sim_param <- function(seed = 666) {
  set.seed(seed)
  list(
    n_genes = 10000,
    pop_sizes = c(sample(20:80, size = 4), 0),
    branch_sds = c(10, runif(6, min = 1, max = 6)),
    indiv_sd = 1
  )
}
sim_data <- do.call(sim_4pops, sim_param(seed = 777))

res <- est_loadings(sim_data)


tib <- make_loadings_tib(sim_data, res)
saveRDS(tib, "../../output/ebmf_methods_loadings.rds")

acc_res <- accuracy_experiment(ntrials, sim_4pops, sim_param)
levels(acc_res$Method) <- c("bf", "greedy", "admix", "div", "cov", "cov_div", "tree", "cov_tree")
saveRDS(acc_res, "../../output/ebmf_methods_accres.rds")
