if (exists("test") && test) {
  ntrials <- 2
} else {
  ntrials <- 20
}
cat("Generating output for admixture methods comparisons:", ntrials, "trials...\n\n")


tree_fit <- function(dat, prior = prior.point.laplace(), Kmax = Inf, min_pve = 1e-6) {
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
      if (any(next_fl$pve[K + 1:2] > min_pve)) {
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

  fl$ebcovmf_s2 <- s2

  return(fl)
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

  # Tree (no backfit).
  ebmf_tree_res <- tree_fit(sim_data$Y, Kmax = K) # %>% flash.backfit()
  all_LL$ebmf_tree <- ebmf_tree_res$loadings.pm[[1]]
  all_lfsr$ebmf_tree <- ebmf_tree_res$loadings.lfsr[[1]]

  # EBcovMF with tree constraints.
  ebcovmf_tree_res <- cov_tree_fit(covmat, Kmax = K)
  ebcovmf_s2 <- ebcovmf_tree_res$ebcovmf_s2
  # ebcovmf_tree_res <- ebcovmf_tree_res %>%
  #   flash.backfit()
  all_LL$ebcovmf_tree <- ebcovmf_tree_res$loadings.pm[[1]]
  all_lfsr$ebcovmf_tree <- ebcovmf_tree_res$loadings.lfsr[[1]]

  # Relaxed tree.
  ebmf_relax_res <- flash.init(sim_data$Y) %>%
    flash.set.verbose(0) %>%
    flash.init.factors(EF = ebmf_tree_res$flash.fit$EF,
                       EF2 = ebmf_tree_res$flash.fit$EF2,
                       prior = c(prior.point.laplace(), prior.normal())) %>%
    flash.backfit()
  all_LL$ebmf_relax <- ebmf_relax_res$loadings.pm[[1]]
  all_lfsr$ebmf_relax <- ebmf_relax_res$loadings.lfsr[[1]]

  # Relaxed EBcovMF.
  ebcovmf_relax_res <- flash.init(covmat) %>%
    flash.set.verbose(0) %>%
    flash.init.factors(EF = ebcovmf_tree_res$flash.fit$EF,
                       EF2 = ebcovmf_tree_res$flash.fit$EF2,
                       prior = prior.point.laplace()) %>%
    flash.backfit()
  all_LL$ebcovmf_relax <- ebcovmf_relax_res$loadings.pm[[1]]
  all_lfsr$ebcovmf_relax <- ebcovmf_relax_res$loadings.lfsr[[1]]

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

    acc_zero <- colSums(next_lfsr * (true_LL == 0), na.rm = TRUE)
    acc_lambda <- colSums((1 - next_lfsr) * (next_LL > 0) * (true_LL > 0), na.rm = TRUE)
    acc_nu <- colSums((1 - next_lfsr) * (next_LL < 0) * (true_LL < 0), na.rm = TRUE)
    acc_lambda_flip <- colSums((1 - next_lfsr) * (next_LL < 0) * (true_LL > 0), na.rm = TRUE)
    acc_nu_flip <- colSums((1 - next_lfsr) * (next_LL > 0) * (true_LL < 0), na.rm = TRUE)

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
  perm_LL <- t(t(perm_LL) / apply(perm_LL, 2, function(x) {
    if (any(x != 0)) {
      max(abs(x))
    } else {
      1
    }
  }))

  return(list(perm_LL = perm_LL, acc = best_acc))
}

accuracy_experiment <- function(ntrials, sim_fn, sim_param_fn, K = 4) {
  acc_tib <- tibble()
  admix_tib <- tibble()
  for (i in 1:ntrials) {
    cat("  Trial", i, "\n")
    sim_par <- sim_param_fn(i)
    sim_par$seed <- i
    sim_data <- do.call(sim_fn, sim_par)
    res <- est_loadings(sim_data, K = K)
    acc <- sapply(res, function(x) sum(x$acc) / sum(!is.na(sim_data$divmat)))
    admix_idx <- which(sim_data$pops == "E")
    admix_acc <- sapply(res, function(x) {
      mean(sign(sim_data$divmat[admix_idx, ]) == sign(x$perm_LL[admix_idx, ]),
           na.rm = TRUE)
    })
    acc_tib <- acc_tib %>%
      bind_rows(acc)
    admix_tib <- admix_tib %>%
      bind_rows(admix_acc)
  }
  cat("\n")

  acc_tib <- acc_tib %>%
    mutate(Trial = row_number()) %>%
    pivot_longer(-Trial, names_to = "Method", values_to = "Accuracy") %>%
    mutate(Method = factor(Method, levels = names(res))) %>%
    add_column(Type = "Underlying Tree")
  admix_tib <- admix_tib %>%
    mutate(Trial = row_number()) %>%
    pivot_longer(-Trial, names_to = "Method", values_to = "Accuracy") %>%
    mutate(Method = factor(Method, levels = names(res))) %>%
    add_column(Type = "Admixed Population")

  return(acc_tib %>% bind_rows(admix_tib))
}


make_loadings_tib <- function(sim_data, res) {
  all_loadings <- do.call(rbind, (lapply(res, `[[`, "perm_LL")))

  colnames(all_loadings) <- paste0("k=", 1:ncol(all_loadings))

  admix_idx <- which(sim_data$pops == "E")

  acc_df <- tibble(Method = names(res),
                   Accuracy = sapply(res, function(x) sum(x$acc) / sum(!is.na(sim_data$divmat))),
                   AdmixAcc = sapply(res, function(x) {
                     mean(sign(sim_data$divmat[admix_idx, ]) == sign(x$perm_LL[admix_idx, ]),
                          na.rm = TRUE)
                   }))

  n <- nrow(sim_data$Y)
  tib <- as_tibble(all_loadings) %>%
    add_column(Idx = rep(1:n, length(res)),
               Pop = rep(sim_data$pops, length(res)),
               Method = rep(names(res), each = n)) %>%
    pivot_longer(-(Idx:Method), names_to = "Factor", values_to = "Loading") %>%
    left_join(acc_df, by = "Method") %>%
    mutate(Method = factor(Method, levels = names(res)))

  return(tib)
}


# Admixture simulation:

sim_param <- function(seed = 666) {
  set.seed(seed)
  list(
    n_genes = 10000,
    pop_sizes = c(rep(40, 4), 20),
    branch_sds = c(10, runif(6, min = 3, max = 6)),
    indiv_sd = 1,
    admix_prop = runif(1, min = 0.5, max = 0.9)
  )
}
sim_data <- do.call(sim_4pops, sim_param())

res <- est_loadings(sim_data)


tib <- make_loadings_tib(sim_data, res)
saveRDS(tib, "../../output/admix_methods_loadings.rds")

acc_res <- accuracy_experiment(ntrials, sim_4pops, sim_param)
levels(acc_res$Method) <- c("tree", "cov_tree", "relax", "cov_relax")
saveRDS(acc_res, "../../output/admix_methods_accres.rds")
