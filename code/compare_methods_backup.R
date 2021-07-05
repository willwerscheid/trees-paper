library(tidyverse)
library(flashier)
library(PMA)
library(ssvd)
library(gtools)

source("./code/sim_fns.R")

do_fits <- function(sim_data, K = 4, fns) {
  res <- list()

  svd_res <- svd(sim_data$Y, nu = K, nv = K)
  if ("svd" %in% fns) {
    res <- list(svd = svd_res$u)
  }

  if ("varimax" %in% fns) {
    varimax_loadings <- loadings(varimax(svd_res$u))
    class(varimax_loadings) <- "matrix"
    res <- c(res, list(varimax = varimax_loadings))
  }

  if ("ssvd" %in% fns) {
    ssvd_res <- ssvd.iter.thresh(sim_data$Y, method = "theory", gamma.v = 0,
                                 u.old = svd_res$u, v.old = svd_res$v, r = K)
    res <- c(res, list(ssvd = ssvd_res$u %*% diag(ssvd_res$d)))
  }

  if ("pmd" %in% fns) {
    pmd_res <- PMD(sim_data$Y, K = K, trace = FALSE)
    res <- c(res, list(pmd = pmd_res$u %*% diag(pmd_res$d)))
  }

  if ("pmd_cv" %in% fns) {
    pmdcv_res <- PMD.cv(
      sim_data$Y,
      type = "ordered",
      lambda = 0,
      sumabsus = 2^seq(0, log2(sqrt(nrow(sim_data$Y))), length.out = 10),
      v = svd_res$v[, 1],
      trace = FALSE,
      center = FALSE
    )
    pmdcv_out <- PMD(
      sim_data$Y,
      type = "ordered",
      sumabsu = pmdcv_res$bestsumabsu,
      sumabsv = 0,
      lambda = 0,
      K = 1,
      v = pmdcv_res$v.init,
      trace = FALSE,
      center = FALSE
    )
    pmdcv_loadings <- pmdcv_out$u * pmdcv_out$d
    pmdcv_resid <- sim_data$Y - (pmdcv_out$u * pmdcv_out$d) %*% t(pmdcv_out$v)

    for (k in 2:K) {
      pmdcv_res <- PMD.cv(
        sim_data$Y,
        type = "ordered",
        lambda = 0,
        trace = FALSE,
        center = FALSE
      )
      pmdcv_out <- PMD(
        pmdcv_resid,
        type = "ordered",
        sumabsu = pmdcv_res$bestsumabsu,
        sumabsv = 0,
        lambda = 0,
        K = 1,
        v = pmdcv_res$v.init,
        trace = FALSE,
        center = FALSE
      )
      pmdcv_loadings <- cbind(pmdcv_loadings, pmdcv_out$u * pmdcv_out$d)
      pmdcv_resid <- pmdcv_resid - (pmdcv_out$u * pmdcv_out$d) %*% t(pmdcv_out$v)
    }
    res <- c(res, list(pmd_cv = pmdcv_loadings))
  }

  if ("ebmf_pn" %in% fns) {
    ebmf_pn_res <- flash.init(sim_data$Y) %>%
      flash.init.factors(svd_res, prior.family = c(prior.point.normal(), prior.normal())) %>%
      flash.backfit(verbose.lvl = 0)
    res <- c(res, list(ebmf_pn = ebmf_pn_res$loadings.pm[[1]] %*% diag(ebmf_pn_res$loadings.scale)))
  }

  if ("ebmf_pl" %in% fns) {
    ebmf_pl_res <- flash.init(sim_data$Y) %>%
      flash.init.factors(svd_res, prior.family = c(prior.point.laplace(), prior.normal())) %>%
      flash.backfit(verbose.lvl = 0)
    res <- c(res, list(ebmf_pl = ebmf_pl_res$loadings.pm[[1]] %*% diag(ebmf_pl_res$loadings.scale)))
  }

  if ("ebmf_longrun" %in% fns) {
    ebmf_pl_res <- flash.init(sim_data$Y) %>%
      flash.init.factors(svd_res, prior.family = c(prior.point.laplace(), prior.normal())) %>%
      flash.backfit(verbose.lvl = 3, tol = 1e-8, maxiter = 10000)
    res <- c(res, list(ebmf_longrun = ebmf_pl_res$loadings.pm[[1]] %*% diag(ebmf_pl_res$loadings.scale)))
  }

  if ("ebmf_greedy" %in% fns) {
    ebmf_greedy_res <- flash.init(sim_data$Y) %>%
      flash.set.verbose(0) %>%
      flash.add.greedy(Kmax = K, prior.family = c(prior.point.laplace(), prior.normal()), tol = 1e-8)
    res <- c(res, list(ebmf_greedy = ebmf_greedy_res$loadings.pm[[1]] %*% diag(ebmf_greedy_res$loadings.scale)))
  }

  if ("ebmf_cov" %in% fns) {
    covmat <- sim_data$Y %*% t(sim_data$Y) / ncol(sim_data$Y)

    ebmf_cov_res <- flash.init(covmat) %>%
      flash.set.verbose(0) %>%
      flash.add.greedy(Kmax = K, prior.family = prior.point.laplace()) %>%
      flash.backfit(verbose.lvl = 3, maxiter = 20)
    s2 <- max(0, mean(diag(covmat) - diag(fitted(ebmf_cov_res))))
    while(abs(s2_diff - 1) > 1e-4) {
      covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
      ebmf_cov_res <- flash.init(covmat_minuss2) %>%
        flash.set.verbose(0) %>%
        flash.add.greedy(Kmax = K, prior.family = prior.point.laplace()) %>%
        flash.backfit(verbose.lvl = 3, maxiter = 20)
      old_s2 <- s2
      s2 <- max(0, mean(diag(covmat) - diag(fitted(ebmf_cov_res))))
      s2_diff <- s2 / old_s2
    }
    res <- c(res, list(ebmf_cov = ebmf_cov_res$loadings.pm[[1]] %*% diag(ebmf_cov_res$loadings.scale)))
  }

  return(res)
}

perm_est_to_true <- function(true_LL, est_LL) {
  K <- ncol(true_LL)
  all_perms <- permutations(K, K)
  perm_crossprods <- apply(all_perms, 1, function(perm) {
    LL_scale <- sqrt(colSums(true_LL^2) * colSums(est_LL[, perm]^2))
    LL_crossprod <- diag(crossprod(true_LL, est_LL[, perm])) / LL_scale
    return(LL_crossprod)
  })
  which_perm <- which.max(colSums(abs(perm_crossprods)))

  best_perm <- all_perms[which_perm, ]
  flips <- sign(perm_crossprods[, which_perm])
  best_crossprod <- mean(abs(perm_crossprods[, which_perm]))

  perm_LL <- est_LL[, best_perm]
  perm_LL <- t(t(perm_LL) * flips)

  return(list(perm_LL = perm_LL, crossprod = best_crossprod))
}

est_loadings <- function(sim_data, K = 4, fns) {
  all_est <- do_fits(sim_data, K = K, fns)
  all_perms <- lapply(all_est, function(LL) perm_est_to_true(sim_data$divmat, LL))

  return(all_perms)
}

make_tib <- function(sim_data, res) {
  all_loadings <- do.call(rbind, (lapply(res, `[[`, "perm_LL")))

  colnames(all_loadings) <- paste0("k=", 1:ncol(all_loadings))

  crossprod_df <- tibble(Method = names(res),
                         Crossprod = sapply(res, `[[`, "crossprod"))

  n <- nrow(sim_data$Y)
  tib <- as_tibble(all_loadings) %>%
    add_column(Idx = rep(1:n, length(res)),
               Pop = rep(sim_data$pops, length(res)),
               Method = rep(names(res), each = n)) %>%
    pivot_longer(-(Idx:Method), names_to = "Factor", values_to = "Loading") %>%
    left_join(crossprod_df, by = "Method") %>%
    mutate(Method = factor(Method, levels = names(res)),
           Crossprod = round(Crossprod, 3))

  return(tib)
}

do_plot <- function(tib) {
  plot(
    ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(Method, Crossprod), cols = vars(Factor), scales = "free_y") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.spacing = unit(1, "lines"))
  )
}

crossprod_experiment <- function(ntrials, sim_fn, sim_param, K = 4, fns) {
  tib <- tibble()
  for (i in 1:ntrials) {
    cat(i, " ")
    sim_param$seed <- i
    sim_data <- do.call(sim_fn, sim_param)
    res <- est_loadings(sim_data, K = K, fns)
    crossprods <- sapply(res, `[[`, "crossprod")
    tib <- tib %>%
      bind_rows(crossprods)
  }

  tib <- tib %>%
    mutate(Trial = row_number()) %>%
    pivot_longer(-Trial, names_to = "Method", values_to = "Crossprod") %>%
    mutate(Method = factor(Method, levels = names(res)))

  return(tib)
}


# Balanced tree:

sim_param <- list(
  pop_sizes = c(rep(40, 4), 0),
  branch_sds = rep(2, 7),
  indiv_sd = 1
)

sim_data <- do.call(sim_4pops, sim_param)

res <- est_loadings(
  sim_data,
  fns = c("svd", "varimax", "ssvd", "pmd", "pmd_cv", "ebmf_pn", "ebmf_pl")
)
tib <- make_tib(sim_data, res)
do_plot(tib)

ggsave("./figs/methodcomp_balanced.png", height = 6, width = 5)

cp_res <- crossprod_experiment(
  20,
  sim_4pops,
  sim_param,
  fns = c("svd", "varimax", "ssvd", "pmd", "pmd_cv", "ebmf_pn", "ebmf_pl")
)
ggplot(cp_res, aes(x = Method, y = Crossprod)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

ggsave("./figs/methodcompboxplot1_balanced.png", height = 6, width = 6)

cp_res <- cp_res %>%
  group_by(Trial) %>%
  mutate(Diff = Crossprod - sum(Crossprod * (Method == "svd")))
ggplot(cp_res %>% filter(Method != "svd"), aes(x = Method, y = Diff)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "vs. SVD")

ggsave("./figs/methodcompboxplot2_balanced.png", height = 6, width = 6)

sim_param <- list(
  pop_sizes = c(rep(40, 4), 0),
  branch_sds = rep(2, 7),
  indiv_sd = 1,
  seed = 555
)

sim_data <- do.call(sim_4pops, sim_param)

res <- est_loadings(
  sim_data,
  fns = c("svd", "ebmf_pl", "ebmf_longrun", "ebmf_greedy")
)
tib <- make_tib(sim_data, res)
do_plot(tib)

cp_res <- crossprod_experiment(
  20,
  sim_4pops,
  sim_param,
  fns = c("svd", "ebmf_pl", "ebmf_greedy")
)
ggplot(cp_res, aes(x = Method, y = Crossprod)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18))

# Decreasing branch lengths:

sim_param <- list(
  pop_sizes = c(rep(40, 4), 0),
  branch_sds = c(1, rep(0.5, 2), rep(0.2, 4)),
  indiv_sd = 0.1
)

sim_data <- do.call(sim_4pops, sim_param)

res <- est_loadings(sim_data)
tib <- make_tib(sim_data, res)
do_plot(tib)

cp_res <- crossprod_experiment(10, sim_4pops, sim_param)
ggplot(cp_res, aes(x = Method, y = Crossprod)) +
  geom_boxplot()


# Admixture (50-50):

sim_param <- list(
  pop_sizes = c(rep(40, 4), 10),
  branch_sds = c(1, rep(0.5, 2), rep(0.2, 4)),
  indiv_sd = 0.1
)

sim_data <- do.call(sim_4pops, sim_param)

res <- est_loadings(sim_data, K = 4)
tib <- make_tib(sim_data, res)
do_plot(tib)

cp_res <- crossprod_experiment(20, sim_4pops, sim_param)


# Star:

sim_param <- list(
  pop_sizes = rep(30, 6),
  branch_sds = rep(0.5, 6),
  indiv_sd = 0.1
)

sim_data <- do.call(sim_star, sim_param)

res <- est_loadings(sim_data, K = 6)
tib <- make_tib(sim_data, res)
do_plot(tib)


# 6 pops:

sim_param <- list(
  pop_sizes = c(rep(30, 6), 0),
  branch_sds = c(1, rep(0.5, 4), rep(0.25, 6)),
  indiv_sd = 0.1
)

sim_data <- do.call(sim_6pops, sim_param)

res <- est_loadings(sim_data, K = 6)
tib <- make_tib(sim_data, res)
do_plot(tib)
