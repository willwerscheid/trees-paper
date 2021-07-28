if (exists("test") && test) {
  ntrials <- 2
} else {
  ntrials <- 20
}
cat("Generating output for sparse methods comparisons:", ntrials, "trials...\n\n")


est_loadings <- function(sim_data,
                         K = 4,
                         fns = c("svd", "varimax", "ssvd",
                                 "pmd", "pmd_cv", "ebmf_pn", "ebmf_pl")) {
  all_est <- do_fits(sim_data, K = K, fns)
  all_perms <- lapply(all_est, function(LL) perm_est_to_true(sim_data$divmat, LL))

  return(all_perms)
}

do_fits <- function(sim_data, K, fns) {
  res <- list()

  svd_res <- svd(sim_data$Y, nu = K, nv = K)
  if ("svd" %in% fns) {
    res <- list(svd = svd_res$u)
  }

  if ("varimax" %in% fns) {
    varimax_loadings <- loadings(
      varimax(svd_res$u %*% diag(sqrt(svd_res$d[1:K])), normalize = FALSE)
    )
    class(varimax_loadings) <- "matrix"
    res <- c(res, list(varimax = varimax_loadings))
  }

  if ("ssvd" %in% fns) {
    ssvd_res <- ssvd.iter.thresh(sim_data$Y, method = "theory", gamma.v = 0,
                                 u.old = svd_res$u, v.old = svd_res$v, r = K)
    res <- c(res, list(ssvd = ssvd_res$u))
  }

  if ("pmd" %in% fns) {
    pmd_res <- PMD(sim_data$Y, K = K, trace = FALSE)
    res <- c(res, list(pmd = pmd_res$u))
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
      flash.init.factors(
        svd_res,
        prior.family = c(prior.point.normal(), prior.normal())
      ) %>%
      flash.backfit(verbose.lvl = 0)
    res <- c(res, list(ebmf_pn = ebmf_pn_res$loadings.pm[[1]]))
  }

  if ("ebmf_pl" %in% fns) {
    ebmf_pl_res <- flash.init(sim_data$Y) %>%
      flash.init.factors(
        svd_res,
        prior.family = c(prior.point.laplace(), prior.normal())
      ) %>%
      flash.backfit(verbose.lvl = 0)
    res <- c(res, list(ebmf_pl = ebmf_pl_res$loadings.pm[[1]]))
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
  perm_LL <- t(t(perm_LL) * flips / apply(perm_LL, 2, function(x) max(abs(x))))

  return(list(perm_LL = perm_LL, crossprod = best_crossprod))
}

crossprod_experiment <- function(ntrials,
                                 sim_fn,
                                 sim_param,
                                 K = 4,
                                 fns = c("svd", "varimax", "ssvd",
                                         "pmd", "pmd_cv", "ebmf_pn", "ebmf_pl")) {
  tib <- tibble()
  for (i in 1:ntrials) {
    cat("  Trial", i, "\n")
    sim_param$seed <- i
    sim_data <- do.call(sim_fn, sim_param)
    res <- est_loadings(sim_data, K = K, fns)
    crossprods <- sapply(res, `[[`, "crossprod")
    tib <- tib %>%
      bind_rows(crossprods)
  }
  cat("\n")

  tib <- tib %>%
    mutate(Trial = row_number()) %>%
    pivot_longer(-Trial, names_to = "Method", values_to = "Crossprod") %>%
    mutate(Method = factor(Method, levels = names(res)))

  return(tib)
}


make_loadings_tib <- function(sim_data, res) {
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


# Balanced tree simulation:

sim_param <- list(
  pop_sizes = c(rep(40, 4), 0),
  branch_sds = rep(2, 7),
  indiv_sd = 1
)
sim_data <- do.call(sim_4pops, sim_param)

res <- est_loadings(sim_data)


tib <- make_loadings_tib(sim_data, res)
saveRDS(tib, "../../output/sparse_methods_loadings.rds")

cp_res <- crossprod_experiment(ntrials, sim_4pops, sim_param)
saveRDS(cp_res, "../../output/sparse_methods_cpres.rds")
