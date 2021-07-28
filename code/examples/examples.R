cat("Running examples...\n\n")


tree_fit <- function(dat, prior = prior.point.laplace(), Kmax = Inf, min_pve = 0, verbose.lvl = 0) {
  ones <- matrix(1, nrow = nrow(dat), ncol = 1)
  ls.soln <- t(solve(crossprod(ones), crossprod(ones, dat)))

  fl <- flash.init(dat) %>%
    flash.set.verbose(verbose.lvl) %>%
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
    if (verbose.lvl > 0) {
      cat("K:", K, "\n")
    }
  }

  fl$loadings.lfsr[[1]][, 1] <- 0
  fl$loadings.lfsr[[1]][is.na(fl$loadings.lfsr[[1]])] <- 1

  return(fl)
}

cov_fit <- function(covmat, prior = prior.point.laplace(), Kmax = 1000, verbose.lvl = 0) {
  fl <- flash.init(covmat) %>%
    flash.set.verbose(verbose.lvl) %>%
    flash.add.greedy(Kmax = Kmax, prior.family = prior)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- flash.init(covmat_minuss2) %>%
      flash.set.verbose(verbose.lvl) %>%
      flash.add.greedy(Kmax = Kmax, prior.family = prior)
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  return(fl)
}

cov_tree_fit <- function(covmat, prior = prior.point.laplace(), Kmax = 1000, verbose.lvl = 0) {
  fl <- tree_fit(covmat, prior, Kmax, verbose.lvl = verbose.lvl)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- tree_fit(covmat_minuss2, prior, Kmax, verbose.lvl = verbose.lvl)
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }

  fl$ebcovmf_s2 <- s2

  return(fl)
}

do_fits <- function(dat, Kmax, tree_Kmax = 100, is_cov = FALSE, verbose.lvl = 0) {
  res <- list()

  svd_res <- svd(dat, nu = Kmax, nv = Kmax)
  res$svd <- list(LL = svd_res$u,
                  pve = svd_res$d[1:ncol(svd_res$u)]^2 / sum(svd_res$d^2))

  if (!is_cov) {
    tree_res <- tree_fit(dat, Kmax = tree_Kmax, verbose.lvl = verbose.lvl)
    tree_res <- tree_res %>%
      flash.nullcheck(remove = TRUE)
    tree_relax_res <- flash.init(dat) %>%
      flash.set.verbose(verbose.lvl) %>%
      flash.init.factors(EF = tree_res$flash.fit$EF,
                         EF2 = tree_res$flash.fit$EF2,
                         prior = prior.point.laplace()) %>%
      flash.backfit() %>%
      flash.nullcheck(remove = TRUE)
    res$tree <- list(LL = tree_res$loadings.pm[[1]],
                     pve = tree_res$pve)
    res$tree_relax <- list(LL = tree_relax_res$loadings.pm[[1]],
                           pve = tree_relax_res$pve)

    covmat <- dat %*% t(dat) / ncol(dat)
  } else {
    covmat <- dat
  }

  covtree_res <- cov_tree_fit(covmat, Kmax = tree_Kmax, verbose.lvl = verbose.lvl)
  covtree_res <- covtree_res %>%
    flash.remove.factors(
      which(covtree_res$pve < sort(covtree_res$pve, decreasing = TRUE)[Kmax])
    )
  covtree_relax_res <- flash.init(covmat) %>%
    flash.set.verbose(verbose.lvl) %>%
    flash.init.factors(EF = covtree_res$flash.fit$EF,
                       EF2 = covtree_res$flash.fit$EF2,
                       prior = prior.point.laplace()) %>%
    flash.backfit() %>%
    flash.nullcheck(remove = TRUE)
  res$covtree <- list(LL = covtree_res$loadings.pm[[1]],
                      pve = covtree_res$pve)
  res$covtree_relax <- list(LL = covtree_relax_res$loadings.pm[[1]],
                            pve = covtree_relax_res$pve)

  return(res)
}


plot_loadings <- function(res, pops, filename = NULL) {
  LL <- lapply(res, `[[`, "LL")

  pve <- lapply(res, `[[`, "pve")
  pve_order <- lapply(pve, order, decreasing = TRUE)

  LL <- mapply(function(x, y) x[, y], LL, pve_order,
               SIMPLIFY = FALSE)
  LL <- lapply(LL, function(x) {
    flip <- 2 * (sum(x > 0) > sum(x < 0)) - 1
    t(t(x) * flip)
  })

  n <- nrow(LL[[1]])
  LL <- lapply(LL, as.vector)
  LL <- mapply(function(x, y) {
    cbind(x,
          rep(1:length(y), each = n),
          rep(y, each = n),
          rep(1:n, length(y)))
  }, LL, pve, SIMPLIFY = FALSE)

  tib <- do.call(rbind, LL)
  colnames(tib) <- c("Loading", "Factor", "PVE", "Idx")

  tib <- as_tibble(tib) %>%
    mutate(Method = rep(names(res), times = sapply(LL, nrow)),
           Pop = pops[Idx],
           Factor = paste0("k = ", Factor)) %>%
    mutate(Method = factor(Method, levels = names(res)))

  plt <- ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(rows = vars(Method), cols = vars(Factor), scales = "free_y") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "lines"))

  if (!is.null(filename)) {
    ggsave(paste0("../../figs/", filename, ".png"), plt, height = 6, width = 8)
  } else {
    plot(plt)
  }
}

plot_wolves <- function(LL, pve, lat, long, filename = NULL) {
  df <- t(t(LL) / apply(LL, 2, function(x) max(abs(x))))
  df <- df[, order(pve, decreasing = TRUE)]
  flip <- apply(df, 2, function(x) 2 * (sum(x == 1) > sum(x == -1)) - 1)
  df <- t(t(df) * flip)
  colnames(df) <- paste0("k = ", 1:ncol(df))
  df <- as_tibble(df)
  df <- df %>% mutate(lat = lat, long = long, idx = row_number())
  df <- df %>% pivot_longer(-(lat:idx), names_to = "factor", values_to = "loading")

  jit <- 2
  buf <- 5
  pt_size <- 3

  world <- map_data("world")

  plt <- ggplot() +
    geom_polygon(data = world,
                 aes(long, lat, group = group),
                 fill = NA,
                 color = "black",
                 size = 0.2) +
    geom_jitter(data = df,
                aes(long, lat, fill = loading),
                width = jit, height = jit, shape = 21, size = pt_size) +
    scale_fill_gradient2(midpoint = 0) +
    coord_map(projection = "simpleconic", lat0 = 45.86036851, lat1 = 78.3747588) +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    xlim(min(df$long) - buf, max(df$long) + buf) +
    ylim(min(df$lat) - buf, max(df$lat) + buf) +
    labs(fill = "loading") +
    facet_wrap(~factor)

  if (!is.null(filename)) {
    ggsave(paste0("../../figs/wolves_", filename, ".png"), plt, height = 6, width = 8)
  } else {
    plot(plt)
  }
}

plot_tgp <- function(LL, pve, superpop, pop, save = FALSE) {
  df <- t(t(LL) / apply(LL, 2, function(x) max(abs(x))))
  df <- df[, order(pve, decreasing = TRUE)]
  flip <- apply(df, 2, function(x) 2 * (sum(x == 1) > sum(x == -1)) - 1)
  df <- t(t(df) * flip)
  colnames(df) <- paste0("k = ", formatC(1:ncol(df), width = 2, format = "d", flag = "0"))
  df <- as_tibble(df)
  df <- df %>% mutate(idx = row_number(), pop = pop, superpop = superpop)
  df <- df %>% pivot_longer(-(idx:superpop), names_to = "factor", values_to = "loading")
  df <- df %>% mutate(pve = rep(sort(pve, decreasing = TRUE), nrow(LL)))
  df <- df %>% mutate(factor = factor(factor))

  levels(df$factor) <- c(
    "Common",
    "Div: OOA",
    "Div: EAS/EUR",
    "Drift: SAS",
    "Drift: AMR (PEL)",
    "Drift: AMR (CLM/MXL)",
    "Div: intra-EAS",
    "Div: intra-EUR",
    "Div: intra-AFR",
    "Drift: EAS",
    "Noise"
  )

  plt <- ggplot(df, aes(x = idx, y = loading, col = superpop)) +
    geom_point(size = 0.5, shape = 21) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    facet_wrap(~factor, ncol = 3)
  if (save) {
    ggsave("../../figs/tgp_allpops.png", height = 8, width = 6.5)
  } else {
    plot(plt)
  }

  df <- df %>%
    group_by(superpop, factor) %>%
    mutate(loading_sd = sd(loading)) %>%
    ungroup() %>%
    filter(loading_sd > .05)

  superpop_plots <- list()
  for (sp in unique(superpop)) {
    sub_df <- df %>%
      filter(superpop == sp) %>%
      mutate(pop = factor(pop, levels = unique(pop)))

    if (nrow(sub_df) > 0) {
      superpop_plots[[sp]] <- ggplot(sub_df, aes(x = idx, y = loading, col = pop)) +
        geom_point(size = 1, shape = 21) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_brewer(palette = "Dark2") +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 10)) +
        theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 8),
              axis.title.y = element_blank()) +
        facet_wrap(~factor) +
        ggtitle(paste("Super-pop.:", sp))
    }
  }

  g1 <- arrangeGrob(
    superpop_plots[["EAS"]],
    superpop_plots[["EUR"]],
    nrow = 1,
    widths = c(5, 4)
  )
  g2 <- arrangeGrob(
    superpop_plots[["AFR"]],
    superpop_plots[["AMR"]],
    g1,
    superpop_plots[["SAS"]],
    ncol = 1,
    heights = c(7, 6.5, 4.25, 4.25)
  )

  if (save) {
    ggsave("../../figs/tgp_subpops.png", g2, height = 10, width = 8)
  } else {
    plot(g2)
  }

}


# 6 population simulation:
cat("  Six-population simulation...\n\n")

sim_data_6pop <- sim_6pops(pop_sizes = c(rep(100, 6), 40),
                           branch_sds = c(20, 10, 6, 4, 4, 4, 2, 2, 2, 1, 1),
                           indiv_sd = 1,
                           admix_prop = 0.7,
                           n_genes = 50000)
sim_data_6pop$pops[sim_data_6pop$pops == "G"] <- "Mix (CD)"
res_6pop <- do_fits(sim_data_6pop$Y, Kmax = 6)

plot_loadings(res_6pop[c(1, 3, 5)], sim_data_6pop$pops, filename = "sim_6pops")


# Star simulation:
cat("  Star simulation...\n\n")

sim_data_star <- sim_star(pop_sizes = rep(100, 4),
                          branch_sds = 2,
                          indiv_sd = 1,
                          n_genes = 50000)
res_star <- do_fits(sim_data_star$Y, Kmax = 4)

plot_loadings(res_star[c(1, 3, 5)], sim_data_star$pops, filename = "sim_star")


# Wolves:
cat("  Wolves dataset...\n\n")

wolvesadmix <- readRDS("../../data/wolvesadmix_covmat.rds")
wolves_coord <- read_delim("../../data/wolvesadmix_coord.tsv", delim = " ", col_names = FALSE)

res_wolves <- do_fits(wolvesadmix, is_cov = TRUE, Kmax = 7, tree_Kmax = 16)

plot_wolves(
  res_wolves$svd$LL,
  pve = res_wolves$svd$pve,
  long = wolves_coord$X1,
  lat = wolves_coord$X2,
  filename = "svd"
)
plot_wolves(
  res_wolves$covtree_relax$LL,
  pve = res_wolves$covtree_relax$pve,
  long = wolves_coord$X1,
  lat = wolves_coord$X2,
  filename = "ebcovmf"
)


# Thousand Genomes Project:
cat("  1000 Genomes Project data...\n\n")

tgp <- readRDS("../../data/tgp_covmat.rds")
tgp_meta <- readRDS("../../data/tgp_meta.rds")

res_tgp <- do_fits(tgp, is_cov = TRUE, Kmax = 11, tree_Kmax = 32)

plot_tgp(
  res_tgp$covtree_relax$LL,
  pve = res_tgp$covtree_relax$pve,
  superpop = tgp_meta$super_pop,
  pop = tgp_meta$pop,
  save = TRUE
)
