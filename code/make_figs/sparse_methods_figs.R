cat("Generating figures for sparse methods comparisons...\n\n")


plot_loadings <- function(tib) {
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

plot_crossprod_res <- function(cp_res) {
  plot(
    ggplot(cp_res, aes(x = Method, y = Crossprod)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18))
  )
}

plot_crossprod_diff <- function(cp_res) {
  cp_res <- cp_res %>%
    group_by(Trial) %>%
    mutate(Diff = Crossprod - sum(Crossprod * (Method == "svd")))
  plot(
    ggplot(cp_res %>% filter(Method != "svd"), aes(x = Method, y = Diff)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(y = "vs. SVD")
  )
}


LL <- readRDS("../../output/sparse_methods_loadings.rds")
plot_loadings(LL)
ggsave("../../figs/methodcomp_balanced.png", height = 6, width = 5)


cp_res <- readRDS("../../output/sparse_methods_cpres.rds")
plot_crossprod_res(cp_res)
ggsave("../../figs/methodcompboxplot1_balanced.png", height = 6, width = 6)

plot_crossprod_diff(cp_res)
ggsave("../../figs/methodcompboxplot2_balanced.png", height = 6, width = 6)
