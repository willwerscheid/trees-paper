cat("Generating figures for EBMF methods comparisons...\n\n")


plot_loadings <- function(tib) {
  plot(
    ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(Method, Accuracy), cols = vars(Factor), scales = "free_y") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.spacing = unit(1, "lines"))
  )
}

plot_accuracy_res <- function(acc_res) {
  plot(
    ggplot(acc_res, aes(x = Method, y = Accuracy)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18))
  )
}

plot_accuracy_diff <- function(acc_res) {
  acc_res <- acc_res %>%
    group_by(Trial) %>%
    mutate(Diff = Accuracy - sum(Accuracy * (Method == "bf")))
  plot(
    ggplot(acc_res %>% filter(Method != "bf"), aes(x = Method, y = Diff)) +
      geom_boxplot() +
      theme_minimal() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 18)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(y = "vs. ebmf_bf")
  )
}


LL <- readRDS("../../output/ebmf_methods_loadings.rds")
plot_loadings(LL)
ggsave("../../figs/methodcomp_unbalanced.png", height = 8, width = 6.5)


acc_res <- readRDS("../../output/ebmf_methods_accres.rds")
plot_accuracy_res(acc_res)
ggsave("../../figs/methodcompboxplot1_unbalanced.png", height = 6, width = 6)

plot_accuracy_diff(acc_res)
ggsave("../../figs/methodcompboxplot2_unbalanced.png", height = 6, width = 6)
