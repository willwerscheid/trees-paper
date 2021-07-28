cat("Generating figures for admixture methods comparisons...\n\n")


plot_loadings <- function(tib) {
  tib <- tib %>%
    mutate(Accuracy = paste("Tree Acc:", round(Accuracy, 3)),
           AdmixAcc = paste("Admix Acc:", round(AdmixAcc, 3)))

  plot(
    ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(Method, Accuracy, AdmixAcc), cols = vars(Factor), scales = "free_y") +
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
            axis.title = element_text(size = 18)) +
      facet_wrap(~Type)
  )
}


LL <- readRDS("../../output/admix_methods_loadings.rds")
plot_loadings(LL)
ggsave("../../figs/methodcomp_admix.png", height = 6, width = 6)


acc_res <- readRDS("../../output/admix_methods_accres.rds")
plot_accuracy_res(acc_res)
ggsave("../../figs/methodcompboxplot_admix.png", height = 4, width = 6)
