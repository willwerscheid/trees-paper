library(tidyverse)

source("./code/sim_fns.R")

sim_param <- list(
  list(
    pop_sizes = c(rep(40, 4), 0),
    branch_sds = rep(2, 7),
    indiv_sd = 1
  ),
  list(
    pop_sizes = c(rep(40, 4), 0),
    branch_sds = rep(c(2, sqrt(2.5), 1, 1, 1, 2, 2), 7),
    indiv_sd = 1
  ),
  list(
    pop_sizes = c(10, 50, 30, 70, 0),
    branch_sds = rep(2, 7),
    indiv_sd = 1
  ),
  list(
    pop_sizes = c(rep(40, 4), 0),
    branch_sds = c(2, 2, 4, 1, 4, 3, 2),
    indiv_sd = 1
  )
)

set.seed(666)
sim_data <- lapply(sim_param, function(par) do.call(sim_4pops, par)$Y)
sim_svds <- lapply(sim_data, function(Y) svd(Y, nu = 4, nv = 0)$u)

all_svds <- do.call(rbind, sim_svds)
colnames(all_svds) <- paste("Factor", 1:4)

tree_names <- c("(a) Balanced", "(b) Exact", "(c) Unequal Pops", "(d) Unequal Drift")

tib <- as_tibble(all_svds) %>%
  mutate(Tree = rep(tree_names, times = sapply(sim_param, function(par) sum(par$pop_sizes))),
         Pop = rep(rep(LETTERS[1:5], 4), times = sapply(sim_param, `[[`, "pop_sizes"))) %>%
  group_by(Tree) %>%
  mutate(Idx = row_number()) %>%
  ungroup() %>%
  pivot_longer(`Factor 1`:`Factor 4`, names_to = "Factor", values_to = "Loading")

ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_grid(rows = vars(Tree), cols = vars(Factor)) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  labs(x = "")

ggsave("./figs/tree_svds.png", height = 6, width = 6)
