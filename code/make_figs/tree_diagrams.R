cat("Creating tree diagrams...\n\n")


# 4 population tree, no admixture

tib <- tibble(y = c(12, 12, rep(8, 2), rep(4, 4)),
              yend = c(12, 8, rep(4, 2), rep(0, 4)),
              x = c(0, 0, rep(0, 2), rep(-3, 2), rep(3, 2)),
              xend = c(0, 0, c(-3, 3), c(-5, -1, 1, 5)),
              v_label = c("v[root]", "v[ABCD]", "v[AB]", "v[CD]",
                          "", "", "", ""),
              e_label = c("", "alpha", "beta[AB]", "beta[CD]",
                          "gamma[A]", "gamma[B]", "gamma[C]", "gamma[D]"),
              e_posy = c(12, 10, rep(6, 2), rep(2, 4)),
              e_posx = c(0, 0.4, -2.4, 2.5, -4.6, -1.4, 1.4, 4.6))

plt <- ggplot(tib, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment() +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 6,
            nudge_x = c(0.7, 0.9, -0.6, 0.6, rep(-0.5, 2), rep(0.5, 2))) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 8, parse = TRUE) +
  annotate("text", x = c(-5, -1, 1, 5), y = rep(-0.7, 4),
           label = c("A", "B", "C", "D"), size = 10) +
  theme_void() +
  coord_fixed()

ggsave("../../figs/tree_4pop.png", height = 6, width = 5)


# 4 population tree, with admixture

tib <- tibble(y = c(12, 12, rep(8, 2), rep(4, 4)),
              yend = c(12, 8, rep(4, 2), rep(0, 4)),
              x = c(0, 0, rep(0, 2), rep(-3, 2), rep(3, 2)),
              xend = c(0, 0, c(-3, 3), c(-6, -2, 2, 6)),
              v_label = c("v[root]", "v[ABCD]", "v[AB]", "v[CD]",
                          "", "", "", ""),
              e_label = c("", "alpha", "beta[AB]", "beta[CD]",
                          "gamma[A]", "gamma[B]", "gamma[C]", "gamma[D]"),
              e_posy = c(12, 10, rep(6, 2), rep(2, 4)),
              e_posx = c(0, 0.4, -2.4, 2.5, -5.2, -2, 2, 5.2))

plt <- ggplot(tib, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment() +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 6,
            nudge_x = c(0.7, 0.9, -0.6, 0.6, rep(-0.5, 2), rep(0.5, 2))) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 8, parse = TRUE) +
  annotate("text", x = c(-6, -2, 2, 6), y = rep(-0.7, 4),
           label = c("A", "B", "C", "D"), size = 10) +
  annotate("segment", x = c(-2, 2), xend = c(-0.2, 0.2),
           y = c(0, 0), yend = c(0, 0),
           linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("point", x = 0, y = 0, size = 3) +
  annotate("text", x = 0, y = -0.7, label = "BC", size = 10) +
  annotate("text", x = c(-1.1, 1.1), y = c(0.4, 0.4), label = c("pi", "1 - pi"),
           size = 6, parse = TRUE) +
  theme_void() +
  coord_fixed()

ggsave("../../figs/tree_4pop_admix.png", height = 6, width = 6)


# 6 population tree, with admixture

tib <- tibble(y = c(12, 12, rep(8, 2), rep(4, 4), rep(2, 4)),
              yend = c(12, 8, rep(4, 2), c(0, 2, 0, 2), rep(0, 4)),
              x = c(0, 0, rep(0, 2), rep(-3, 2), rep(3, 2), c(-2.5, 4.5), c(-2.5, 4.5)),
              xend = c(0, 0, c(-3, 3), c(-6, -2.5, 2, 4.5), c(-2, 6), c(-4, 4)),
              v_label = c("v[root]", "v[ABCDEF]", "v[ABC]", "v[DEF]",
                          "", "", "", "", "v[BC]", "v[DE]", "", ""),
              e_label = c("", "alpha", "beta[ABC]", "beta[DEF]",
                          "gamma[A]", "gamma[BC]", "gamma[D]", "gamma[EF]",
                          "delta[B]", "delta[E]", "delta[C]", "delta[F]"),
              e_posy = c(12, 10, rep(6, 2), c(2, 3, 2, 3), rep(1.2, 4)),
              e_posx = c(0, 0.7, -2.7, 2.75, -5.15, -2.1, 2, 4.5, -3.8, 3.75, -1.75, 5.7), )

plt <- ggplot(tib, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment(size = c(0, 20, 10, 6, rep(4, 3), rep(2, 2), 1, 2, 1) / 2) +
  #geom_point(aes(x = xend, y = yend), size = 3) +
  #geom_text(aes(x = xend, y = yend), parse = TRUE, size = 5,
  #          nudge_x = c(0.6, 1, -0.7, 0.65, rep(-0.5, 2), rep(0.5, 2), c(-0.55, 0.55), rep(0, 2))) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 8, parse = TRUE) +
  annotate("text", x = c(-6, -4, -2, 2, 4, 6), y = rep(-0.7, 6),
           label = c("A", "B", "C", "D", "E", "F"), size = 10) +
  annotate("segment", x = c(-2, 2), xend = c(-0.2, 0.2),
           y = c(0, 0), yend = c(0, 0),
           linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("point", x = 0, y = 0, size = 3) +
  annotate("text", x = 0, y = -0.7, label = "CD", size = 10) +
  annotate("text", x = c(-1.1, 1.1), y = c(0.4, 0.4), label = c("pi", "1 - pi"),
           size = 5, parse = TRUE) +
  theme_void() +
  coord_fixed()
plot(plt)

ggsave("../../figs/tree_6pop_admix.png", height = 6, width = 6)


# 4 population star

tib <- tibble(y = rep(0, 6),
              yend = c(0, 2 * sqrt(2), 2 * sqrt(2) / sqrt(5), -2, -2, 2 * sqrt(2) / sqrt(5)),
              x = rep(0, 6),
              xend = c(0, 0, -4 * sqrt(2) / sqrt(5), -2, 2, 4 * sqrt(2) / sqrt(5)),
              v_label = c("", "v[root]", "", "", "", ""),
              e_label = c("", "alpha", "beta[A]", "beta[B]", "beta[C]", "beta[D]"),
              e_posy = c(0, sqrt(2), sqrt(2) / sqrt(5), -1, -1, sqrt(2) / sqrt(5)),
              e_posx = c(0, 0.2, -1.9, -1.4, 1.5, 2))

plt <- ggplot(tib, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment() +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 6,
            nudge_x = c(0, 0.4, 0, 0, 0, 0)) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 8, parse = TRUE) +
  annotate("text",
           x = c(-4 * sqrt(2) / sqrt(5) - 0.4, -2, 2, 4 * sqrt(2) / sqrt(5) + 0.4),
           y = c(2 * sqrt(2) / sqrt(5), -2.4, -2.4, 2 * sqrt(2) / sqrt(5)),
           label = c("A", "B", "C", "D"), size = 10) +
  theme_void() +
  coord_fixed()

plot(plt)

ggsave("../../figs/star_4pop.png", height = 6, width = 5)


# Figure illustrating the proof of the existence of the divergence factorization

tib1 <- tibble(y = c(12, 12, rep(8, 2), rep(4, 2)),
               yend = c(12, 8, rep(4, 2), rep(0, 2)),
               x = c(0, 0, rep(0, 2), rep(-3, 2)),
               xend = c(0, 0, c(-3, 3), c(-6, 0)),
               v_label = c("", "v[anc]", "v[AB]", "", "", ""),
               e_label = c("", "", "gamma[AB]", "", "gamma[A]", "gamma[B]"),
               e_posy = c(12, 10, rep(6, 2), rep(2, 2)),
               e_posx = c(0, 0.4, -2.4, 2.5, -5.4, -0.8))

p1 <- ggplot(tib1, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment(linetype = c("solid", "dashed", "solid", "dashed", "solid", "solid")) +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 6,
            nudge_x = c(0.7, 0.9, -0.8, 0.6, rep(-0.5, 2))) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 8, parse = TRUE) +
  annotate("text", x = c(-6, 0), y = rep(-0.7, 2),
           label = c("A", "B"), size = 10) +
  theme_void() +
  coord_fixed() +
  ggtitle("")

ggsave("../../figs/div_proof_Tk.png", height = 6, width = 6)

tib2 <- tib1[1:4, ]
tib2$v_label <- c("", "v[anc]", "", "")

p2 <- ggplot(tib2, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment(linetype = c("solid", "dashed", "solid", "dashed")) +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 8,
            nudge_x = c(0.7, 0.8, -0.8, 0.6)) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 10, parse = TRUE) +
  annotate("text", x = -3, y = 3.3, label = c("AB"), size = 10) +
  theme_void() +
  coord_fixed()

ggsave("../../figs/div_proof_Tk-1.png", height = 6, width = 6)

tib3 <- tib2
tib3$e_label <- c("", "", 'gamma[AB] + delta[AB]^"+"', "")
tib3$e_posx <- tib3$e_posx - 0.6

p3 <- ggplot(tib3, aes(x = x, xend = xend, y = y, yend = yend, label = v_label)) +
  geom_segment(linetype = c("solid", "dashed", "solid", "dashed")) +
  geom_point(aes(x = xend, y = yend), size = 3) +
  geom_text(aes(x = xend, y = yend), parse = TRUE, size = 8,
            nudge_x = c(0.7, 0.8, -0.8, 0.6)) +
  geom_text(aes(x = e_posx, y = e_posy, label = e_label), size = 10, parse = TRUE) +
  theme_void() +
  coord_fixed()

ggsave("../../figs/div_proof_tildeTk-1.png", height = 6, width = 6)
