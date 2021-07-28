args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c(
  "no-output",
  "no-figs",
  "no-examples",
  "out-to-file",
  "test"
)
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

if ("out-to-file" %in% args) {
  out_file <- file("../output/code_output.txt", open = "wt")
  sink(out_file)
  sink(out_file, type = "message")
}

test <- "test" %in% args


library(tidyverse)
library(flashier)
library(PMA)
library(ssvd)
library(gtools)
library(RColorBrewer)
library(gridExtra)

source("./sim_fns.R")


# Output:

if (!("no-output" %in% args)) {
  cat("\nGenerating output...\n\n")

  setwd("./make_output/")

  source("./sparse_methods_output.R")
  source("./ebmf_methods_output.R")
  source("./admix_methods_output.R")

  setwd("../")
}


# Figures:

if (!("no-figs" %in% args)) {
  cat("\nGenerating figures...\n\n")

  setwd("./make_figs/")

  source("./tree_diagrams.R")
  source("./tree_svds.R")
  source("./sparse_methods_figs.R")
  source("./ebmf_methods_figs.R")
  source("./admix_methods_figs.R")

  setwd("../")
}


# Examples:

if (!("no-examples" %in% args)) {
  setwd("./examples/")

  source("./examples.R")

  setwd("./")
}


cat("\n\n")
sessionInfo()
