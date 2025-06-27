# This codes is sent to a cluster because the computations are time sonsuming on
# a regular computer. See /sendaway/Compute_Edward_distance.R for conversion to
# genpop. This file needs a folder /data/ and /output/ at root

library(tidyverse)
library(adegenet)

here("data", "EJA_lcwgs_gen.rds") |>
  readRDS() |>
  makefreq() |> # Generate allele frequencies
  t() |>
  as_tibble(rownames = "POSITION") |>
  separate_wider_delim(POSITION, delim = "_", names = c("CHROM", "TBD", "POSALLELE")) |>
  separate_wider_delim(POSALLELE, delim = ".", names = c("POS", "ALLELE")) |>
  select(-c(TBD)) |>
  mutate(CHROM = paste0(CHROM, ".1")) |>
  saveRDS(here("output", "EJA_allele_frequencies_tb.rds"))