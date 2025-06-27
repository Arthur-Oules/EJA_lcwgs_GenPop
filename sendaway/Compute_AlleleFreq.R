# This codes is sent to a cluster because the computations are time sonsuming on
# a regular computer. See EJA_lcwgs_population_structure lines 528-870 for
# conversion to genlight

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