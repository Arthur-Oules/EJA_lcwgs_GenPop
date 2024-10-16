library(here)
library(dartR)
library(tidyverse)

lcwgs_gl <- readRDS(file = here("data", "lcwgs", "Embiotoca_lcwgs_gl.rds"))

lcwgs_gl_islanders <- lcwgs_gl[!lcwgs_gl$pop %in% c("guadalupe_island", "catalina_island", "san_clemente_island"), ] |>
  gl.filter.monomorphs(lcwgs_gl_islanders, verbose = 5) |> 
  gl.recalc.metrics(lcwgs_gl_islanders, verbose = 5)

lcwgs_gl_islanders@other$latlon <- here("data", "low_coverage_sampling.csv") |> 
  read_csv2() |>
  select(Popmap, Latitude, Longitude) |>
  right_join(tibble(Popmap = pop(lcwgs_gl_islanders))) |> 
  rename(lat = Latitude, lon = Longitude) |> 
  select(-c(Popmap))

saveRDS(
  lcwgs_gl_islanders,
  file = here("data", "lcwgs", "Embiotoca_lcwgs_islanders_gl.rds")
)

rm(lcwgs_gl)

IBD <- gl.ibd(lcwgs_gl_islanders, coordinates = "latlon")

IBD |> saveRDS(file = here("output", "Embiotocidae_islanders_IBD.rds"))