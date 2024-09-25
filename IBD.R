library(here)
library(dartR)

lcwgs_gl <- readRDS(file = here("data", "lcwgs", "Embiotoca_lcwgs_gl.rds"))

IBD <- gl.ibd(lcwgs_gl, coordinates = "latlon")

IBD |> saveRDS(file = here("output", "Embiotocidae_IBD.rds"))