library(vcfR)
library(adegenet)

EJA_gi <- read.vcfR(here("data", "8-Embiotoca_filtered.vcf.gz")) |>
  vcfR2genind()

popmap_lcwgs <- c(
  rep("bodega_bay",            5), # BB
  rep("big_creek",            13), # BIGC -> BCR
  rep("catalina_island",      11), # CAT
  rep("elkhorn",               6), # ELK
  rep("guadalupe_island",     10), # GUA
  rep("pacific_grove",        13), # HOP -> Pacific Grove PGR 
  rep("isla_san_jeronimo",    11), # ISJ
  rep("laguna_beach",         13), # LB
  rep("point_dume",            5), # PD
  rep("santa_barbara",        13), # SB
  rep("san_clemente_island",  10), # SCL
  rep("la_jolla_san_diego",   13), # SD
  rep("tomales_bay",           1), # TB
  rep("redondo_beach",        10), # RB
  rep("santa_cruz_harbour",   10) # SCH
)

EJA_gi@pop <- factor(popmap_lcwgs)

EJA_gen <- genind2genpop(EJA_gi)

saveRDS(EJA_gen, here("output", "EJA_lcwgs_gen.rds"))
saveRDS(EJA_gi, here("output", "EJA_lcwgs_gi.rds"))
rm(EJA_vcf, EJA_gi)

EJA_chord_distance <- dist.genpop(EJA_gen, method = 2) #Chose method = 2 for Edward's distance

saveRDS(EJA_chord_distance, here("output", "EJA_lcwgs_Edward_distance.rds"))