# Surfperch_GenPop
Scripts for the Black Surfperch population low coverage genomics project conducted at the [Bernardi lab](https://bernardi.eeb.ucsc.edu/) at the University of California, Santa Cruz as part of the [California Conservation Genomics Program](https://www.ccgproject.org/).

See the main article for full details - [PH].

All necessary data to reproduce results can be downloaded on Dryad at <http://doi.org/10.5061/dryad.nk98sf85h>

The repo structure is as such:

EJA_lcwgs_GenPop  
├─ BayeScan2.1/  
|  └─ R functions/  
|     └─ plot_R.R ----------------- Plotting function from [BayeScan2.1](https://cmpg.unibe.ch/software/BayeScan/)  
├─ functions/  
|  ├─ plot_functions.R ------------ Contains helper function for plotting purpose  
|  └─ tidy_functions.R ------------ Contains helper function for analyses purpose  
├─ sendaway.R  
|  ├─ Compute_AlleleFreq.R -------- Contains code sent to a acluster to compute allele frequencies  
|  └─ Compute_Edward_distance ----- Contains code sent to a acluster to compute Cavalli−Sforza and Edward Chord Distance.  
├─ EJA_lcwgs_Extra_Figures.qmd ---- Contains code to generate extra analyses figures  
├─ EJA_lcwgs_Figures.qmd ---------- Contains code that generated all figures in main article  
├─ EJA_lcwgs_outliers.qmd --------- Contains analyses to find outliers loci using [pcadapt 4.4.0](https://cloud.r-project.org/web/packages/pcadapt/index.html), [OutFLANK 0.2](https://github.com/whitlock/OutFLANK) and [BayeScan2.1](https://cmpg.unibe.ch/software/BayeScan/)  
└─ EJA_lcwgs_population_structure - Contains analyses related to isolation by distance and ancestry

If you need any help with the scripts contact me here -> arthur.oules\[at]ibcp\[dot]fr
