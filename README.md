# BayesOptimalTreatment
R data and code for the paper:

"A Bayesian Decision Framework for Optimizing Sequential Combination Antiretroviral Therapy in People with HIV".

## Data

The data analyzed in Section 6 “Application: WIHS Data Analysis” of the paper are from The Women's Interagency HIV Study (WIHS), which is a multisite, longitudinal cohort study of the natural and treated history of women living with HIV and women at-risk for HIV in the United States.
The data are publicly available. However, one need to fill in a request form for access. Full details of the data are available at https://statepi.jhsph.edu/wihs/wordpress/. R data for the simulation study and the WIHS data analysis of the paper are available at the following link: https://drive.google.com/file/d/1svezgKtrlniVX0hgfE8-V2gpUVy3FcQT/view?usp=sharing.

## Code 

The R scripts in the folder “Bayes_Optimal_Treatment_Simulation” are for Section 5 “Simulation Study”, and the R scripts in the folder “Bayes_Optimal_Treatment_WIHS” are for Section 6 “Application: WIHS Data Analysis”.

Libraries and Version Numbers: R version 4.0.0, Rcpp 1.0.7, RcppArmadillo 0.9.880.1.0, MCMCpack 1.4-8, coda 0.19-3, TruncatedNormal 2.2, Matrix 1.2-18, BiasedUrn 1.07, ggplot2 3.3.0, gridExtra 2.3, viridis 0.5.1.

### Instructions for Use

In the folder "Bayes_Optimal_Treatment_Simulation":

* The R script “Simu_Main.R” reproduces Figure 3, Figure 4, Figure 5 in the manuscript, and Figure S1, Figure S2, Table S3 in the supplementary material;

* The R script "Simu_Data_Generate.R" generates the simulated dataset, “MCMC_R_Functions.R” provides R functions for MCMC, “MCMC_Rcpp_Functions.cpp” provides Rcpp functions for MCMC, “MCMC_Main.R” provides the MCMC main function, "SGD_Data.R" processes data for SGD, "Drug_Similarity.R" calculates the subset-tree similarity between different cART regimens, “SGD_R_Functions.R” provides R functions for SGD, “SGD_Rcpp_Functions.cpp” provides Rcpp functions for SGD, and “SGD_Main.R” provides the SGD main function;

* The R data file "Simu_Data_Preprocess.Rdata" saves the preprocessed data from the WIHS dataset for generating simulation truths, "Simu_Truths.Rdata" saves the simulation truths, "Simu_MCMC_Results.Rdata" saves the MCMC posterior samples, and "Simu_SGD_Results.Rdata" saves the SGD results.

In the folder "Bayes_Optimal_Treatment_WIHS":

* The R script “WIHS_Main.R” reproduces Figure 6, Figure 7, Figure 8, Figure 9, and Table 1 in the manuscript;

* The R script "WIHS_Data_Preprocess.R" preprocesses the WIHS dataset for data analysis, “MCMC_R_Functions.R” provides R functions for MCMC, “MCMC_Rcpp_Functions.cpp” provides Rcpp functions for MCMC, “MCMC_Main.R” provides the MCMC main function, "SGD_Data.R" processes data for SGD, "Drug_Similarity.R" calculates the subset-tree similarity between different cART regimens, “SGD_R_Functions.R” provides R functions for SGD, “SGD_Rcpp_Functions.cpp” provides Rcpp functions for SGD, and “SGD_Main.R” provides the SGD main function;

* The R data file "WIHS_Data_Preprocess.Rdata" saves the preprocessed WIHS dataset for data analysis, "WIHS_MCMC_Results.Rdata" saves the MCMC posterior samples, and "WIHS_SGD_Results.Rdata" saves the SGD results.
