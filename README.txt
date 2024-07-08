README

#############################################################################################################################
# Source code and data for the manuscript "A marginalized zero-inflated negative binomial                                   #
#                    model for spatial data: modeling COVID-19 deaths in Georgia"                                           #
#   Authors     :    Fedelis Mutiso, Hong Li, John L. Pearce, Sara E. Benjamin-Neelon, Noel T. Mueller, and Brian Neelon    #
#############################################################################################################################

 
Comments: 

1) The data subfolder contains the case study, sensitivity analysis, and simulation study data files

2) The results subfolder contains results output from the simulation and case studies

3) ordered_adjmat_ga.txt in the data subfolder contains adjacency matrix for US state of Georgia used in all scripts

4) covid_data_ga.csv in the data subfolder is the dataset for COVID-19 outcomes and other covariates used in the case study and sensitivity analysis

5) beta.inits.csv in the data subfolder contains the initial values for the case study and sensitivity analysis model parameters

6) To reproduce Table 1 and Figures 3 and 4 in the manuscript, run the simulation_study.R script. The run time for 20k iterations with burnin of 10k and a thin of 10 is approximately 1.06 hours on a Windows 11 desktop computer with an Intel i7-9700 processor running at 3 GHZ with 32 GB RAM

7) To reproduce Figure 1, Figure 2, Table 2, Figure 5 and Figure 6, run the case_study.R script. The run time for 100k iterations with burnin of 50k and a thin of 25 is approximately 6 hours

8) Finally, to reproduce Table S2 of the Supplementary material, run the sensitivity_analysis.R script. The run time for 100k iterations with burnin of 50k and a thin of 25 is approximately 6 hours

9) After saving the `bimj_files` folder containing the r scripts and data/results subfolders, user will need to specify the appropriate working directory to this folder in the R programs

10) The results presented in the manuscript were produced using the following versions of R and R packages. Results may vary slightly if using versions other than those listed below


####################
# Simulation Study #
####################

> sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3 tmap_3.3-3         sf_1.0-8           rgdal_1.5-32       sp_1.5-0           tigris_1.6.1       data.table_1.14.2  forcats_0.5.1     
 [9] stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.2   
[17] fields_14.0        viridis_0.6.2      viridisLite_0.4.0  splines2_0.4.5     spam_2.9-0         msm_1.6.9          MCMCpack_1.6-3     MASS_7.3-58       
[25] coda_0.19-4        mvtnorm_1.1-3      BayesLogit_2.1    

loaded via a namespace (and not attached):
 [1] leafem_0.2.0        googledrive_2.0.0   colorspace_2.0-3    ellipsis_0.3.2      class_7.3-20        leaflet_2.1.1       base64enc_0.1-3     fs_1.5.2           
 [9] dichromat_2.0-0.1   rstudioapi_0.13     proxy_0.4-27        farver_2.1.1        MatrixModels_0.5-0  fansi_1.0.3         lubridate_1.8.0     xml2_1.3.3         
[17] codetools_0.2-18    splines_4.2.2       jsonlite_1.8.0      tmaptools_3.1-1     mcmc_0.9-7          broom_1.0.0         dbplyr_2.2.1        png_0.1-7          
[25] compiler_4.2.2      httr_1.4.3          backports_1.4.1     assertthat_0.2.1    Matrix_1.5-1        fastmap_1.1.0       gargle_1.2.0        cli_3.3.0          
[33] s2_1.1.0            htmltools_0.5.2     quantreg_5.94       tools_4.2.2         dotCall64_1.0-1     gtable_0.3.0        glue_1.6.2          wk_0.6.0           
[41] maps_3.4.0          rappdirs_0.3.3      Rcpp_1.0.9          raster_3.5-21       cellranger_1.1.0    vctrs_0.4.1         leafsync_0.1.0      crosstalk_1.2.0    
[49] lwgeom_0.2-9        rvest_1.0.2         lifecycle_1.0.1     XML_3.99-0.11       googlesheets4_1.0.0 terra_1.5-34        scales_1.2.0        hms_1.1.1          
[57] parallel_4.2.2      expm_0.999-6        SparseM_1.81        curl_4.3.2          gridExtra_2.3       stringi_1.7.8       maptools_1.1-4      e1071_1.7-11       
[65] rlang_1.0.3         pkgconfig_2.0.3     lattice_0.20-45     labeling_0.4.2      htmlwidgets_1.5.4   tidyselect_1.1.2    magrittr_2.0.3      R6_2.5.1           
[73] generics_0.1.3      DBI_1.1.3           pillar_1.8.0        haven_2.5.0         foreign_0.8-83      withr_2.5.0         units_0.8-0         stars_0.5-6        
[81] survival_3.4-0      abind_1.4-5         modelr_0.1.8        crayon_1.5.1        uuid_1.1-0          KernSmooth_2.23-20  utf8_1.2.2          tzdb_0.3.0         
[89] grid_4.2.2          readxl_1.4.0        reprex_2.0.1        digest_0.6.29       classInt_0.4-7      munsell_0.5.0  




##############################
# Case study and Sensitivity #
##############################

> sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3 tmap_3.3-3         sf_1.0-8           rgdal_1.5-32       sp_1.5-0           tigris_1.6.1       data.table_1.14.2  pscl_1.5.5        
 [9] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
[17] tidyverse_1.3.2    splines2_0.4.5     spam_2.9-0         msm_1.6.9          MCMCpack_1.6-3     MASS_7.3-58        coda_0.19-4        mvtnorm_1.1-3     
[25] BayesLogit_2.1    

loaded via a namespace (and not attached):
 [1] googledrive_2.0.0   leafem_0.2.0        colorspace_2.0-3    ellipsis_0.3.2      class_7.3-20        leaflet_2.1.1       base64enc_0.1-3     fs_1.5.2           
 [9] dichromat_2.0-0.1   rstudioapi_0.13     proxy_0.4-27        farver_2.1.1        MatrixModels_0.5-0  bit64_4.0.5         fansi_1.0.3         lubridate_1.8.0    
[17] xml2_1.3.3          codetools_0.2-18    splines_4.2.2       jsonlite_1.8.0      tmaptools_3.1-1     mcmc_0.9-7          broom_1.0.0         dbplyr_2.2.1       
[25] png_0.1-7           compiler_4.2.2      httr_1.4.3          backports_1.4.1     assertthat_0.2.1    Matrix_1.5-1        fastmap_1.1.0       gargle_1.2.0       
[33] cli_3.3.0           s2_1.1.0            htmltools_0.5.2     quantreg_5.94       tools_4.2.2         dotCall64_1.0-1     gtable_0.3.0        glue_1.6.2         
[41] wk_0.6.0            rappdirs_0.3.3      Rcpp_1.0.9          raster_3.5-21       cellranger_1.1.0    vctrs_0.4.1         leafsync_0.1.0      crosstalk_1.2.0    
[49] lwgeom_0.2-9        rvest_1.0.2         lifecycle_1.0.1     XML_3.99-0.11       googlesheets4_1.0.0 terra_1.5-34        scales_1.2.0        vroom_1.5.7        
[57] hms_1.1.1           parallel_4.2.2      expm_0.999-6        SparseM_1.81        curl_4.3.2          stringi_1.7.8       maptools_1.1-4      e1071_1.7-11       
[65] rlang_1.0.3         pkgconfig_2.0.3     lattice_0.20-45     labeling_0.4.2      htmlwidgets_1.5.4   bit_4.0.4           tidyselect_1.1.2    magrittr_2.0.3     
[73] R6_2.5.1            generics_0.1.3      DBI_1.1.3           pillar_1.8.0        haven_2.5.0         foreign_0.8-83      withr_2.5.0         units_0.8-0        
[81] stars_0.5-6         survival_3.4-0      abind_1.4-5         modelr_0.1.8        crayon_1.5.1        uuid_1.1-0          KernSmooth_2.23-20  utf8_1.2.2         
[89] tzdb_0.3.0          grid_4.2.2          readxl_1.4.0        reprex_2.0.1        digest_0.6.29       classInt_0.4-7      munsell_0.5.0       viridisLite_0.4.0  
> 