# A Comparison of Bayesian Methods for the Stryker Family of Vehicles
This repository contains R code for implementing Bayesian methods for estimating the reliability for the Stryker family of vehicles using data from developmental testing (DT) and operational testing (OT).

The main code is in the following R notebooks:
* `stryker_otonly.rmd`: OT Data with Weakly-informative Priors (Method 1)
* `stryker_dw.rmd`: DT/OT Fixed Downweighting (Methods 2 & 3). The choice of prior can be set by setting the `pr.type` variable to `refpr` (Weakly Informative, Method 2) or `infpr` (Informative, Method 3).
* `stryker_infpr_npp.rmd`: DT/OT with Normalized Power Priors (Method 4) 

For completeness, the results of our runs are provided in `results/` and associated plots are provided in `plots/` (generated using `plots.rmd`).
