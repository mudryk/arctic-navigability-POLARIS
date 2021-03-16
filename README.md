# arctic-navigability-POLARIS
Calculates local grid cell navigability from simulated sea ice thickness (model output). 4 composite "routes" also considered

The scripts in this directory can be used to reproduce figs 2-4 and supplementary figs 1-2 

Procedure:
0. Obtain and process NCAR hi_d daily sea ice thickness output, 1960-2100, 40 realizations. 
Assumed format is 1 file per year of daily data/blocks of 10 realizations (Nr=40). 
Subroutine also provided to read Nr=2 sample data: 1 file per year of daily data/block of 2 realizations
Sample data linked for 1 year of 10 realizations and a 2-realization sample for 1960-2100. 

1. run calculate_output.pro  (e.g. from IDL command line)
computes first/last day of transit and season length for Nr=40 or Nr=2 data as specified
./data/out/: d0d1Ls_WL, d0d1Ls_SR, Nad (data used for figs 2-4)

2. run calculate_d0d1Ls_obs.pro
computes observed values over shipping routes
./data/out/: d0d1Ls_obs (used for supply fig 2 along with d0d1Ls_SR)

3. run plotting routines, eps figures produced in ./figs/
a. plot_d0d1Ls_maps.pro (fig 2 using d0d1Ls_WL)
b. plot_d0d1Ls_SR.pro   (fig 3 using d0d1Ls_SR)
c. plot_Nad.pro         (fig 4 using Nad)
d. plot_bias.pro        (SM-fig 1 using d0d1Ls_obs, d0d1Ls_SR)
