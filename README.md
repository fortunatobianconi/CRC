# CRC
Conditional Robust Calibration 
(Fundamental code)

#################################
# 1) Files description
#################################


This is the main code folder. It contains the following files:

1. CRC_step1_lotkavolterra.m: script for parameter perturbation, simulation of the Lotka-Volterra model and computation of the distance functions
2. CRC_step1_esynthetic.m: script for parameter perturbation, simulation of the E_Synthetic.m model and computation of the distance functions
3. CRC_step1_MM_model.m: script for parameter perturbation, simulation of the p38MAPK MM model and computation of the distance functions


4. BaseHist_multiplearea_piuNr.m: script for performing the intersection between distance functions and for estimating conditional parameter densities. 
   It runs the following functions: upperLowerSet_Nr.m and intersection.m
5. upperLowerSet_Nr.m: selection of the samples of the parameter vector for which all distance functions have values under or above the 
   corresponding tolerance value chosen
6. intersection.m: intersection between the selected values of distance functions of all observables 
7. BaseHist_MIRI.m: script for the computation of the Moment Independent Robustness Indicator (MIRI)


The Models folder contains the following ODE models:
1. LotkaVolterra.m: ODEs of the Lotka-Volterra predator-prey model
2. E_Synthetic.m: ODEs implementation of the E_Synthetic model 
3. MM_model.m: ODE model of the p38 MAPK signaling pathway 

The Data folder contains the datasets for calibration in .mat format:
1. data_LotkaVolterra.mat: two observables measured at 8 time points
2. data_E_Synthetic.mat: two observables measured at 11 time points
3. data_MM_model.mat: 16 observables measured at 6 time points
4. exp_data_index.mat: array of the observables indexes in the experimental dataset for the MM model
5. sim_data_index.mat: array of the observables indexes in the ODE MM model 

#################################
# 2) Usage
#################################

1. Run Robustness_time_series.m (CRC_step1.m), setting the following parameters:
* time_points: time points where data are sampled;
* parameters of the model;
* LBpi and UBpi: lower and upper boundaries of the perturbed parameter space;
* Nr and NSample: number of realizations and number of parameter samples for each realization;
* model: pointer to the model to calibrate;
* Nobs: number of measured output variables.
* definition of the chosen distance function

This script creates two .mat files where all the results are stored: one with the ODE model simulations 
and one with the computed distance functions. 

2. Run CRC_step2.m, setting the following parameters:
* protein_name: names of the proteins selected for calibration
* lowThr or/and highThr: array of chosen thresholds for each distance function.

To perform model calibration, lowThr has to be defined. To run robustness analysis, we also run CRC
using highThr, i.e. we repeat the entire calibration process but in order to identify those parameters
whose corresponding distance functions are larger than highThr. Then, run BaseHist_MIRI.m for generation 
of Moment Independent Robustness Indicator (MIRI). 

This script calls the following functions: upperLowerSet_Nr.m and intersection.m.
It computes: xbin_p and ks_p, cell arrays for storing conditional parameter densities and 
CloudCondL, i.e. modes of the parameters density functions conditioned to the region of interest. 

3. (Optional) Run CRA_MIRI.m: if both the calibration toward low values and high values of distance functions 
have been performed, MIRI can be computed, after loading results obtained in the previous steps. 
This script returns: MIRIT, i.e. array containing, for each realization, MIRI indexes for each parameter and
figures of the parameter density functions, conditioned to both opposite regions.
           

#################################
# 3) OS and Matlab version
#################################

This code has been tested on Ubuntu 16.04 LTS (64bit) using Matlab R2018a.
