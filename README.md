# CRC
Conditional Robust Calibration 
#################################
# 1) Files description
#################################


This is the main code folder. It contains the following file:

1. E_Synthetic.m in the Models folder: ODEs implementation of the E_Synthetic model  
2. data_E_Synthetic.mat in the Data folder: dataset according to which the calibration is performed
3. Robustness_time_series.m: script for parameter perturbation, simulation of the ODE model and computation of the distance functions
4. BaseHist_multiplearea_piuNr.m: script for performing the intersection between distance functions and for estimating conditional parameter 
   densities. It runs the following functions: upperLowerSet_Nr.m and intersection.m
5. upperLowerSet_Nr.m: selection of the samples of the parameter vector for which all distance functions have values under or above the 
   corresponding tolerance value chosen
6. intersection.m: intersection between the selected values of distance functions of all observables 
7. BaseHist_MIRI.m: script for the computation of the Moment Independent Robustness Indicator (MIRI)



#################################
# 2) Usage
#################################

1. Run Robustness_time_series.m, setting the following parameters:
            - time_points: time points where data are sampled
            - parameters of the model. 
            - LBpi and UBpi: lower and upper bound of the perturbed parameter space 
            - Nr and NSample: number of realizations and number of parameter samples for each realization
            - model: pointer to the model to calibrate
            - Nobs: number of measured output variables

This script creates two .mat file where all the results are stored.

2. Run BaseHist_multiplearea_piuNr.m, setting the following parameters:
           - name of the target region chosen for calibration
           - protein_name: names of the proteins selected for calibration
           - lowThr or/and highThr: array of chosen thresholds for each distance function
This script calls upperLowerSet_Nr.m and intersection.m.
It computes: xbin_p and ks_p, cell arrays for storing conditional parameter densities
            CloudCondL: modes of the parameters density functions conditioned to the region of interest. 

3. (Optional) Run BaseHist_MIRI.m: if both the calibration toward low values and high values of distance functions have been
   performed, MIRI can be computed, after loading results obtained in the previous steps.
This script returns: MIRIT: array containing, for each realization, MIRI indexes for each parameter 
                     figures of the parameter density functions, conditioned to both opposite regions
           
        


#################################
# 3) OS and Matlab version
#################################

This code has been tested on Ubuntu 16.04 LTS (64bit) using Matlab R2014a.
