# ClinicianInitiatedDelivery
Code and examples for the method used in Health Outcomes Associated with Clinician-Initiated Delivery for Hypertensive Disorders at 34-38 Weeks’ Gestation

Because of confidentiality reasons, we cannot provide the original data that was analyzed in the paper. Thus, we provide synthetized data for week 38 that display the methods that were used, so that they can be implemented by other researchers. 

The main file that run the program is 3_main.R. You will need to use the setwd command to ensure that all of the files are in the same folder. 

The mitss.R file implements the MITSS method as described in Gutman, R. & Rubin, D. B. (2015), Estimation of causal effects of binary treatments in unconfounded studies. Statistics in medicine, 34(26), 3381-3398.

The mitss_func.R file includes functions are helper functions for the MITSS algorithm.

The synth_week_38.RDS file includes the synthesized data that enables to see the implementation of the code.

The synth_ps_week_38.RDS file includes the variable used to generate the propensity score.

