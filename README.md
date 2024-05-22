# VFA-optimization

INSTRUCTIONS FOR USE
This is the companion codes for the manuscript: Tay et al. Biosens. Bioelectron. 2023, 222.

By Dousabel May Yi Tay (dousabel@mit.edu)
_____________________________________________________________________

>> Base code
- This folder contains models to predict signals for the three different assay formats.
  > Premix_VFA_parfor_nHRP.m: Uses the Premix model reported and primarily used in the reported paper. 
	Takes in sample concentration, cells of operating conditions (see code) to output expected cyan intensity
  > Sequential_VFA_parfor_nHRP_only.m: Applies a sequential load of antigen and label to predict cyan intensity.
	Refer to Fig. 4 in associated paper.
  > Top_Readout_Soln_Premix_parfor_nHRP_only: Applies a solution mixing of capture, reporter and antigen to predict cyan intensity.
	Refer to Fig. 4 in associated paper. This uses Top_Readout_Rxn.m, which contains the relevant set of ODEs to predict the complexes formed by solution mixing.

- Note that these base codes need to be in the same folder as the simulation codes for prediction

>> Image Processing
- This folder contains the codes for MATLAB-ImageJ automated data processing
  >Process_Cyan_WhiteBkgrd.m: Uses imfindcircles to determine location of the VFA wells and sends it to ImageJ to extract Cyan intensities
	"Sample_Image.JPG" is provided as a sample image to demonstrate code application
  >ImageJ-Code: This folder contains the written macro to install into ImageJ to facilitate data processing

>> Calibration
- This folder contains a sample code to determine the fitting parameters for
	1. The variation in Cyan intensity with amount of SA-HRP on the paper
	2. Effective Capture molecule fraction & Non-specific binding parameters
  > Fit_Cyan_White_SAHRP.m: Uses "Processed_Cyan_White_data.mat" to determine the coefficients for a saturating exponential fit to Cyan intensities
  > rcSso_NP_Saliva_Calibration_EffCBD_NSB.m: Uses fmincon to determine the effective capture molecule (CBD) fraction on the paper, and subsequently freezes it
	to determine the non-specific binding parameters of HRP to cellulose.

- Relevant data should be collected, and these calibration codes should be run to determine the parameters to be used in subsequent simulations

>> Parameter Investigation
- This folder contains codes to examine the effect of Antigen concentration, HRP and Reporter concentration and sample volume on resulting signals
	These codes use the Base code Premix_VFA_parfor_nHRP.m
  > rcSso_NP_Saliva_vary_XXX_Expt_Theory.m, where XXX is (Antigen)Conc, HRPConc, RepConc, or SampleVol:
	These codes vary the parameter specified in XXX to determine the effect of signal and compares it against experimental data ("Sorted_NP_Saliva_Data.mat",See folder "Sample Data")
	Note that relevant data for flow rate and variation of signal with cyan intensity are included in the .mat files 
	"White_Saliva_Flow_rates.mat", "White_Saliva_HRP_Calibration_beta_Cleaned.mat","opt_NP_Saliva_params_all_1layer_Cellulose_NSB.mat" (See folder "Sample Data")

>> Sample Data
- This folder contains sample data to compare experimental against theoretical results for "Parameter Investigation". Refer to Figures 5G-I in manuscript 

>> Performance Metrics
- This folder contains all codes used to investigate the variation of relevant performance metrics - LoD, Sensitivity, Signal-to-Noise Ratio and Difference across multiple parameters
  > Calibration Curves
	- This folder contains the codes to obtain the Calibration curves (Figure 6A, D, G, K in manuscript)
	> rcSso_NP_Saliva_LoD_XXX.m, where XXX is HRP_Conc, IncubTime, k_on, k_off, RepConc, SampleVol: Plots calibration curves while varying the parameter XXX
	> run_LoD_range.m: runs all codes "rcSso_NP_Saliva_LoD_XXX.m" and records time taken for each
  > Determine_LoD_Sensitivity
	- This folder contains the codes to obtain the Calibration curves (Figure 6B, E, H, L in manuscript)
	> rcSso_NP_Saliva_LoD_variation_XXX.m, where XXX is HRP_Conc, IncubTime, k_on, k_off, RepConc, SampleVol: Plots calibration curves while varying the parameter XXX
	> rcSso_NP_Saliva_LoD_variation_XXX.m, where XXX is HRP_Conc, IncubTime, k_on, k_off, RepConc, SampleVol: Plots calibration curves while varying the parameter XXX
	> rcSso_NP_Saliva_LoD_variation_fun.m: Helper function to solve to obtain LoD
	> rcSso_NP_Saliva_LoD_variation_Solver.m: Solver function to obtain LoD
	> run_LoD_Sensitivity.m: runs all codes "rcSso_NP_Saliva_LoD_variation_XXX.m" & "rcSso_NP_Saliva_Sensitivty_variation_XXX.m" and records time taken for each
  > Determine_S_Nmetrics
	- This folder contains the codes to obtain the Calibration curves (Figure 6C, F, J, M in manuscript)
	> rcSso_NP_Saliva_SN_variation_XXX.m, where XXX is HRP_Conc, IncubTime, k_on, k_off, RepConc, SampleVol: Plots calibration curves while varying the parameter XXX
	> run_S_Nratio.m: runs all codes "rcSso_NP_Saliva_SN_variation_XXX.m" and records time taken for each

>> VFA Format Comparison
- This folder contains codes to compare the three different VFA formats with simulation (See Figure 4 in manuscript)
  > Model_Expt_TBH4_Formats.m: Main code to plot all predictions of cyan intensity for the three formats: Premix, Sequential and TopReadout (Mixing Capture and reporter)
	This file uses the data in the .mat files as uploaded for flow rates and relevant calibration data located in the same folder.
