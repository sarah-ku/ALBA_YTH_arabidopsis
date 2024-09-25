# Code for Reichel et al: ALBA proteins facilitate cytoplasmic YTHDF-mediated reading of m6A in plants

Data objects required for running the code and reproducing our results will be uploaded to Zenodo: 10.5281/zenodo.11241987

---

Description of code for individual analyses:

### Scripts for ALBA and ECT2 binding analysis using iCLIP and HyperTRIBE

### Scripts for motif analysis

### Scripts for deep learning modelling and subsequent analysis

### Scripts for HyperTRIBE interaction analysis

1. HyperTRIBE_interaction_statistical_modelling.R

This script is the core statistical modelling for the robust approach based on INLA (see Reichel et al. 2024 for methodological explanation). It takes as input data which has been preprocessed by the HyperTRIBE analysis pipeline HyperTRIBER (Rennie et al. 2021). Data objects which can be used as input to the above script are available on Zenodo (link above).

2. HyperTRIBE_interaction_plot_results.R

This script produces plots presented in Figure 7 and Supplementary Figure 11 from the paper, basd on the output of script 1.

3. HyperTRIBE_interaction_SUP_3_WAY_VENN.R

Code for producing plots comparing the original HyperTRIBER output with the robust method (script 1.) and the two-sample based method (for matching ADAR between sets of fusion samples).
