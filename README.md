# Code for Reichel et al: ALBA proteins facilitate cytoplasmic YTHDF-mediated reading of m6A in plants

Data objects required for running the code and reproducing our results are uploaded to Zenodo: 10.5281/zenodo.11241987

Supplementary processed data together with figure source data is available at Reichel et al. 2024.

---

Description of code for individual analyses:

### Scripts for ALBA and ECT2 binding analysis using iCLIP and HyperTRIBE

### Scripts for motif analysis

### Scripts for deep learning modelling and subsequent analysis

1. ALBA_ECT2_neural_network_training_data.R

This script takes as input the atlas of ~48K m6A sites in Arabidosis, annotates them according to presence or absence of ALBA4 (iCLIP, Reichel et al. 2024) and ECT2 (iCLIP, Arribas-Hernandaz et al. 2021) and outputs files necessary for training the dual-output deep learning model. Features are 601 nt sequences (saved in .fasta format) and all sites on the same gene are kept to the same fold to prevent possible leakage of features between different folds.

2. m6a_model_5fcv.py

Script for processing the m6A input data, training the m6A neural network and making held out predictions in the vicinity of miCLIP m6A sites. Processed sites are available as supplementary in Reichel et al. 2024.

3. iclip_model_5fcv.py

Script for training the multi-output neural network for predicting RNA-protein binding of ECT2 and ALBA4 in the vicinity of m6A sites.

4. iclip_model_motif_analysis.R

Script which 

### Scripts for HyperTRIBE interaction analysis

1. HyperTRIBE_interaction_statistical_modelling.R

This script is the core statistical modelling for the robust approach based on INLA (see Reichel et al. 2024 for methodological explanation). It takes as input data which has been preprocessed by the HyperTRIBE analysis pipeline HyperTRIBER (Rennie et al. 2021). Data objects which can be used as input to the above script are available on Zenodo (link above).

2. HyperTRIBE_interaction_plot_results.R

This script produces plots presented in Figure 7 and Supplementary Figure 11 from the paper, basd on the output of script 1.

3. HyperTRIBE_interaction_SUP_3_WAY_VENN.R

Code for producing plots comparing the original HyperTRIBER output with the robust method (script 1.) and the two-sample based method (for matching ADAR between sets of fusion samples).
