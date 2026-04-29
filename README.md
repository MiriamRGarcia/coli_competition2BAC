# coli_competition2BAC


**CompExp_WvsT_YFP_mCherry_BAC_M9**: Repository with experimental data and codes presented in: “N. Martínez-López, A. Pedreira, M. R. García, F. Schreiber and N. Nordholt - A growth survival trade-off quantitatively predicts microbial selection under periodic disinfection.”

## The repository consists of the following folders:

**Data**: Experimental data used in the study, obtained from a preliminary time-kill assay and a competition experiment performed with strains: 

- Wild-type E. coli tagged with YFP (W-YFP), 
- Wild-type E. coli tagged with mCherry (W-mCh),
- Tolerant E. coli tagged with YFP (T-YFP), 
- Tolerant E. coli tagged with mCherry (T-mCh).

The preliminary time-kill assay was performed for the strains separately, following the steps: 

- (i) Overnight cultures of the strains starting from an inoculum of approximately 1e5 CFU/mL were grown for 24 hours in M9 medium (n=4 biological replicates per strain) till the stationary phase.
- (ii) Benzalkonium chloride (BAC) was applied to the cultures at concentrations: BAC = 30,40,50,60,75,150 µg/mL.
- (iii) Cell counts were recorded for the different BAC concentrations at sampling times: ts = 0,5,10,20 min.

The competition experiment was performed for the combinations: W-YFP vs T-mCh, and W-mCh vs T-YFP. Strains W and T were fluorescently labeled with YFP or mCherry to separate cells using flow cytometry. The competition experiment was performed following the steps:

- (i) Overnight cultures of W-YFP were mixed with T-mCh, and W-mCh with T-YFP at different initial mixing ratios (i.e., varying the initial proportion of T to W cells) to a total inoculum of approximately 1e6 CFU/mL (initial mixing ratios: T0:W0 = 1:1, 1:102, 1:104). Growth periods lasted 24 hours in M9 medium, and disinfection periods lasted 10 minutes.
- (ii) The mixed cultures were subjected to successive cycles alternating growth with disinfection periods, testing different BAC concentrations (BAC = 0,30,40,50 µg/mL).
- (iii) Subpopulation frequencies for YFP and mCherry cells were determined before each BAC disinfection period using flow cytometry. Additionally, total cell counts were also determined.

**Codes**: MATLAB codes for data processing and calculations.


### Folder **Data** contains:

- **PreliminaryTKC_WvsT_YFP_mCherry_BAC_M9.xlsx**: Excel tables with cell counts recorded in the preliminary time-kill assay.
- **CellCounts_CompExp_WvsT_YFP_mCherry_BAC_M9.xlsx**: Excel tables with cell counts recorded in the competition experiment (before each disinfection period).
- **FlowCytometry_CompExp_WvsT_YFP_mCherry_BAC_M9**: Folder containing raw flow cytometry data recorded during the competition experiment (.fcs format). The folder is organized into subfolders that differentiate the cycles of the competition experiment during which data were recorded.

### Folder **Codes** contains:
- Main functions to run calculations:
    - Run_Model.m: MATLAB script to run model predictions of the competition experiment. The script calculates and plots: (i) The microbial traits (growth rates and survival fractions) and confidence intervals for the competing strains, (ii) Model-based subpopulation frequencies of the competing strains during the competition experiment vs flow cytometry data, (iii) Model-based predictions of the subpopulation dynamics for the competing strains, (iv) Model-based predictions of the Selection-Extinction (SE) planes, (v) Analytical approximation of the model-based predictions for the extinction cycles for the competing strains.
    - Run_ANOVA_TKC.m: MATLAB script to run two-way repeated measures ANOVA to validate the hypothesis that cell survival for strains W and T under BAC exposure is independent of the fluorescent label (YFP or mCherry).
    - Run_Fit_FCdata.m: MATLAB script to fit the microbial traits (growth rates and survival fractions) of the competing strains against the flow cytometry data recorded during the competition experiment and compare with the experimental values used to parametrize the model (determined from preliminary experiments).




- Subfolder Functions: Auxiliary MATLAB functions to assist the main scripts. Contains:
    - **Sim_Model.m: Function implementing the model of the competition experiment.
    - **Plot_Model_vs_Data.m**: Function to plot model-based predictions for the competition experiment and compare with flow cytometry data.
    - **Plot_SE_Plane**: Function to run and plot model-based predictions of Selection-Extinction (SE) planes for the competition experiment.

- Subfolder ProcessData: Auxiliary MATLAB codes to read and process experimental data. Contains:
    - **Read_TKC.m**: MATLAB script to read and process cell count data for the preliminary time-kill assay. The results are stored in the MATLAB data file Data_TKC.mat
    -  **Read_Traits.m**: MATLAB script to read and process microbial traits (growth rates and survival fractions) for the competing strains and obtain confidence intervals. The results are stored in the MATLAB data file Data_Traits.mat
    -  **Read_CompExp.m**: MATLAB script to read and process cell counts and flow cytometry data recorded during the competition experiment. The script uses the auxiliary function Process_FCdata.m. The results are stored in the MATLAB data files: Data_CompExp_CellCounts.mat and Data_CompExp_FC.mat. 
    - **Process_FCdata.m**: Auxiliary MATLAB function to process raw cytometry data recorded during the competition experiment. The function cleans background noise and debris based on FSC-A vs SSC-A, and doublets based on SSC-H vs SSC-A.
    - **Train_Clusters_FCdata.m**: MATLAB script to train a Gaussian Mixture Model (GMM) algorithm for clustering flow cytometry data recorded during the competition experiment. The results are stored in the MATLAB data file Data_CompExp_TrainedGMM.
    - **Process_FCclusters.m**: MATLAB script to calculate subpopulation frequencies of YFP and mCherry during the competition experiment from clustered flow cytometry data. The script uses the auxiliary function ClusteringGMM.m The results are stored in the MATLAB data file Data_CompExp_FC_Clustered
    - **ClusteringGMM.m**: Auxiliary MATLAB function to perform clustering of the flow cytometry data recorded during the competition experiment using the trained GMM model. 


Note! The MATLAB function fca_readfcs is necessary in Read_CompExp.m to read flow cytometry data format (.fcs), and ESS toolbox is necessary in Run_Fit_FCdata.m for fitting the model to flow cytometry data.
