# Optimizing species distribution models for extrapolation – guidance through algorithm, variable, and ensemble modelling decisions

### [Eugene Yu Hin Yau](https://scholar.google.com.hk/citations?user=dQG4FA4AAAAJ), [Timothy C Bonebrake](https://scholar.google.com.hk/citations?user=B-_yE1YAAAAJ), [Alice C Hughes](https://scholar.google.com.hk/citations?user=jlkGFIIAAAAJ)*
*Correspondence: alice.c.hughes@unimelb.edu.au

[![](https://img.shields.io/badge/Citation-Ecological_Informatics-blue)](https://doi.org/link) 
[![](https://img.shields.io/badge/Data-Figshare-blue)](https://figshare.com/s/feef1f9467edabf71a97) 

Species distribution models (SDMs) are useful tools to predict species’ potential distributions in under-sampled, novel environments. However, available data often fails to represent the full range of viable environments for various reasons, including sampling biases, biogeographic barriers, and the availability of potentially suitable habitats. Therefore, many SDMs are likely based on truncated niches and will extrapolate in unsampled environments, ultimately over- or under-estimating species’ distributions, especially under climate change. To identify optimal modelling decisions to maximize SDM extrapolation performance, we generated 1007 virtual species, then subsampled virtual occurrences from the northernmost 80% virtual presences to deliberately generate a dataset representing truncated niche. This dataset was used to train SDMs which extrapolate in hot, non-analog environments as in climate change assessments. This workflow allows us to test how different SDM algorithms, SDM ensemble methods, and variable selection strategies improve SDM extrapolation performance. Niche truncation, caused by the deletion of southernmost 20% occurrence data, increased inaccurate predictions by more than 50%. SDMs may overpredict species distributions when extrapolating in unsampled environments, distorting our estimation of climate change response and invasion risk. When data representing full range of species’ viable environments is unavailable, decisions on algorithm, variable, and ensemble methods can improve SDM (extrapolation) performance across time and space to estimate climate change consequences and manage invasion risks. We suggest including all potentially relevant variables, reporting variations in ensemble model predictions caused by different ensemble methods, and using different sets of SDM algorithms to minimize either overprediction or underprediction depending on modelling goals.

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/readme_figs/Updated%20diversity%20diff%20EMmean.png" width=1000>    
<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

## Table of Contents

- R scripts used in this study:
  - Step 0 -- Prepare input data for running SDMs [`Code/Step0`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step0)
  
  - Step 1 -- Generate virtual species [`Code/Step1`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step1)
     - [Calculate VIF](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step1/SEAbfyTrun_analysis_REV_noise_VIF.R) for all variables
     - Generate [variable-based](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step1/SEAbfyTrun_BASE_VIF5_hpc18_1.R) or [PCA-based](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step1/SEAbfyTrun_BASE_VIF5_PC_hpc18_1.R) virtual species
     - [Generate sample size](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/Code/Step1/Supp_Vsp_gen_%20occ_sample_size.R) for variable-based/PCA-based virtual species
     - [Generate PCs](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/Code/Step1/Supp_Vsp_getPCs.R) for building PCA-based virtual species
   
  - Step 2,3,4 -- Sample virtual occurrence datasets from different latitudinal ranges, and construct SDMs from these data, project constructed SDMs to the full latitudinal range [`Code/Step2_3_4`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step2_3_4)
    - Latitudinal data truncation by percentage of data [`Code/Step2_3_4/Trunc_by_percentage`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage)<br />
      SDM algorithm (tuned) -- [`/Test_SDM_algorithms_tuned`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage/Test_SDM_algorithms_tuned)
      <br />
      SDM algorithm (untuned) -- [`/Test_SDM_algorithms_untuned`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage/Test_SDM_algorithms_untuned)
      <br />
      Choice of environmental variables -- [`/Test_env_variables`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage/Test_env_variables)
      <br />
      Ensemble methods -- [`/Test_ensemble_methods`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage/Test_ensemble_methods)
      <br />
    - Latitudinal data truncation by fixed latitude range [`Code/Step2_3_4/Trunc_by_percentage`](https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/tree/main/Code/Step2_3_4/Trunc_by_percentage)
  
   - Step 5 -- Quantify SDM performance [`Code/Step5`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step5)
     - [SDM algorithm](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step5/SEAbfyTrun_analysis_VIF5_algo.R)
     - [Choice of environmental variables](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step5/SEAbfyTrun_analysis_VIF5_noise.R)
     - [Ensemble methods (all algorithms)](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step5/SEAbfyTrun_analysis_VIF5_EM.R)
     - [Ensemble methods (selected algorithms)](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step5/SEAbfyTrun_analysis_VIF5_get_selectedalgo_EM.R)
     - [Binarisation method](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step5/SEAbfyTrun_analysis_VIF5_BINmethod.R)

  - Step 6 -- Visualise results [`Code/Step6`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step6)
     - Figure 2. Alpha diversity trends under different levels of niche truncation -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_trunlevel.R)
     - Figure 3. Extrapolation performance of different SDM algorithms -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_algo_FPFN_barplot.R)
     - Figure 4. Extrapolation performance of selected (4 best-performing) SDM algorithms -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_EM_FPFN_barplot.R)
     - Supplementary<br />
       Density plot for performance of different SDM algorithms -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_algo_ggplot.R)<br />
       Density plot for performance of different sets of environmental variables -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_noise_ggplot.R)<br />
       Density plot for performance of different ensemble methods -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_EM_ggplot.R)<br />
       Bar chart for performance of different binarisation methods -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_BINmethod_FPFN_barplot.R)<br />
       Relationship between training data sample size and SDM extrapolation performance -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step6/SEAbfyTrun_analysis_VIF5_occ_n.R)<br />
       Pairwise p-value matrices for bar charts -- [`Code/Step6/Pairwise_p-value_matrices`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step6/Pairwise_p-value_matrices)<br />
   
   - Step 7 -- Additional analyses [`Code/Step7`](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/tree/main/Code/Step7)
     - Calculate variable importance for Table 2 -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step7/SEAbfyTrun_analysis_VIF5_noise_var_imp.R)
     - Calculate VIF for environmental variables (including noise variables) -- [R script](https://github.com/eugeneyau/Optimizing-species-distribution-models-for-extrapolation/blob/main/Code/Step7/SEAbfyTrun_analysis_REV_noise_VIF.R)
  
- All files essential for running our R scripts can be downloaded from our [figshare repository](https://figshare.com/s/feef1f9467edabf71a97)

- Tropical Asian butterfly occurrence records (v1.1) from [Yau et al. (2025)](https://doi.org/10.1038/s41597-025-05333-w) can be downloaded from their [figshare repository](https://doi.org/10.6084/m9.figshare.25037645). Visit their [GitHub page](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution) to learn more.

- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean) can be accessed from [Google Earth Engine](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from [`Code/Variables`](https://github.com/eugeneyau/Enhancing-SDM-transferability/blob/main/Code/Variables/GEE_NDVImean.txt)

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/readme_figs/Methods%20flowchart%20github%20v2.png" width=800>    


