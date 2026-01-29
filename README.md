# Enhancing extrapolation capability of species distribution models – guidance through algorithm, variable, and ensemble modelling decisions

### Eugene Yu Hin Yau, Alice C Hughes, Timothy C Bonebrake*
*Correspondence: tbone@hku.hk

[![](https://img.shields.io/badge/Citation-Journal_name-blue)](https://doi.org/link)

Species distribution models (SDMs) are useful tools to predict species’ potential distributions in unsampled, novel environments (e.g., future climates). However, data available often fail to represent the full range of viable environments for various reasons, including sampling biases, biogeographic barriers, and the availability of potentially suitable habitats. Therefore, many SDMs are likely based on truncated fundamental niche and will have to extrapolate in unsampled environments, ultimately over- or under-estimating species’ distributions, especially under climate change.
<br /><br />
To identify optimal modelling decisions to maximize SDM extrapolation performance, we generated 984 virtual species, then subsampled virtual occurrences from the northernmost 80% virtual presences to deliberately generate a dataset representing truncated niche. This dataset was used to train SDMs which extrapolate in unsampled environments. This workflow allows us to test how different SDM algorithms, SDM ensemble methods, and variable selection strategies improve SDM extrapolation performance.
<br /><br />
Niche truncation, caused by the deletion of southernmost 20% occurrence data, increased inaccurate predictions by 57.26% in average. We found that SDMs tend to overpredict species distributions when extrapolating in unsampled environments, distorting our estimation of climate change response and invasion risk. Depending on modelling goals, we suggest using different sets of SDM algorithms to minimize either overprediction or underprediction when extrapolating SDMs. We also recommend including all potentially relevant variables for maximal extrapolation performance. Moreover, models ensembled differently can significantly differ in their SDM extrapolation performance, therefore such variation should be reported. While niche truncation is best dealt by using data more representative of species’ full range of viable environments, decisions related to algorithm, variable, and ensemble methods can significantly improve model performance. 

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/readme_figs/Updated%20diversity%20diff%20EMmean.png" width=1000>    

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />






## Table of Contents

- R scripts used in this study:
  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Step1`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)
  
- All files essential for running our R scripts can be downloaded from our [`Figshare repository (TrunSDMsupp_files.zip)`](https://figshare.com/s/feef1f9467edabf71a97) before running our codes
- Tropical Asian butterfly occurrence records from [Yau et al. (2025)](https://doi.org/10.1038/s41597-025-05333-w) can be downloaded [`here`](https://doi.org/10.6084/m9.figshare.25037645)
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean) can be accessed from [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:[`Code/Variables`](https://github.com/eugeneyau/Enhancing-SDM-transferability/blob/main/Code/Variables/GEE_NDVImean.txt)

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/readme_figs/Methods%20flowchart%20github.png" width=800>    


