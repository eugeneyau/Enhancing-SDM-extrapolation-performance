# Enhancing transferability of species distribution models – guidance through algorithm, variable, and ensemble modelling decisions

### Eugene Yu Hin Yau, Alice C Hughes, Timothy C Bonebrake
[![](https://img.shields.io/badge/Citation-Ecography-blue)](https://doi.org/link)

Species distribution models (SDMs) are useful tools to estimate species’ distributions in unsampled, novel environments (e.g., future climates). However, data available often fail to represent the full range of viable environments for various reasons, including sampling biases, biogeographic barriers, and the availability of potentially suitable habitats. Therefore, many SDMs are likely based on truncated niche spaces and will have to extrapolate in unsampled environments, ultimately over- or under-estimating species’ distributions, especially under climate change. 
<br /><br />
We generated four occurrence datasets for 1,567 species, deliberately truncating the input data by subsampling four different latitudinal ranges. These datasets were used to explore the consequences of niche truncation on SDM transferability (ability to accurately extrapolate in unsampled conditions) by comparing these truncated ranges with species distributions generated from the complete dataset. We also tested the ability of various algorithms, the use of various thresholds, and variable selection to counteract truncated data on model accuracy.
<br /><br />
Niche truncation resulted in a 41.5% increase in mean SDM inaccuracy when 21.7% areas were unsampled. SDMs tend to overpredict species distributions when extrapolating in unsampled environments, distorting our estimation of climate change response and invasion risk. We suggest using a set of SDM algorithms tailored for modelling goals when extrapolating SDMs. We also recommend an inclusive approach in variable selection when ecological knowledge is insufficient to precisely identify the most relevant variables. Moreover, the choice of ensemble algorithms significantly impacts SDM transferability, which means that variation in predictions resulting from different methods should be clearly reported. While niche truncation is best dealt with using data more representative of species’ full range of viable environments, decisions related to algorithm, variable, and ensemble choices can improve model transferability.
<br />

 

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-extrapolation-performance/blob/main/readme_figs/Updated%20diversity%20diff%20EMmean.png" width=900>    

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





# Table of Contents

- Tropical Asian butterfly occurrence records from [Yau et al. (2024)](https://ecoevorxiv.org/repository/view/7470/). 
  - [`Tropical Asian Butterfly Occurrence Database`](https://doi.org/10.6084/m9.figshare.25037645)
- R scripts used in this paper:
  - [`Code/Virtual species`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)

  - [`Code/Truncation`](https://github.com/)
  
  - [`Sample`](https://github.com/)
     - [`R Markdown file`](https://github.com/)
     - [`R script`](https://github.com/)
- All files essential for running our R scripts can be downloaded from our [`Figshare repository (TrunSDMsupp_files.zip)`](https://figshare.com/s/feef1f9467edabf71a97) before running our codes.
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean), view in [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:[`Code/Variables`](https://github.com/eugeneyau/Enhancing-SDM-transferability/blob/main/Code/Variables/GEE_NDVImean.txt)

