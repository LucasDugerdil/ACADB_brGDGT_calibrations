# BRT calibrations for Arid Central Asian for paleo brGDGT 

## Overview
This GitHub project is associated to the publication of *Boosted Regression Trees machine-learning method drastically improves the brGDGT-based climate reconstruction in drylands.* in *Paleoceanography and Paleoclimatology* in 2025 (Dugerdil et al., 2025c) and accessible at **DOI:** [10.1029/2025PA005214](https://doi.org/10.1029/2025PA005214)

**Author**: **Lucas Dugerdil**<sup>1,2</sup>

**Affiliations**: 
1. Key Laboratory of Western China's Environmental Systems (Ministry of Education), College of Earth and Environmental Sciences, Lanzhou University, Lanzhou, China
2. Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE, F-69364, Lyon, France  
3. Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM, Montpellier, France  
**ORCID:** [0000-0003-0266-564X](https://orcid.org/0000-0003-0266-564X)

## Description
This R script permits to easily apply the BRT calibration trained on the ACADB and the two subsets *K-warm/arid* and *K-cold/wet* for your own paleo brGDGT datas. The GDGT data from the XRD loess-palaeosol sequence are provided in the repository as a data exemple.

## How to install/run the ACADB brGDGT calibrations?
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button `<> Code` beyond. 

## Use it for your paleo brGDGT data 
1. Import your paleo data:
	- The brGDGT paleo matrix should be in `.csv` using comma as separator and . as decimal
	- The name of the column should be similar to the one of `XRD.csv`, the paleo exemple
	- Import your `core.csv` file into `./Import/Paleo/` folder
	- (Optional) build a metadata matrix including the core name, Latitude, Longitude and surface climate parameters

2. Open the `ACADB_brGDGT_calibrations.Rproj` file in Rstudio

3. Package installation
	- Verify you have the needed packages installed (i.e. `reshape2, ggplot2, patchwork, randomForest, dplyr, palaeoSIG, caret, gbm`)
	- If not use the command `install.packages("reshape2")` for instance
	
4. Modify the scrip with your own paleo brGDGT data
	- In the section `#### Import paleo data ####` replace the `GDGT.XRD` path to your own `core.csv` matrix
	- (Optional) same for the `metadata_core.csv` file

5. Run the script
	- The function `GDGT.histo.plot.surf.core()` allows you to verify that no mistake appear in your raw brGDGT fractional abundances distributions
	- The function `FT.core()` make the BRT prediction for your brGDGT paleo data. It is runned three times for each dataset (i.e. ACADB, *K-warm/arid* and *K-cold/wet*)
	- The function `Combine.ML.cluster()` build the weighted combined model and plot the BRT results
	- (Optional) to apply the `randomTF()` test on your different models, turn to `TRUE` the test change line 53 into `test.randomTF = T`, then the script will launch the function `Plot.randomTF()`
	- (Optional) many option and settings for each function (e.g. change the `BRT` model for the `RF`, export figure in `plotly`, etc.) can be discovered when look at the `./Import/Script/BRT_script.R` file
