## BRT calibration for Arid Central Asian brGDGT studies

# Project description
This GitHub project is associated to

# How to install/run the ACADB brGDGT calibration ?
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button <> Code beyond. 
# Use it for your paleo brGDGT data 
1. Import your paleo data:
	- The brGDGT paleo matrix should be in `.csv` using comma as separator and . as decimal
	- The name of the column should be similar to the one of `XRD.csv`, the paleo exemple
	- Import your `core.csv` file into `./Import/Paleo/` folder
	- (Optional) build a metadata matrix including the core name, Latitude, Longitude and surface climate parameters

2. Open the main.R file in Rstudio

3. Package installation
	- Verify you have the needed packages installed (i.e. `reshape2, ggplot2, patchwork, randomForest, dplyr, palaeoSIG, caret, gbm`)
	- If not use the command `install.packages("reshape2")` for instance
	
4. Modify the scrip with your own paleo brGDGT data
	- In the section `#### Import paleo data ####` replace the `GDGT.XRD` path to your own `core.csv` matrix
	- (Optional) same for the `metadata_core.csv` file

5. Run the script
