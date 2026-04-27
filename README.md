# malaria_serpepi_sampling
Supporting materials for Ko YK, Li S, Britton T, Kagaya W, Kaneko A, 2025. "Optimizing Age-Structured Sampling for Estimating the Seroconversion Rate in Malaria Seroepidemiology: A Simulation Study" <br>

## Main scripts
- [functions.jl](https://github.com/KoKYura/malaria_seroepi_sampling/blob/main/functions.jl) : Julia code for the simulations in the paper 
- [app.R](https://github.com/KoKYura/malaria_seroepi_sampling/blob/main/app.R) : The Shiny app code 

## Instructions for the shiny application
To run the shiny application locally on your computer, follow the steps below:

#### 0. Install R, Julia, and RStudio:
If you have not already installed R, Julia, and RStudio, download and install them.
- R: https://cran.r-project.org/
- Julia: https://julialang.org/
- RStudio: https://posit.co/download/rstudio-desktop/

#### 1. Download the shiny app files:
You can download all the necessary files from here (https://github.com/KoKYura/malaria_seroepi_sampling)

#### 2. Install required R packages:
Open Rstudio and install the necessary packages by running the following commands in the R console.
```r
install.packages(“shiny”)
install.packages(“shinydashboard”)
install.packages(“showtext”)
install.packages(“DT”)
install.packages(“JuliaCall”)
```
#### 3. Configure Julia from R:
After installing the R packages, set up the Julia environment through JuliaCall by running:
```r
library(JuliaCall)
julia_setup()
```
If Julia is not detected automatically, specify the path to your Julia installation manually within julia_setup().
The application additionally depends on the following Julia packages: [Distributions, DataFrames, JLD2, Random, Optim]. These should be installed in the local Julia environment before launching the app by running:

#### 4. Run the application:
Open the “app.R” file in Rstudio, and launch the application either clicking the “Run App” button in Rstudio or running the following command in the R console:
shiny::runApp()

#### 5. Select a scenario from the left panel:
Choose the simulation scenario based on your interest: 1) Stable SCR or 2) Change in SCR.

#### 6. Set parameters and run simulations:
Modify the following parameters:<br>
<br>
For the Stable SCR scenario,
- Number of samples
- Seroconversion rate (SCR)
- Seroreversion rate (SRR)
- Estimating parameter (SCR only, SCR and SRR)<br>
<br>
For the Change in SCR scenario,

- Number of samples
- SCR before the change point
- SCR after the change point
- SRR
- Change point (Years)
- Estimating parameter (SCRs only, SCRs and SRR)<br>
<br>
After setting all parameters, click “Run” to start the simulation. Once the simulation is complete, you will obtain the the 95%CI for the SCR(s) along with the age-based sampling structures that yield the highest and lowest precision.

