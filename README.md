# Shape-based Partially Linear Single-index Cox Model for Alzheimer’s Disease Conversion
The R codes for implementing shape-based partially linear Single-index Cox model (SPLS-Cox) in simulation studies and real data analysis. The aim of this paper is to propose a
SPLS-Cox model that can handle both scalar and shape predictors. This new development is motivated by establishing the likelihood of conversion to AD in 372 patients with mild
cognitive impairment (MCI) enrolled in the Alzheimer's Disease Neuroimaging Initiative and the early shape-based markers of conversion extracted from the brain white matter region, corpus callosum (CC). These 372 MCI patients
were followed over 48 months, with 161 MCI participants progressing to AD at 48 months. Our SPLS-Cox model establishes both the estimation procedure and the point-wise confidence band. Simulation studies are conducted to
evaluate the finite sample performance of our SPLS-Cox. The real application reveals that the CC contour shape is a significant predictor for AD conversion.

## Data
* `ccdata.csv`: The CC shape data of 372 MCI patients, a 372*200 matrix.
* `CCinfo.csv`: The CC subregion volumetric data of 372 MCI patients, we segmented the CC into five regions including  anterior, central, mid-anterior, mid-posterior, and posterior. We also calculated the volume of each subregion across all patients
* `clinical.dat`: The clinical scalar covariates of these patients, including demographic and APOE information.

## Simulation studies
### code for proposed method in simulation study
* `simulation_shape.R`:  Simulation studies in the main text, which synthetic data that mimics the CC shape data. The function *generator* generates the data based on different configurations. The function ***parasimulation*** returns the estimation results, which can be implemented via parallel computing.
* `shape_RMSE_result.R`:  Summarizing the estimation result.
* `shape_pointwise_result.R`: Developing the point-wise confidence band and summarizing the result. The function *get_PWCB* is related to the algorithm in Section 2.4.

## Real data application
* `hcc-simplified.csv` : synthetic hepatocellular carcinoma data with noise, for illustration purpose.
* `ITR-covgpsmatch-funcs.R` : the functions for implementation, including matching with covariates and generalized propensity scores.
* `hcc-covgpsmatch.R` : the main function.
