# Shape-based Partially Linear Single-index Cox Model for Alzheimer’s Disease Conversion
The R codes for implementing shape-based partially linear Single-index Cox model (SPLS-Cox) in simulation studies and real data analysis. The aim of this paper is to propose a
SPLS-Cox model that can handle both scalar and shape predictors. This new development is motivated by establishing the likelihood of conversion to AD in 372 patients with mild
cognitive impairment (MCI) enrolled in the Alzheimer's Disease Neuroimaging Initiative and the early shape-based markers of conversion extracted from the brain white matter region, corpus callosum (CC). These 372 MCI patients
were followed over 48 months, with 161 MCI participants progressing to AD at 48 months. Our SPLS-Cox model establishes both the estimation procedure and the point-wise confidence band. Simulation studies are conducted to
evaluate the finite sample performance of our SPLS-Cox. The real application reveals that the CC contour shape is a significant predictor for AD conversion.

## Data
* `ccdata.csv`: The CC shape data of 372 MCI patients, a 372*200 matrix.
* `CCinfo.csv`: The CC subregion volumetric data of 372 MCI patients, we segmented the CC into five regions including  anterior, central, mid-anterior, mid-posterior, and posterior. We also calculated the volume of each subregion across all patients
* `clinical.dat`: The clinical scalar covariates of these patients, including demographic and APOE information. We exclude the 331 row due to the missingness of the corresponding ccdata and CCinfo.

## Simulation studies
### code for proposed method in simulation study
* `simulation_shape.R`:  Simulation studies in the main text, which synthetic data that mimics the CC shape data. The function *generator* generates the data based on different configurations. The function ***parasimulation*** returns the estimation results, which can be implemented via parallel computing.
* `shape_RMSE_result.R`: Summarizing the estimation result.
* `shape_pointwise_result.R`: Developing the point-wise confidence band and summarizing the result. The function ***get_PWCB*** is related to the algorithm in Section 2.4.

### code for proposed method in additional simulation study in supplementary material
Additional simulation results based on general simulation settings, which can be found in the Supplementary Material. The codes are analogous to that in shape-based simulation.

### code for proposed method in additional simulation study in supplementary material
The code for implementing and summarizing functional linear cox regression model (FLCRM, Kong et al., 2018) in the scenarios above. For comparison the C-index with the proposed method.

## Real data application
* `ccdata_proposed.R`: The implementation of SPLS-Cox on real data via function ***CC_estimate*** and parallel computing.
* `ccdata_result.R`: Developing the point-wise confidence band and summarizing the result. The function ***get_ci*** is analogous to the ***get_PWCB*** previously mentioned.
* `ccdata_FLCRM.R`: implementing and summarizing the FLCRM with both scalar and shape predictors.
* `ccdata_onlyscalar.R`: implementing and summarizing standard Cox methods based on only scalar predictor.
