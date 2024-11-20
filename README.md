# Shape-based Partially Linear Single-index Cox Model for Alzheimer’s Disease Conversion
The R codes for implementing shape-based partially linear Single-index Cox model (SPLS-Cox) in simulation studies and real data analysis. The aim of this paper is to propose a shape-based
SPLS-Cox model that can handle both scalar and shape predictors. This new development is motivated by establishing the likelihood of conversion to AD in 372 patients with mild
cognitive impairment (MCI) enrolled in the Alzheimer's Disease Neuroimaging Initiative and the early shape-based markers of conversion extracted from the brain white matter region, corpus callosum (CC). These 372 MCI patients
were followed over 48 months, with 161 MCI participants progressing to AD at 48 months. Our SPLS-Cox model establishes both the estimation procedure and the point-wise confidence band. Simulation studies are conducted to
evaluate the finite sample performance of our SPLS-Cox. The real application reveals that the CC contour shape is a significant predictor for AD conversion.

## Simulation studies
* `simulation_continuous.R` :  scenarios with continuous data for the proposed and compared methods.
* `simulation_survival.R` :  scenarios with survival data for the proposed and compared methods.

## Real data application
* `hcc-simplified.csv` : synthetic hepatocellular carcinoma data with noise, for illustration purpose.
* `ITR-covgpsmatch-funcs.R` : the functions for implementation, including matching with covariates and generalized propensity scores.
* `hcc-covgpsmatch.R` : the main function.
