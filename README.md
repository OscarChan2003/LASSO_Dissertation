# A Comparative Study of Variable Selection Methods for Multiply Imputed Data

This repository contains the R code and simulation scripts used for the dissertation:

**"A Comparative Study of Local Quadratic Approximation and Group Coordinate Descent for Variable Selection with Multiply Imputed Data"**

*Author: Shing Him Chan*
*Institution: Durham University*
*Date: April 24, 2025*

---

## Overview

Missing data is a common challenge in statistical analysis. Multiple Imputation (MI) is a standard technique to handle missingness, but performing variable selection consistently across imputed datasets requires specialized methods. This study compares two main approaches for penalized regression-based variable selection in the context of multiply imputed data:

1.  **Local Quadratic Approximation (LQA)**: Based on the MI-LASSO framework proposed by Chen & Wang (2013).
2.  **Group Coordinate Descent (GCD)**: Based on the GLASSO/GaLASSO framework adapted for MI data by Du et al. (2022).

The study evaluates the performance (variable selection accuracy, estimation accuracy, computational efficiency) of standard and adaptive versions of these methods under various simulation conditions, including different numbers of predictors (`p=10`, `p=20`), levels of missingness, and signal-to-noise ratios (SNR).

Furthermore, this work introduces:
*   An **adaptive version of the LQA-based MI-LASSO (MI-aLASSO)** using BIC tuning.
*   **Cross-validation (CV)** tuning procedures for the LQA-based methods (MI-LASSO and MI-aLASSO) to compare against the original BIC tuning and the CV tuning used in GCD methods.

---

## Repository Structure

This repository contains the primary R scripts for running the simulations described in the dissertation.

**LQA-based Methods (MI-LASSO/MI-aLASSO with BIC or CV tuning):**

*   `cv_MI_aLASSO_p10.R`: Cross-validated MI-aLASSO (LQA) for p=10.
*   `cv_MI_aLASSO_p20.R`: Cross-validated MI-aLASSO (LQA) for p=20.
*   `cv_MI_LASSO_p10.R`: Cross-validated MI-LASSO (LQA) for p=10.
*   `cv_MI_LASSO_p20.R`: Cross-validated MI-LASSO (LQA) for p=20.
*   `MI_aLASSO_p10.R`: MI-aLASSO (LQA, BIC-tuned) for p=10.
*   `MI_aLASSO_p20.R`: MI-aLASSO (LQA, BIC-tuned) for p=20.
*   `MI_LASSO_p10.R`: MI-LASSO (LQA, BIC-tuned) for p=10.
*   `MI_LASSO_p20.R`: MI-LASSO (LQA, BIC-tuned) for p=20.

**GCD-based Methods (GLASSO/GaLASSO with CV tuning):**

*   `GaLASSO_p10.R`: GaLASSO (GCD) for p=10.
*   `GaLASSO_p20.R`: GaLASSO (GCD) for p=20.
*   `GLASSO_p10.R`: GLASSO (GCD) for p=10.
*   `GLASSO_p20.R`: GLASSO (GCD) for p=20.

**Note:** Files ending in `_p10.R` are for simulations with p=10 predictors, and `_p20.R` are for p=20 predictors. The separation into multiple files facilitates running different methods and parameter settings in parallel.

---

## Dependencies

The required R packages differ depending on the method group:

### For LQA-based Methods (MI-LASSO, MI-aLASSO, and their CV versions)

*   **`mice`**: Used for performing Multiple Imputation using Chained Equations.
    ```R
    install.packages("mice")
    ```
    *(No other external packages are required for these specific scripts, beyond base R functionality).*

### For GCD-based Methods (GLASSO, GaLASSO)

*   **`mice`**: Also used here for the initial multiple imputation step.
    ```R
    # If not already installed
    # install.packages("mice")
    ```
*   **`devtools`**: Needed to install the `miselect` package from GitHub.
    ```R
    install.packages("devtools")
    ```
*   **`miselect`**: Implements the grouped coordinate descent methods for multiply imputed data (from Du et al., 2022).
    ```R
    devtools::install_github("umich-cphds/miselect")
    ```

---

## Methods Compared & Developed

### Existing Methods Evaluated

*   **MI-LASSO (BIC-tuned)**: Chen & Wang (2013), LQA approach.
*   **GLASSO (CV-tuned)**: Du et al. (2022) adaptation of Group LASSO using GCD for MI data.
*   **GaLASSO (CV-tuned)**: Du et al. (2022) adaptation of Group Adaptive LASSO using GCD for MI data.

### Methods Developed in this Dissertation

*   **MI-aLASSO (BIC & CV tuned)**: An adaptive version of MI-LASSO using the LQA framework, incorporating adaptive weights inspired by Zou (2006) and Wang & Leng (2008), with weight calculation adapted from Du et al. (2022) for the MI context. Implemented here with both BIC and CV tuning.
*   **MI-LASSO (CV-tuned)**: Implementation of k-fold cross-validation for selecting the tuning parameter (`lambda`) for the original MI-LASSO method, providing a direct comparison to BIC and the CV used in GCD methods.

---

## Usage

1.  Ensure all required R packages for the desired method(s) are installed (see Dependencies).
2.  Clone the repository.
3.  Open the desired R script (e.g., `cv_MI_aLASSO_p20.R`).
4.  Modify simulation parameters within the script if needed (e.g., number of simulation `replicates`, `st.seed`, `SNR_list`, `missing_props`, `n`, `p`, number of imputations `D`).
5.  Run the script in R or RStudio.

**Parallel Execution on HPC:**
*   The simulation study is divided into multiple distinct R script files (separated by method, tuning approach, and predictor size `p`). This structure is designed to allow for parallel execution of different simulation settings, which was utilized on the Durham University HPC (Hamilton) to complete the study efficiently. Each script can be submitted as an independent job.

*Note: Simulations, especially with multiple imputations and cross-validation, can be computationally intensive.*

---

## Output

Each simulation script performs the following steps at the end of its run:

1.  **Aggregation across Imputations:** Calculates average performance metrics (coefficients, MSE, computation time) and determines variable selection consistency (SEN, SPE, MCC) across the `D` multiple imputations for each specific simulation seed and parameter setting.
2.  **Aggregation across Seeds:** Calculates the final average and standard deviation of the performance metrics (SEN, SPE, MCC, MSE, computation time) across all simulation `replicates` (seeds) for each unique combination of simulation parameters (missing proportion, SNR).
3.  **Saving Results:** Creates a local directory named according to the date, predictor size (`p`), and method (e.g., `YYYY-MM-DD_p10_MI-LASSO`). Within this directory, it saves three `.csv` files containing:
    *   `res_... .csv`: Raw results for every imputation within every seed.
    *   `aggregated_res_... .csv`: Results aggregated across imputations for each seed.
    *   `combined_res_... .csv`: Final results aggregated across all seeds for each parameter combination (this typically contains the main summary statistics reported in the study).

---

## Citation

If you use the code or findings from this study, please cite the original dissertation:

*   Chan, S. H. (2025). *A Comparative Study of Local Quadratic Approximation and Group Coordinate Descent for Variable Selection with Multiply Imputed Data*. Durham University.

Please also cite the relevant papers for the methods used:

*   Chen, Q., & Wang, S. (2013). Variable selection for multiply-imputed data with application to dioxin exposure study. *Statistics in Medicine*, *32*(21), 3646–3659.
*   Du, J., Boss, J., Han, P., Beesley, L. J., Kleinsasser, M., Goutman, S. A., Batterman, S., Feldman, E. L., & Mukherjee, B. (2022). Variable selection with multiply-imputed datasets: choosing between stacked and grouped methods. *Journal of Computational and Graphical Statistics*, *31*(4), 1063–1075.
*   Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software*, *45*(3), 1–67.
*   Kleinsasser, M., Rix, A., & Du, J. (2024). *miselect: Variable selection for multiply imputed data*. R package version 0.9.2. [https://github.com/umich-cphds/miselect](https://github.com/umich-cphds/miselect)

---

## Acknowledgement

This work made use of the Hamilton HPC Service of Durham University to generate the simulation studies.
