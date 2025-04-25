rm(list = ls())
library(mice)

setwd("~/LASSO_Code/results")

# Install necessary packages if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("miselect", quietly = TRUE)) {
  devtools::install_github("umich-cphds/miselect", build_opts = c())
}
library(miselect)

#### Parameters
n <- 100
replicate <- 100
st.seed <- 1
p <- 10
D <- 5

missing_props <- c(0.2,0.4,0.6)   # Proportion of missing data
SNR_list <- c(2,4,6)              # Desired Signal-to-Noise Ratios


# True coefficients (signal in first 5 predictors)
beta.true <- c(3, 2, 1.5, 1, 0.5, rep(0, p - 5))
beta.check <- beta.true != 0


#### Data generation
generate_and_process_data <- function(seed, n, p, D, desired_SNR, missing_proportion, beta) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  cov.X <- (t(X) %*% X) / n
  y_signal <- scale(X %*% beta, center = TRUE, scale = FALSE)
  var_signal <- var(as.vector(y_signal))
  var_noise <- var_signal / desired_SNR
  sigma <- sqrt(var_noise)
  y <- as.vector(y_signal + rnorm(n, sd = sigma))
  
  x_missing_indices <- sample(1:(n * p), size = floor(missing_proportion * n * p))
  X[x_missing_indices] <- NA
  #y_missing_indices <- sample(1:(n), size = floor(missing_proportion * n))
  #X[y_missing_indices] <- NA
  data <- data.frame(X, Y = y)
  
  imputed_data <- mice(data, m = D, method = 'pmm', maxit = 5, seed = seed)
  newdata <- vector("list", D)
  x_list <- vector("list", D)
  y_list <- vector("list", D)
  
  for (d in 1:D) {
    newdata[[d]] <- complete(imputed_data, d)
    x_list[[d]] <- as.matrix(newdata[[d]][, 1:(ncol(newdata[[d]]) - 1)])
    y_list[[d]] <- newdata[[d]][, ncol(newdata[[d]])]
  }
  
  return(list(cov.X = cov.X, imputed_data = imputed_data,
              newdata = newdata, x = x_list, y = y_list))
}

#### GLASSO formatting function
process_galasso_output <- function(coef_list, p, D, tol = 1e-8) {
  # Input checks
  if (!is.list(coef_list)) {
    warning("coef_list is not a list in process_galasso_output. Returning NULL.")
    return(NULL)
  }
  # Sometimes coef() might return a single matrix if D=1 or internally simplified
  if (is.matrix(coef_list) && ncol(coef_list) == D && nrow(coef_list) == (p + 1)) {
    coef_matrix <- coef_list
  } else if (is.list(coef_list) && length(coef_list) == D) {
    # Check if elements are correct
    all_elements_valid = all(sapply(coef_list, function(el) is.numeric(el) && length(el) == (p + 1)))
    if (!all_elements_valid) {
      warning("coef_list elements have unexpected format in process_galasso_output. Returning NULL.")
      return(NULL)
    }
    # Combine list elements into a matrix
    coef_matrix <- do.call(cbind, coef_list)
    if (nrow(coef_matrix) != (p+1) || ncol(coef_matrix) != D) {
      warning("Failed to bind coef_list into correct matrix dimensions. Returning NULL.")
      return(NULL)
    }
  } else {
    warning("Unexpected format for coef_list in process_galasso_output. Returning NULL.")
    return(NULL)
  }
  
  # Check for non-finite values in the coefficient matrix
  if(any(!is.finite(coef_matrix))){
    warning("Non-finite coefficients found in process_galasso_output. Returning NULL.")
    return(NULL)
  }
  
  # Calculate varsel based on non-intercept coefficients
  # Use a tolerance check: is the sum of squares across imputations > tol?
  final_b_equiv <- coef_matrix[-1, , drop = FALSE] # p x D matrix
  varsel <- rowSums(final_b_equiv^2, na.rm = TRUE) > tol # Calculate rowSums first
  
  # Ensure varsel is a logical vector of length p
  if(length(varsel) != p || !is.logical(varsel)) {
    warning("varsel calculation failed in process_galasso_output. Returning NULL.")
    return(NULL)
  }
  
  # Return structure similar to MI.lasso/MI.alasso
  return(list(coefficients = coef_matrix, varsel = varsel, final_b = final_b_equiv))
}

#### MSE function
MSE <- function(Hat.beta, beta, cov) {
  diff <- Hat.beta - beta
  MSE_value <- as.numeric(t(diff) %*% cov %*% diff)
  return(MSE_value)
}

###############################################################################
# SET UP RESULT STORAGE
###############################################################################
# Total simulation runs: (nn) x (missing_props) x (SNR_list) x (seeds) x (methods) x (D imputations)
total_runs <- length(missing_props) * length(SNR_list) * replicate * D
res_cols <- c("seed", "missing_prop", "SNR", "imputation",
              paste0("O", c("Int", paste0("b", 1:p))),
              paste0("Var", paste0("X", 1:p)),
              "MSE", "comp_time")
res <- array(NA, dim = c(total_runs, length(res_cols)),
             dimnames = list(NULL, res_cols))

#### Simulation
counter <- 1
Cov.X.store <- list()
lamvec <- (2^(seq(-1, 4, by = 0.02)))^2 / 2

time_initial <- proc.time()[["elapsed"]]
for (miss_prop in missing_props) {
  for (snr in SNR_list) {
    for (seed in st.seed:(st.seed + replicate - 1)) {
      set.seed(seed)
      generated.data <- generate_and_process_data(seed = seed, n = n, p = p, D = D,
                                                  desired_SNR = snr,
                                                  missing_proportion = miss_prop,
                                                  beta = beta.true)
      Cov.X.store[[length(Cov.X.store) + 1]] <- generated.data$cov.X
      newdata <- generated.data$newdata
      x_list <- generated.data$x
      y_list <- generated.data$y
      
      pf <- rep(1, p)
      adWeight <- rep(1, p)
      if (!exists("cv.galasso")) stop("Function 'cv.galasso' not found. Is 'miselect' loaded correctly?")

      ## fit the model with the optimal lambda
      time_start <- proc.time()
      cv_fit <- cv.galasso(x_list, y_list, adWeight = adWeight, pf = pf) 
      best_lambda_val <- cv_fit$lambda.1se
      raw_coefs <- coef(cv_fit, lambda = best_lambda_val)
      comp_time <- (proc.time() - time_start)[["elapsed"]]
      
      model_fit <- process_galasso_output(raw_coefs, p, D)
      
      for (d in 1:D) {
        res[counter, "seed"] <- seed
        res[counter, "missing_prop"] <- miss_prop
        res[counter, "SNR"] <- snr
        res[counter, "imputation"] <- d
        
        coef_vec <- as.vector(model_fit$coefficients[, d])
        res[counter, 5:(4 + p + 1)] <- coef_vec[1:(p + 1)]
        res[counter, (5 + p + 1):(4 + 2 * p + 1)] <- as.logical(model_fit$varsel)
        
        res[counter, (4 + 2 * p + 2)] <- MSE(Hat.beta = as.numeric(res[counter, 6:(4 + p + 1)]),
                                             beta = beta.true,
                                             cov = as.matrix(generated.data$cov.X))
        res[counter, (4 + 2 * p + 3)] <- comp_time
        
        if ((counter/D) %% 5 == 0) {
          cat(counter/D, "/", replicate * length(missing_props) * length(SNR_list)," completed, Time elapsed:", (proc.time() - time_initial)[["elapsed"]], "\n")
        }
        
        counter <- counter + 1
      }
    }
  }
}

# Define functions for sensitivity and specificity:
SEN <- function(model.pick, true.pick) {
  tp <- sum(true.pick == TRUE & model.pick == TRUE)
  fn <- sum(true.pick == TRUE & model.pick == FALSE)
  if (tp + fn > 0) tp / (tp + fn) else NA
}
SPE <- function(model.pick, true.pick) {
  tn <- sum(true.pick == FALSE & model.pick == FALSE)
  fp <- sum(true.pick == FALSE & model.pick == TRUE)
  if (tn + fp > 0) tn / (tn + fp) else NA
}
MCC <- function(model.pick, true.pick) {
  tp <- sum(true.pick == TRUE  & model.pick == TRUE)
  tn <- sum(true.pick == FALSE & model.pick == FALSE)
  fp <- sum(true.pick == FALSE & model.pick == TRUE)
  fn <- sum(true.pick == TRUE  & model.pick == FALSE)
  
  numerator   <- (tp * tn) - (fp * fn)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  if (denominator > 0) numerator / denominator else NA
}


#### Aggregated results (combined imputations)
total_runs.withoutD <- length(missing_props) * length(SNR_list) * replicate
aggregated_res_cols <- c("missing_prop", "SNR", "seed", 
              paste0("O", c("Int", paste0("b", 1:p))),
              paste0("Var", paste0("X", 1:p)),
              "SEN", "SPE", "MCC",
              "MSE.avg", "comp_time")
aggregated_res <- array(NA, dim = c(total_runs.withoutD, length(aggregated_res_cols)),
             dimnames = list(NULL, aggregated_res_cols))

agg_counter <- 1

for (miss_prop in missing_props) {
  for (snr in SNR_list) {
    for (seed in st.seed:(st.seed + replicate - 1)) {
      aggregated_res[agg_counter, "missing_prop"] <- miss_prop
      aggregated_res[agg_counter, "SNR"] <- snr
      aggregated_res[agg_counter, "seed"] <- seed
      
      aggregated_res[agg_counter, "OInt"] <- mean(res[(D*(agg_counter-1)):D,"OInt"])
      
      for (i in 1:p) {
        aggregated_res[agg_counter, paste0("Ob", i)] <- 
          mean(res[(D*(agg_counter-1)+1):(D*agg_counter), paste0("Ob", i)])
        aggregated_res[agg_counter, paste0("VarX", i)] <- 
          res[agg_counter*D, paste0("VarX", i)]
      }
      
      aggregated_res[agg_counter, "SEN"] <- SEN(as.logical(aggregated_res[agg_counter,paste0("VarX", 1:p)]),beta.check)
      aggregated_res[agg_counter, "SPE"] <- SPE(as.logical(aggregated_res[agg_counter,paste0("VarX", 1:p)]),beta.check)
      aggregated_res[agg_counter, "MCC"] <- MCC(as.logical(aggregated_res[agg_counter,paste0("VarX", 1:p)]),beta.check)
      
      aggregated_res[agg_counter, "MSE.avg"] <- mean(res[(D*(agg_counter-1)):D,"MSE"])
      aggregated_res[agg_counter, "comp_time"] <- mean(res[(D*(agg_counter-1)):D,"comp_time"])
      
      agg_counter <- agg_counter + 1
    }
  }
}

#### Final results (combined seeds)
total_runs.withoutD.withoutS <- length(missing_props) * length(SNR_list)
combined_res_cols <- c("missing_prop", "SNR", 
                         paste0("O", c("Int", paste0("b", 1:p))),
                         "SEN.avg", "SEN.sd",
                         "SPE.avg", "SPE.sd",
                         "MCC.avg", "MCC.sd",
                         "MSE.avg", "MSE.sd",
                          "comp_time.avg", "comp_time.sd")
combined_res <- array(NA, dim = c(total_runs.withoutD.withoutS, length(combined_res_cols)),
                        dimnames = list(NULL, combined_res_cols))

combined_counter <- 1

for (miss_prop in missing_props) {
  for (snr in SNR_list) {
    combined_res[combined_counter, "missing_prop"] <- miss_prop
    combined_res[combined_counter, "SNR"] <- snr
    
    combined_res[combined_counter, "OInt"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"OInt"])
    
    for (i in 1:p) {
      combined_res[combined_counter, paste0("Ob", i)] <- 
        mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter), paste0("Ob", i)])
    }
    
    combined_res[combined_counter, "SEN.avg"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"SEN"])
    combined_res[combined_counter, "SEN.sd"] <- sd(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"SEN"])
    
    combined_res[combined_counter, "SPE.avg"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"SPE"])
    combined_res[combined_counter, "SPE.sd"] <- sd(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"SPE"])
    
    combined_res[combined_counter, "MCC.avg"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"MCC"])
    combined_res[combined_counter, "MCC.sd"] <- sd(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"MCC"])
    
    combined_res[combined_counter, "MSE.avg"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"MSE.avg"])
    combined_res[combined_counter, "MSE.sd"] <- sd(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"MSE.avg"])
    
    combined_res[combined_counter, "comp_time.avg"] <- mean(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"comp_time"])
    combined_res[combined_counter, "comp_time.sd"] <- sd(aggregated_res[((replicate*combined_counter)-replicate+1):(replicate*combined_counter),"comp_time"])
    
    combined_counter <- combined_counter + 1
  }
}

### Exporting these results
# Define today's date and method
today <- format(Sys.Date(), "%Y-%m-%d")
method <- "GLASSO"

# Create directory name
folder_name <- paste0(today, "_p", p, "_", method)

# Create the folder (only if it doesn't already exist)
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

# Convert objects to data frames if needed
res_df <- as.data.frame(res)
aggregated_res_df <- as.data.frame(aggregated_res)
combined_res_df <- as.data.frame(combined_res)

# Build file paths within the new folder
file_res <- file.path(folder_name, paste0("res_", today, "_p", p, "_", method, ".csv"))
file_agg <- file.path(folder_name, paste0("aggregated_res_", today, "_p", p, "_", method, ".csv"))
file_combined <- file.path(folder_name, paste0("combined_res_", today, "_p", p, "_", method, ".csv"))

# Write the CSVs to the folder
write.csv(res_df, file = file_res, row.names = FALSE)
write.csv(aggregated_res_df, file = file_agg, row.names = FALSE)
write.csv(combined_res_df, file = file_combined, row.names = FALSE)

cat("Results written to folder:", normalizePath(folder_name), "\n")
