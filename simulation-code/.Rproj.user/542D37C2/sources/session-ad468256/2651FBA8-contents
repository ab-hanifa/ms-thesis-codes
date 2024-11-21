gen_two_level_slope_data <- function(m, n, fixed_coeff, sigma_mat, data_seed) {
  # m => number of cluster
  # n => size of each cluster
  N = m*n 
  
  set.seed(data_seed)
  
  x1c <- rnorm(N)
  x2 <- rnorm(N)
  x2b <- ifelse(x2 <= 0.5, 0, 1)
  
  u <- MASS::mvrnorm(m, mu = c(0, 0), Sigma = sigma_mat)
  u0j <- rep(u[, 1], each = n)
  u1j <- rep(u[, 2], each = n)
  
  beta0 <- fixed_coeff[1]
  beta1 <- fixed_coeff[2]
  beta2 <- fixed_coeff[3]
  
  eta_ij <- beta0 + u0j + beta1*x1c + u1j*x1c + beta2*x2b 
  pi_ij <- exp(eta_ij) / (1 + exp(eta_ij))
  
  yij <- rbinom(N, size = 1, prob = pi_ij)
  
  multi_data <- data.frame(
    cluster = rep(seq(m), each = n),
    id = rep(seq(n), times = m),
    u0j = u0j,
    u1j = u1j,
    X1c = x1c,
    X2b= x2b,
    Yij = yij
  )
  return(multi_data)
}


calc_two_lvl_slp_mor_stats <- function(sigma_mat, sigma_u_sq_hat,
                                       var_sigma_u_sq_hat, x1_val) {
  
  # sigma mat components -----------
  sigma_u1_sq <- sigma_mat[1, 1]
  sigma_u12_sq <- sigma_mat[1, 2]
  sigma_u2_sq <- sigma_mat[2, 2]
  
  # calc mor -----------------------
  mor_hat <- exp(
    sqrt(2*sigma_u_sq_hat[1] + 2*x1_val^2*sigma_u_sq_hat[3] + 4*x1_val*sigma_u_sq_hat[2]) 
    * qnorm(0.75)
  )
  log_mor_hat <- log(mor_hat)
  
  
  # delta method -------------------
  log_mor_expr <- function(x, x1_val) {
    # here ran-effect is parameterized as function of variance term
    # instead of standard dev.
    sqrt((2*x[1]) + (2*x1_val^2*x[3]) + (4*x1_val*x[2])) * qnorm(0.75)
  }
  
  J <- numDeriv::jacobian(log_mor_expr, sigma_u_sq_hat, x1_val = x1_val)
  log_se_mor_hat <- sqrt(as.numeric(J %*% var_sigma_u_sq_hat %*% t(J)))
  se_mor_hat <- exp(log_se_mor_hat)
  
  # coverage calc ------------------
  ci <- log_mor_hat + c(-1, 1) * 1.96 * log_se_mor_hat
  ci_exp <- exp(ci)
  
  true_mor <- exp(
    sqrt(2*sigma_u1_sq + 2*x1_val^2*sigma_u2_sq + 4*x1_val*sigma_u12_sq) * qnorm(0.75)
  )
  coverage <- as.numeric(ci_exp[1] <= true_mor && ci_exp[2] >= true_mor)
  
  rel_bias <- (mor_hat - true_mor) / true_mor
  
  return(
    c(true_mor = true_mor, mor_hat = mor_hat, se_mor_hat = se_mor_hat, 
      coverage = coverage, rel_bias = rel_bias)
  )
  
}



est_two_lvl_slope_mor <- function(m, n, fixed_coeff, sigma_mat, data_seed) {
  
  # data generation ----------------
  multi_data_slope <- gen_two_level_slope_data(m, n, fixed_coeff,
                                               sigma_mat, data_seed)
  
  # sigma mat components -----------
  sigma_u1_sq <- sigma_mat[1, 1]
  sigma_u12_sq <- sigma_mat[1, 2]
  sigma_u2_sq <- sigma_mat[2, 2]
  
  # prevalence ---------------------
  prevalence <- mean(multi_data_slope$Yij, na.rm = TRUE)
  
  # model fitting ------------------
  multi_model_slope <- mixed_model(fixed = Yij ~ X1c + X2b, 
                                   random = ~ X1c | cluster, 
                                   data = multi_data_slope, 
                                   nAGQ = 11, family = binomial("logit"))
  
  # extracting value from model ----
  fixed_effect <- unname(multi_model_slope$coefficients)
  beta0_hat <- fixed_effect[1]
  beta1_hat <- fixed_effect[2]
  beta2_hat <- fixed_effect[3]
  sigma_u_sq_hat <- sym_mat_to_vec(multi_model_slope$D) # D_11 D_12 D_22
  sigma_u1_sq_hat <- sigma_u_sq_hat[1]
  sigma_u12_sq_hat <- sigma_u_sq_hat[2]
  sigma_u2_sq_hat <- sigma_u_sq_hat[3]
  
  is_model_converged <- multi_model_slope$converged
  
  # gets variance of rand effect
  # D_11 D_12 D_22 
  var_sigma_u_sq_hat <- var_slope_rand_effect(multi_model_slope) 
  
  # mor calc. ----------------------
  x1_q1 <- unnamed_quantile(multi_data_slope$X1c, probs = 0.25)
  x1_q2 <- mean(multi_data_slope$X1c)
  x1_q3 <- unnamed_quantile(multi_data_slope$X1c, probs = 0.75)
  
  mor_stats_q1 <- calc_two_lvl_slp_mor_stats(sigma_mat, sigma_u_sq_hat,
                                             var_sigma_u_sq_hat, x1_q1)
  names(mor_stats_q1) <- paste0(names(mor_stats_q1), "_q1")
  
  mor_stats_q2 <- calc_two_lvl_slp_mor_stats(sigma_mat, sigma_u_sq_hat,
                                             var_sigma_u_sq_hat, x1_q2)
  names(mor_stats_q2) <- paste0(names(mor_stats_q2), "_q2")
  
  mor_stats_q3 <- calc_two_lvl_slp_mor_stats(sigma_mat, sigma_u_sq_hat,
                                             var_sigma_u_sq_hat, x1_q3)
  names(mor_stats_q3) <- paste0(names(mor_stats_q3), "_q3")
  
  # returned values by above three is 
  # true_mor, mor_hat, se_mor_hat, coverage, rel_bias
  
  out_vec <- c(mor_stats_q1, mor_stats_q2, mor_stats_q3,
              sigma_u1_sq_hat = sigma_u1_sq_hat,
              sigma_u12_sq_hat = sigma_u12_sq_hat,
              sigma_u2_sq_hat = sigma_u2_sq_hat,
              beta0_hat = beta0_hat, beta1_hat = beta1_hat,
              beta2_hat = beta2_hat,
              prevalence = prevalence,
              converged = is_model_converged)
  # ***reminder `converged` needs to be the last column****
  
  # getting mean_mor and mean_se_mor to catch unstable estimate
  mean_mor_hat <- mean(
    c(mor_stats_q1[2], mor_stats_q2[2], mor_stats_q3[2]), na.rm = FALSE
      )
  
  mean_se_mor_hat <- mean(
    c(mor_stats_q1[3], mor_stats_q2[3], mor_stats_q3[3]), na.rm = FALSE
  )
  
  if(is.na(is_model_converged) || !is_model_converged || mean_mor_hat > 20 || 
     mean_se_mor_hat > 20) {
    out_vec_names <- names(out_vec)
    out_vec <- c(rep(NA, length(out_vec) - 1), FALSE)
    names(out_vec) <- out_vec_names
  }
  
  return(out_vec)
}


simulate_two_lvl_slope <- function(m, n, fixed_coeff, sigma_mat, nsims = 1000, 
                           log_file, seed, more_iter, ...) {

  # creating extra sims to get nsims after accounting 
  # for non-converged cases ---------------------
  total_sims = nsims + more_iter
  conv_case_num = 0
  runs_required = 0
  
  # creating placeholder matrix for result --------
  out_mat <- matrix(NA, nrow = total_sims, ncol = 23)
  
  # mor stats names -----
  mor_stats_names <- c("true_mor", "mor_hat", "se_mor_hat", "coverage","rel_bias")
  mor_stats_names_q1 <- paste0(mor_stats_names, "_q1")
  mor_stats_names_q2 <- paste0(mor_stats_names, "_q2")
  mor_stats_names_q3 <- paste0(mor_stats_names, "_q3")
  
  out_colnames <- c(mor_stats_names_q1, mor_stats_names_q2,
                    mor_stats_names_q3,
                    "sigma_u1_sq_hat", "sigma_u12_sq_hat", "sigma_u2_sq_hat",
                    "beta0_hat", "beta1_hat", "beta2_hat", 
                    "prevalence", "converged")
  colnames(out_mat) <- out_colnames
  
  # set the random state initiating seed
  set.seed(seed)
  iteration_seeds <- sample(x = 10:1000000, size = total_sims, replace = FALSE)
  # starting simulation ---------------------------
  for (i in 1:total_sims) {
    
    # checking for converged nsims ----------------
    if (conv_case_num == nsims) { 
      break
    } else {
      runs_required = runs_required + 1
    }
    
    data_seed <- iteration_seeds[i]
    model_output <- safe_and_quietly(fun = est_two_lvl_slope_mor,
                                     m = m, n = n, 
                                     fixed_coeff = fixed_coeff,
                                     sigma_mat = sigma_mat, 
                                     data_seed = data_seed)
    
    # extraction ---------------------------------
    model_result <- model_output$result
    model_messages <- model_output$messages
    model_warnings <- model_output$warnings
    model_error <- model_output$error
    
    # account for error -------------------------
    msg_prob_detected <- is_valid_str(paste0(model_messages, collapse = ""))
    warning_detected <- is_valid_str(paste0(model_warnings, collapse = ""))
    error_detected <- is_valid_str(paste0(model_error, collapse = ""))
    
    if(is_na_result(model_result) || 
       msg_prob_detected || warning_detected || error_detected) {
      # setting model_result to be vec of NA in case of 
      # non-convergence
      model_result <- c(rep(NA, ncol(out_mat) - 1), 0)
      names(model_result) <- out_colnames
    }
    
    model_convergence <- as.logical(model_result["converged"])
    
    # checking total number of converged cases ----
    if(model_convergence) {
      conv_case_num = conv_case_num + 1
    }
    
    # logging simulation run ----------------------
    log_sim_run(convergence = model_convergence, message = model_messages, 
                warning = model_warnings, error = model_error,
                data_seed = data_seed, log_file = log_file, iter_no = i)
    
    # storing sim-run ---------------------------
    out_mat[i, ] <- model_result
  }
  
  return(
    list(
      out_mat = out_mat,
      runs_required = runs_required
    )
  )
}

