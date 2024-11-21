source(here::here("exploration.R"))
library(GLMMadaptive)
library(tidyverse)
set.seed(009)

data_sim_fn_int <- function(m, n, fixed_coef, #seed, 
                            sigma2, reff){
  
  
  N <- m*n
  #set.seed(seed = seed)
  x1c <- rnorm(N, mean = 1, sd = 2)
  x2c <- rnorm(N, mean = 0, sd = 1)
  x2b <- ifelse(x2c <= 0.5, 1, 0)
  
  if(reff == "same"){
    uj <- rep(rnorm(m, mean = 0, sd = sqrt(sigma2)), each = n)
    xb_mat <- cbind(1, x1c, x2b, uj)
    etai <- xb_mat %*% rbind(fixed_coef, 1)
  }
  else {
    if(length(sigma2) < 2){
      stop("When `reff != 'same'` sigma \n
           has to be gearter than or equal \n
           to two and positive numbers \n")
    }
    uj1 <- rep(rnorm(m, mean = 0, sd = sqrt(sigma2[1])), each = n)
    uj2 <- rep(rnorm(m, mean = 0, sd = sqrt(sigma2[2])), each = n)
    xb_mat <- cbind(1, x1c, x2b, uj1, uj2)
    etai <- xb_mat %*% rbind(fixed_coef, diag(2))
  }
  
  pi_hat12 <- exp(etai)/(1 + rowSums(exp(etai)))
  
  pi_hat3 <- 1 - rowSums(pi_hat12)
  pi_hat <- cbind(pi_hat12, pi_hat3)
  
  y <- apply(
    pi_hat,
    1, 
    FUN = function(p) which(rmultinom(n = 1, 
                                      size = 1,
                                      prob = p) == 1)
  )
  
  df <- tibble::tibble(
    cluster = rep(seq(m), each = n),
    cluster_size = rep(seq(n), times = m),
    x1c = x1c,
    x2b = x2b,
    yij = relevel(factor(y), ref = "3")
  )
  return(df)
}


estimation_fn_int <- function(datx, sigma2){
  data_gen_int <- {{datx}}
  
  # prevalence <- data_gen_int |> 
  #   dplyr::count(yij) |> 
  #   dplyr::mutate(proportion = n/sum(n))
  # 
  # require(GLMMadaptive)
  
  mod_stat <- sapply(1:2, function(x){
    filt_dat <- c(2, 1)
    mod1 <- mixed_model(
      fixed = I(yij == x) ~ x1c + x2b,
      random = ~ 1 | cluster,
      data = data_gen_int |> dplyr::filter(yij != filt_dat[x]),
      family = binomial("logit"),
      control = list(iter_EM = 0, 
                     max_coef_value = 1000, 
                     tol1 = 1e-5, tol2 = 1e-5)
    )
    
    fixed_efct <- unname(fixef(mod1))
    beta_0 <- fixed_efct[1]
    beta_1 <- fixed_efct[2]
    beta_2 <- fixed_efct[3]
    
    sigma_2_u_hat <- mod1$D[[1]] # random effect parameter
    ### variance of random effect
    var_sigma_u_2_hat <- vcov_orig_scale(mod1)
    
    ## convergence
    is_converged <- mod1$converged
    
    ## MOR1 calculation
    mor_hat <- exp(sqrt(2 * sigma_2_u_hat) * qnorm(0.75))
    
    ## Log mor to find ci
    log_mor <- log(mor_hat)
    
    # delta method 
    log_mor_int_expr <- function(x) {
      sqrt(2 * x) * qnorm(0.75)
    }
    
    J <- numDeriv::jacobian(log_mor_int_expr, x = sigma_2_u_hat)
    log_se_mor <- as.numeric(sqrt(t(J) %*% var_sigma_u_2_hat %*% J))
    se_mor_ht <- exp(log_se_mor)
    # coverage calculation
    ci <- log_mor + c(-1, 1) * 1.96 * log_se_mor
    ci_exp_1 <- exp(ci)
    true_mor <- exp(sqrt(2 * sigma2[x]) * qnorm(0.75))
    coverage <- 1*(ci_exp_1[1] < true_mor && ci_exp_1[2] > true_mor)
    
    ## Handling NA in converged or mor values
    is_converged <- ifelse(is.na(is_converged), FALSE, is_converged)
    mor_hat <- ifelse(is.na(mor_hat), 0, mor_hat)
    se_mor_ht <- ifelse(is.na(se_mor_ht), 0, se_mor_ht)
    
    outcome_vector <- c(
      MOR_hat = mor_hat,
      true_mor = true_mor,
      sigma_2_u_hat = sigma_2_u_hat,
      beta_0 = beta_0,
      beta_1 = beta_1,
      beta_2 = beta_2,
      coverage = coverage
    )
    
    if(!is_converged || mor_hat > 5.5 || se_mor_ht > 5 ) {
      out_vec_names <- names(outcome_vector)
      outcome_vector <- c(rep(NA, length(outcome_vector) - 1), FALSE)
      names(outcome_vector) <- out_vec_names
    }
    
    return(outcome_vector)
  }
  )
  return(as_tibble(t(mod_stat)) |>
           mutate(model = c("Model1", "Model2")))
}
