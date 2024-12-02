library(tidyverse)
library(ordinal)
set.seed(009)

sim_ordinal <- function(n, m, fixed_effects, 
                        threshold, sigma_u, #seed,
                        ran.type){
  
  
  # ## Setting seed
  # set.seed(seed = seed)
  N <- m*n
  x <- stats::rnorm(N, mean = 2, sd = 1.5)
  
  if(any(sigma_u < 0) || n < 0 || m < 0){
    stop("n, m or sigma_u have to be positive numbers")
  }
  
  if(ran.type == "slope" && length(sigma_u) < 2){
    stop("sigma_u have to be greater than 1")
  }
  
  if(ran.type == "intercept"){
    u0j <- rep(stats::rnorm(m, mean = 0, sd = sigma_u[1]), each = n)
    x_mat <- cbind(x, u0j)
    xb <- as.vector(x_mat %*% t(cbind(fixed_effects, 1)))
  }
  else if(ran.type == "slope") {
    u0j <- rep(stats::rnorm(m, mean = 0, sd = sigma_u[1]), each = n)
    u1j <- rep(stats::rnorm(m, mean = 0, sd = sigma_u[2]), each = n)
    x_mat <- cbind(x, u0j, u1j*x)
    xb <-  as.vector(x_mat %*% t(cbind(fixed_effects, 1, 1)))
  }
  else {
    stop("`ran.type` either be 'intercept' or 'slope'")
  }
  
  cumulative_prob <- sapply(threshold, function(p){
    1/(1 + exp(-(p - xb)))
  })
  
  cat_prob <- cbind(cumulative_prob[, 1],
                    cumulative_prob[, -1] - 
                      cumulative_prob[, -ncol(cumulative_prob)],
                    1 - cumulative_prob[, ncol(cumulative_prob)])
  
  y <- apply(cat_prob, 1, FUN = function(c){
    which.max(cumsum(c) >= runif(1))
  })
  
  dat <- dplyr::tibble(
    cluster = as.factor(rep(seq(m), each = n)),
    cluster_size = rep(seq(n), times = m),
    x = x,
    ordinal_y = factor(y)
  )
  
  return(dat)
  
}

#############################################
estim_fn_ord <- function(data_ord, sigma_u, 
                         ran.type){
  
  #dat1 <- {{data_ord}}
  # 
  # prevalence <- dat1 |> 
  #   dplyr::count(ordinal_y) |> 
  #   dplyr::mutate(prop = n/sum(n))
  
  if(ran.type == "intercept"){
    mod1 <- ordinal::clmm2(ordinal_y ~ x,
                           random = cluster,
                           data = data_ord,
                           Hess = TRUE, 
                           control = clmm2.control(
                             method = "nlminb",
                             maxIter = 100, 
                             gradTol = 1e-3),
                           threshold = "flexible",
                           link = "logistic")
    
    sd_intercept <- mod1$stDev[[1]]
   # mod1$data_ord <- data_ord
    is_converged <- mod1$convergence
    
    if(!is_converged || det(mod1$Hessian) <= 0 || sd_intercept < 0.08){
      outcome_vec <- NA
    }
    
    else {
      mor_hat <- exp(sqrt(2)*sd_intercept*qnorm(0.75))
      
      true_mor <- exp(sqrt(2)*sigma_u*qnorm(0.75))
      
      assign("data_ord", data_ord, envir = globalenv())
      pfr_mod <- profile(mod1, nSteps = 5)
      ci_sd <- as.vector(confint(pfr_mod))
      
      # Optionally clean up after profiling
      rm("data_ord", envir = globalenv())
      
      ci_mor <- exp(ci_sd*sqrt(2)*qnorm(0.75)) 
      
      coverage <- 1*(ci_mor[1] < true_mor && ci_mor[2] > true_mor)
      outcome_vec <- c(mor_hat = mor_hat, 
                       true_mor = true_mor,
                       sd = sd_intercept,
                       mod1$coefficients[-5], 
                       coverage = coverage)
    }
  }
  
  else{
    mod1 <- ordinal::clmm(ordinal_y ~ x + (1 + x | cluster),
                          data = data_ord, nAGQ = 7,
                          link = "logit")
    
    outcome_vec <- c(sd1 = mod1$ST[[1]][1],
                     sd2 = mod1$ST[[1]][2],
                     mod1$coefficients)
  }
  
  return(outcome_vec)#list(prevalence, outcome_vec))
  
}
