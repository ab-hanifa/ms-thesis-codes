mat_to_vec <- function(x) x[upper.tri(x, diag = TRUE)]

#mat_to_vec(mod_glm1$D)


D_chol_to_D <- function(x) {
  # transform the log-cholesky parameterized 
  # value back to original scale
  D <- GLMMadaptive:::chol_transf(x)
  D[upper.tri(D, diag = TRUE)]
}

#D_chol_to_D(mod_glm1$D)

## Converting the log cholesky factor to original variance format
vcov_orig_scale <- function(model) {
  D <- model$D
  
  # transform from covariance matrix to entries of 
  # cholesky factor with log-transformed main diagonal
  D_chol_entries <- GLMMadaptive:::chol_transf(D)
  V_chol <- vcov(model, parm = "var-cov")
  
  J <- numDeriv::jacobian(D_chol_to_D, D_chol_entries)
  # delta method (estimated covariance matrix of entries of D)
  V <- J %*% V_chol %*% t(J)
  
  colnames(V) <- colnames(V_chol)
  rownames(V) <- rownames(V_chol)

  return(V)
}

#vcov_orig_scale(model = mod_glm1)

var_slope_rand_effect <- function(model) {
  rand_var <- vcov_orig_scale(model)
  # v <- diag(rand_var)
  v_mat <- diag(diag(rand_var))
  rownames(v_mat) <- c("D_11", "D_12", "D_22")
  colnames(v_mat) <- c("D_11", "D_12", "D_22")
  return(v_mat)
}

#var_slope_rand_effect(mod_glm1)
