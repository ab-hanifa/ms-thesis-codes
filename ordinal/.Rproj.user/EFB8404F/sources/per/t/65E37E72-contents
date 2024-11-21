source("ordinal_fn_covarage_profile.R")

set.seed(009)
cluster_size <- c(5, 10, 30, 50)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number, 
                            cluster_size)

data_diff <- parallel::mclapply(1:4000, function(x){
  map2(.x = cl_nmbr_size$cluster_number,
       .y = cl_nmbr_size$cluster_size, 
       ~ sim_ordinal(
         n = .y, # cluster size
         m = .x, # cluster number
         fixed_effects = 0.75,
         threshold = c(0.40, 0.75, 1.25),
         sigma_u = 0.50,
         ran.type = "intercept"
       ))
},
mc.cores = 4
)

ordinal_estimation_diff <- parallel::mclapply(data_diff, function(x) {
  map(x, ~estim_fn_ord(
    data_ord = .x, 
    sigma_u = 0.5,
    ran.type = "intercept")
  )
},
mc.cores = 4)

save(ordinal_estimation_diff, 
     file = "sim_result/coverage-diff/ordinal_estimation_diff.Rdata")

## I took each time 200 elements of the data for 20 times and then took 
##  to simulate due lack of  processors