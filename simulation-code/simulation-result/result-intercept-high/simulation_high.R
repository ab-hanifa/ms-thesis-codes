
source("intercept-multinom-function-n.R")
set.seed(1009)

data_high <- parallel::mclapply(1:3000, function(x){
  map2(.x = cl_nmbr_size$cluster_number,
       .y = cl_nmbr_size$cluster_size, 
       ~ data_sim_fn_int(m = .x, n = .y, # cluster number and size
                         fixed_coef = cbind(c(-0.80, 1.46, 0.73),
                                            c(-0.95, 1.60, 0.90)),
                         sigma2 = c(0.55, 0.75), 
                         reff = "diff"
       ))
  }, 
mc.cores = 5)



cat_glmadap_all_high <- parallel::mclapply(1:24, function(p) {
  lapply(data_high, function(x) {
    estimation_fn_int(datx = x[[p]], sigma2 = c(0.55, 0.75))
  })
},
mc.cores = 4)

save(cat_glmadap_all_high, 
     file = "cat_glmadap_all_high.Rdata")

#### I used four replication due lack of processor and each time I took
#### 750 elements of the data list