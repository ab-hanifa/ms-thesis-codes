library(tidyverse)
# Create a vector of file names
file_names <- paste0("cat_glmadap_all_high_", 1:4, ".Rdata")

# Use lapply to load each file into the global environment
sapply(file_names, load, envir = .GlobalEnv)

data_list_high_fn <- function(dt_list){
  
  inner_function <- function(dt_list_inner){
  # Function to check for errors
  is_error <- function(x) inherits(x, "try-error")
  
  # Use try() to catch errors and exclude them
  clean_list <- lapply(dt_list_inner, function(x) try(x, silent = TRUE))
  filtered_list <- Filter(function(x) !is_error(x), dt_list_inner)
  
  data_mlt_lst_new <- lapply(filtered_list, function(y) {
    Filter(function(df) !any(is.na(df)),
           lapply(1:500, function(x) y[[x]]))})
  
  data_mlt <- lapply(data_mlt_lst_new, function(.x) .x |> list_rbind())
  
  return(data_mlt)
  
  }
  #inner_function(cat_glmadap_all_high_2)
  dtl <- lapply(dt_list, inner_function) 
  dtlst <- lapply(1:30, function(x){
    lapply(dtl, function(y) y[[x]]) |> list_rbind()
  })
  return(dtlst)
}

high_model_dt <- data_list_high_fn(list(cat_glmadap_all_high_1, 
                                        cat_glmadap_all_high_2, 
                                        cat_glmadap_all_high_3,
                                        cat_glmadap_all_high_4))


cluster_size <- c(10, 20, 30, 50)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number,
                            cluster_size)
index <- c(1:30)[-seq(5, 30, by = 5)]
high_model_result <- lapply(index, function(x) high_model_dt[[x]])
## Exploratory analysis for sigma square

plots <- lapply(1:24, function(x) {
  high_model_result[[x]]|> 
      # filter(sigma_2_u_hat > 0.10 &
      #          sigma_2_u_hat < quantile(sigma_2_u_hat, probs = 0.75) +
      #          1.5*IQR(sigma_2_u_hat),
      #        .by = model) |>
    ggplot(
      mapping = aes(x = MOR_hat,
                    color = model)
           )+
    geom_histogram()+ 
    xlab(label = paste("M = ",
                       cl_nmbr_size$cluster_number[x],
                       "N = ", cl_nmbr_size$cluster_size[x]))+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 6))
  }
)

library(gridExtra)
grid.arrange(grobs = plots, ncol = 6)

########################################################################
sum_fn <- function(dt){
  # dt <- x |>
  #   filter(MOR_hat > 1.01 &
  #            MOR_hat < quantile(MOR_hat, probs = 0.75) + 1.5*IQR(MOR_hat),
  #          .by = model)
  
  if(nrow(dt) < 2000) data_m <- dt
  else data_m <- dt |> slice(1:1000, .by = model)
  
  data_m |> 
    mutate(sigma2 = if_else(model == "Model1", 0.55, 0.75)) |> 
    reframe(
      m_beta_0 = mean(beta_0),
      m_beta_1 = mean(beta_1),
      m_beta_2 = mean(beta_2),
      sigma_2_hat = mean(sigma_2_u_hat),
      mor_hat = mean(MOR_hat),
     rb = (mean(MOR_hat) - mean(true_mor))*100/mean(true_mor),
     sd_mor = sd(MOR_hat),
      m_coverage = mean(coverage),
      .by = "model"
    ) |> 
    mutate(iter = if_else(nrow(dt) < 2000, round(nrow(dt)/2), 1000))
}

cluster_size <- c(10, 20, 30, 50, 100)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number,
                            cluster_size)

sim_sum_mult_glmadap_high <- map(high_model_dt, ~sum_fn(.x)) |> 
  list_rbind() |> 
  arrange(model) |> 
  bind_cols(cl_nmbr_size |>
              bind_rows(cl_nmbr_size)) |> 
  select(model, cluster_number, 
         cluster_size, everything())|> 
  filter(cluster_size != 100) #|>
# mutate(Iter = rep(#1000, 32)
#   sapply(1:16, function(x) 4000 - length(data_mlt_lst_new[[x]])), 2)
#        )

sim_sum_mult_glmadap_high |> view()
# 
save(sim_sum_mult_glmadap_high, file = "../sim_sum_mult_glmadap_high.Rdata")
# 
# 
# calculate_iqr_bounds <- function(x) {
#   # Remove NA values
#   x <- na.omit(x)
#   
#   # Calculate the quartiles
#   q1 <- quantile(x, 0.25)
#   q3 <- quantile(x, 0.75)
#   
#   # Calculate the IQR
#   iqr_val <- q3 - q1
#   
#   # Calculate the lower and upper bounds
#   lower_bound <- q1 - 1.5 * iqr_val
#   upper_bound <- q3 + 1.5 * iqr_val
#   
#   return(c(lower_bound = lower_bound, upper_bound = upper_bound))
# }
# 
# high_model_dt[[1]] |> 
#   reframe(
#     lt = calculate_iqr_bounds(sigma_2_u_hat),
#     .by = model)

###################### Descriptive Stat Multinomial ##########################
cluster_size <- c(10, 20, 30, 50, 100)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number,
                            cluster_size)

descriptive_glmadap_high <- map(high_model_dt, ~ .x |> 
                                   filter(sigma_2_u_hat >= 0.09 & 
                                            sigma_2_u_hat <= quantile(sigma_2_u_hat, 
                                                                      probs = 0.75,
                                                                      na.rm = TRUE) + 
                                            1.5*IQR(sigma_2_u_hat, 
                                                    na.rm = TRUE)) |> 
                                   reframe(min = min(sigma_2_u_hat,
                                                     na.rm = TRUE),
                                           max = max(sigma_2_u_hat,
                                                     na.rm = TRUE),
                                           med = median(sigma_2_u_hat,
                                                        na.rm = TRUE),
                                           mean = mean(sigma_2_u_hat, 
                                                       na.rm = TRUE),
                                           sd_sigma2_hat = sd(sigma_2_u_hat,
                                                              na.rm = TRUE),
                                           .by = model)) |>
  list_rbind() |> 
  arrange(model) |> 
  bind_cols(cl_nmbr_size |>
              bind_rows(cl_nmbr_size)) |> 
  select(model, cluster_number, 
         cluster_size, everything())|> 
  filter(cluster_size != 100)

save(descriptive_glmadap_high, file = "../descriptive_glmadap_high.Rdata")
