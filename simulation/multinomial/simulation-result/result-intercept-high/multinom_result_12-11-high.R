library(tidyverse)
load("simulation-result/result-intercept-high/cat_glmadap_all_high_1.Rdata")
load("simulation-result/result-intercept-high/cat_glmadap_all_high_2.Rdata")
load("simulation-result/result-intercept-high/cat_glmadap_all_high_3.Rdata")
load("simulation-result/result-intercept-high/cat_glmadap_all_high_4.Rdata")

combine_and_filter_data <- function(list_of_data, indices = 1:24) {
  map(indices, function(index) {
    map(list_of_data, function(x) {
      Filter(function(df) !any(is.na(df)), x[[index]]) |> 
        list_rbind()
    }) |> 
      bind_rows() 
  })
}

data_lists <- list(cat_glmadap_all_high_1,
                   cat_glmadap_all_high_2,
                   cat_glmadap_all_high_3,
                   cat_glmadap_all_high_4)
dat_list <- combine_and_filter_data(data_lists)

sum_fn <- function(x){
  dt <- x |>
    filter(MOR_hat > 1.1 &
             MOR_hat < quantile(MOR_hat, probs = 0.75) + 1.5*IQR(MOR_hat),
           .by = model)

  if(nrow(dt) < 2000) data_m <- dt
  else data_m <- dt |> slice(1:1000, .by = model)
  
  data_m |> 
    mutate(sigma = if_else(model == "Model1", 0.55, 0.75)) |> 
    reframe(
      m_beta_0 = mean(beta_0),
      m_beta_1 = mean(beta_1),
      m_beta_2 = mean(beta_2),
      sigma_2 = mean(sigma_2_u_hat),
      mor_hat = mean(MOR_hat),
      rb = (mean(MOR_hat) - mean(true_mor))*100/mean(true_mor),
      se_mor = sd(MOR_hat),
      m_coverage = mean(coverage),
      .by = "model"
    ) |> 
    mutate(iter = round(3000 - nrow(dt)/2 + 1000))
}

cluster_size <- c(10, 20, 30, 50)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number,
                            cluster_size)

multinom_mod_high <- map(dat_list, sum_fn) |> 
  list_rbind() |>
  arrange(model) |> 
  bind_cols(cl_nmbr_size |> 
            bind_rows(cl_nmbr_size)) |> 
  select(cluster_number,
         cluster_size,
         everything())

save(multinom_mod_high, 
     file = "simulation-result/result-intercept-high/multinom_mod_high.Rdata")

plots <- lapply(1:24, function(x) {
  dat_list[[x]]|> 
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

