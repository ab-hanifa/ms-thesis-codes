library(tidyverse)
library(tools)


dat_ord_same_fn <- function(file_lst, path = "ordinal/sim_result/coverage-same/") {
  map(1:24, function(index) { 
    map_dfr(file_lst, function(file) {
      
      dat_list <- get(load(file.path(path, paste0(file, ".Rdata"))))
      
      map_dfr(dat_list, function(x) {
        column_names <- c("mor_hat", "true_mor", "sd",
                          "1|2", "2|3", "3|4", "x", "coverage")
        
        
        if (!is.list(x) || is.null(x[[index]]) || all(is.na(x[[index]]))) {
          return(setNames(rep(NA, length(column_names)), column_names))
        } else {
          
          as_tibble(t(x[[index]])) 
        }
      }, .id = "source_file") |>  
        rename(
          "alpha_h1" = `1|2`,
          "alpha_h2" = `2|3`,
          "alpha_h3" = `3|4`,
          "beta" = x
        )
    })
  })
}

# Get list of file names (without extensions) in the directory
file_lst <- file_path_sans_ext(
  list.files(path = "ordinal/sim_result/coverage-same/", 
             pattern = "\\.Rdata$", 
             full.names = FALSE)
)

# Apply the function
same_cover_dat <- dat_ord_same_fn(file_lst)

cluster_size <- c(5, 10, 30, 50)
cluster_number <- c(10, 30, 50, 100, 150, 200)
cl_nmbr_size <- expand_grid(cluster_number, 
                            cluster_size)

true_mor <- exp(sqrt(2*0.50^2)*qnorm(0.75))
thresholds <- c(0.50, 1, 1.50)
thresholds_diff <- c(0.40, 0.75, 1.25)
true_beta <- 0.75
sigma <- 0.50

rb_fn <- function(x){
  filtered_dt <-  x |> 
    drop_na(mor_hat) #|> 
  # filter(sd >= 0.09 & sd <= quantile(sd, probs = 0.75,
  #                                    na.rm = TRUE) + 
  #          1.5*IQR(sd, na.rm = TRUE))
  iter <- nrow(x) - nrow(filtered_dt) + 1000
  
  if(nrow(filtered_dt) > 1000){
    filtered_dt <- filtered_dt |> 
      dplyr::slice(1:1000)
  }
  filtered_dt |> 
    reframe(
      alpha1 = mean(alpha_h1, na.rm = TRUE),
      alpha2 = mean(alpha_h2, na.rm = TRUE),
      alpha3 = mean(alpha_h3, na.rm = TRUE),
      beta = mean(beta, na.rm = TRUE),
      sigma = mean(sd, na.rm = TRUE),
      mean_mor = mean(mor_hat, na.rm = TRUE),
      se_mor = sd(mor_hat, na.rm = TRUE),
      rb_mor = (mean(mor_hat) - mean(true_mor))/mean(true_mor)*100, 
      m_coverage = mean(coverage, na.rm = TRUE),
      iter = iter
    )
}

ordi_rb_dat <- map(same_cover_dat, rb_fn) |> 
  list_rbind() |> 
  bind_cols(cl_nmbr_size) |> 
  select(cluster_number, 
         cluster_size,
         everything())

ordi_rb_dat |> view()
save(ordi_rb_dat, file = "ordi_rb_dat.Rdata")


plots <- lapply(1:24, function(x) {
  same_cover_dat[[x]]|> 
    ggplot(
      mapping = aes(x = mor_hat))+
    geom_histogram(col = "steelblue")+ 
    xlab(label = paste("M = ",
                       cl_nmbr_size$cluster_number[x],
                       "N = ", cl_nmbr_size$cluster_size[x]))+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 6))
}
)

grid.arrange(grobs = plots, ncol = 6)


library(gridExtra)
