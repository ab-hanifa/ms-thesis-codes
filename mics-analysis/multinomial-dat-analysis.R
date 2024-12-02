library(tidyverse)
library(haven)
library(foreign)
library(GLMMadaptive)
library(gt)
library(broom.mixed)
library(gtsummary)
library(sjPlot)

#view(health_seek_data)
load("mics-analysis/health_seek_data.Rdata")

multinom_summary_tab <- health_seek_data |> 
  select(melevel, HH6, health_care,
         windex5, HH7, #ethnicity,
         HH6) |> 
  tbl_summary(
    by = health_care,
    percent = "row",
    missing = "no",
    label = list(
      melevel ~ "Mother's Education Level"
    )
  ) |> 
  add_p() |> 
  add_overall() 

save(multinom_summary_tab, file = "mics-analysis/multinom_summary_tab.Rdata")


##
health_seek_data |> 
  count(health_care)

source("../ms-thesis/simulation-code/intercept-codes/exploration.R")
ci_mor_fn <- function(mod1, mor_hat){
  var_sigma_u_2_hat <- vcov_orig_scale(mod1)
  sigma_2_u_hat <- mod1$D[[1]]
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
  return(ci_exp_1)
}

#####
model1_null <- mixed_model(fixed = I(health_care == "govt") ~ 1,
                           random = ~ 1|HH7A,
                           family = binomial(),
                           data = health_seek_data |> 
                             filter(health_care != "private"))

null_mor1 <- exp(sqrt(2 * as.vector(model1_null$D)) * qnorm(0.75))
null_mor1
null_ci_mor1 <- ci_mor_fn(model1_null, null_mor1)
null_ci_mor1
##################
model2_null <- mixed_model(fixed = I(health_care == "private") ~ 1,
                           random = ~ 1|HH7A,
                           family = binomial(),
                           data = health_seek_data |> 
                             filter(health_care != "govt"))

null_mor2 <- exp(sqrt(2 * as.vector(model2_null$D)) * qnorm(0.75))
null_mor2
null_ci_mor2 <- ci_mor_fn(model2_null, null_mor2)
null_ci_mor2
##############

model1 <- mixed_model(fixed = I(health_care == "govt") ~ melevel + HH6 + 
                        windex5,
                      random = ~ 1|HH7A,
                      family = binomial(),
                      data = health_seek_data |> 
                        filter(health_care != "private"))

mor1 <- exp(sqrt(2 * as.vector(model1$D)) * qnorm(0.75))
mor1
ci_mor1 <- ci_mor_fn(model1, mor1)
ci_mor1

mod1_summ <- broom.mixed::tidy(model1, exponentiate = TRUE) 
save(mod1_summ, file = "mics-analysis/mod1_summ.Rdata")

confint(model1)
vcov(model1)
summary(model1)
################
model2 <- mixed_model(fixed = I(health_care == "private") ~ 
                        melevel + HH6 + windex5,
                      random = ~ 1|HH7A,
                      family = binomial(),
                      data = health_seek_data |> 
                        filter(health_care != "govt"))


mor2 <- exp(sqrt(2 * as.vector(model2$D)) * qnorm(0.75))
mor2
ci_mor2 <- ci_mor_fn(model2, mor2)
ci_mor2

mod2_summ <- broom.mixed::tidy(model2, exponentiate = TRUE) 

save(mod2_summ, file = "mics-analysis/mod2_summ.Rdata")

summary(model2)

library(sjPlot)
broom.mixed::tidy(model1)

library(ggplot2)

health_seek_data |> 
  select(health_care, HH7) |> 
  count(
    health_care, division = HH7
  ) |> 
  mutate(
    division = fct_reorder(division, n, .desc = TRUE)  
  ) |>
  ggplot(
    mapping = aes(
      x = division,
      y = n, 
      fill = health_care
    )
  ) +
  geom_col(
    position = position_stack(),
    width = 0.90, 
    color = "black") +  
  scale_fill_manual(values = c("steelblue",
                               "coral", 
                               "darkgreen")) +  
  labs(
    title = "Healthcare Access Type by Divisions",
    x = "Division", 
    y = "Frequency",
    caption = "Source: Health Seek Data"
  ) +
  theme_minimal(base_size = 12) + 
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1, 
      vjust = 0.5, 
      face = "italic"),  
    legend.position = "top",  
    legend.title = element_blank(), 
    panel.grid.major = element_blank(), #element_line(color = "grey80",
                          #          linewidth = 0.2),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      lineheight = 12),  
    plot.caption = element_text(
      hjust = 1, 
      size = 8, 
      color = "grey50")  
  )

boot_fn_mlt <- function(y) {
  map_dfr(1:1000, function(x) {
    health_seek_data |> 
      select(HH7A, health_care) |> 
      filter(HH7A == {{y}}) |> 
      sample_n(size = n(), replace = TRUE) |> 
      count(health_care, HH7A) |> 
      mutate(p = n/sum(n))
  }) |> 
    reframe(
      prop = mean(p),
      sd = sd(p),
      .by = c(HH7A, health_care)
    )
}

x <- health_seek_data |> 
  count(y = as.character(HH7A)) |> pull(y)
x

dat_boot_mlt <- map_dfr(x, boot_fn_mlt)

health_seek_data |> 
  count(health_care) |> 
  reframe(original_p = n/sum(n),
          .by = c(HH7A)) |> 
  select(original_p)
fn_err_mlt <- function(x){
  dat_boot_mlt |> 
    bind_cols(health_seek_data |> 
                count(HH7A, health_care) |> 
                reframe(original_p = n/sum(n),
                        .by = c(HH7A)) |> 
                select(original_p)) |> 
    mutate(
      ymin = original_p - 2*sd,
      ymax = original_p + 2*sd
    ) |> 
    pivot_wider(
      names_from = health_care,
      values_from = original_p
    ) |> 
    select(HH7A, sd, ymin, ymax, {{x}}) |>
    drop_na() |> 
    arrange({{x}}) |>  
    mutate(HH7A = factor(HH7A, levels = HH7A)) |> 
    ggplot(
      mapping = aes(
        x = HH7A,
        y = {{x}}
      ))+
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    geom_point(size = 2,
               color = "black") +
    labs(
      title = paste0(rlang::as_name(rlang::enquo(x))),
      x = "Districts",
      y = "Proportion"
    ) +
    theme_minimal(base_size = 12) + 
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1, 
        vjust = 0.5, 
        face = "italic"),  
      axis.title.x = element_blank(),
      legend.position = "top",  
      legend.title = element_blank(), 
      panel.grid.major = element_blank(), 
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        lineheight = 12),  
      plot.caption = element_text(
        hjust = 1, 
        size = 8, 
        color = "grey50")  
    )
}

fn_err_mlt(govt)
fn_err(Stunted)
fn_err(`Severely Stunted`)
gridExtra::grid.arrange(fn_err_mlt(govt),
                        fn_err_mlt(private),
                        fn_err_mlt(Informal))

dat_boot_mlt |> 
  bind_cols(health_seek_data |> 
              count(HH7A, health_care) |> 
              reframe(original_p = n/sum(n),
                      .by = c(HH7A)) |> 
              select(original_p)) |> 
  mutate(
    ymin = original_p - 2*sd,
    ymax = original_p + 2*sd
  ) |> 
  arrange(HH7A, original_p) |>  
 # mutate(HH7A = factor(HH7A, levels = HH7A)) |> 
  ggplot(
    mapping = aes(
      x = HH7A,
      y = original_p,
      color = health_care
    ))+
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  geom_point(size = 2) +
  labs(
    #title = paste0(rlang::as_name(rlang::enquo(x))),
    x = "Districts",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 12) + 
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1, 
      vjust = 0.5, 
      face = "italic"),  
    axis.title.x = element_blank(),
    legend.position = "top",  
    legend.title = element_blank(), 
    panel.grid.major = element_blank(), 
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      lineheight = 12),  
    plot.caption = element_text(
      hjust = 1, 
      size = 8, 
      color = "grey50")  
  )
