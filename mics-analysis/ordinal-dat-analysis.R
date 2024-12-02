library(tidyverse)
library(haven)
library(foreign)
library(sjPlot)
library(gt)
library(gtsummary)
library(broom.mixed)

load("mics-analysis/ordinal_dat_mics.Rdata")

ordinal_dat_mics |> 
  count(stunting)

library(ordinal)
library(brant)

ordinal_model_null <- clmm2(stunting ~ 1,
                           random = HH7A,
                           Hess = TRUE,
                           data = ordinal_dat_mics,
                           control = clmm2.control(
                             method = "nlminb",
                             maxIter = 100, 
                             gradTol = 1e-3),
                           threshold = "flexible",
                           link = "logistic")

pr_null <- profile(ordinal_model_null, nSteps = 10)
confint(pr_null)

ordinal_model <- clmm(
  stunting ~ CAGE + HH7 + HL4 + 
    HH6 + windex5 + melevel + BMI + (1|HH7A),
  data = ordinal_dat_mics,
  threshold = "equidistant",
  link = "logit"
)

ordinal_model_flexible <- clmm2(
  stunting ~ CAGE + HH7 + BMI + 
    HH6 + windex5 + melevel,
  random = HH7A,
  data = ordinal_dat_mics,
  Hess = TRUE,
  control = clmm2.control(method = "nlminb",
                          maxIter = 100, gradTol = 1e-4),
  threshold = "flexible",
  link = "logistic"
)

pr_flex <- profile(ordinal_model_flexible, nSteps = 10)
confint(pr_flex)
ordinal_sum_tab <- ordinal_dat_mics |> 
  select(CAGE, HH7, HL4, HH6, stunting,
         windex5, melevel, BMI) |> 
  tbl_summary(
    by = stunting,
    percent = "row",
    missing = "no",
    label = list(
      CAGE ~ "Child Age",
      melevel ~ "Mother's Education Level",
      windex5 ~ "Wealth Index",
      BMI ~ "Body Mass Index"
    ),
    digits = everything() ~ 2
  ) |> 
  add_p() |> 
  add_overall()

save(ordinal_sum_tab, file = "mics-analysis/ordinal_sum_tab.Rdata")
save(ordinal_model_sum, file = "mics-analysis/ordinal_model_sum.Rdata")
save(ordinal_model, file = "mics-analysis/ordinal_model.Rdata")
save(ordinal_model_null, file = "mics-analysis/ordinal_model_null.Rdata")

tbl_regression(ordinal_model)
ordinal_model_sum <- tidy(ordinal_model_flexible, 
                          exponentiate = TRUE,
                          conf.int = TRUE) |> 
  select(term, estimate, conf.low, conf.high, p.value)
ordinal_model_sum 

ordinal_model$ST


ordinal_dat_mics |> 
  select(stunting, HH7) |> 
  count(
    stunting, division = HH7
  ) |> 
  mutate(
    division = fct_reorder(division, n, .desc = TRUE)  
  ) |>
  ggplot(
    mapping = aes(
      x = division,
      y = n, 
      fill = stunting
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
    title = "Stunting by Division",
    x = "Divisions", 
    y = "Frequency",
    caption = "Source: MICS Data"
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

###################### bootstrap #################
x <- ordinal_dat_mics |> 
  count(y = as.character(HH7A)) |> pull(y)
x

boot_fn <- function(y) {
  map_dfr(1:1000, function(x) {ordinal_dat_mics |> 
      select(HH7A, stunting) |> 
      filter(HH7A == {{y}}) |> 
      sample_n(size = n(), replace = TRUE) |> 
      count(stunting, HH7A) |> 
      mutate(p = n/sum(n))
  }) |> 
    reframe(
      prop = mean(p),
      sd = sd(p),
      .by = c(HH7A, stunting)
    )
}

dat_boot <- map_dfr(x, boot_fn)
save(dat_boot, file = "../mics-analysis/dat_boot_ord.Rdata")
load("mics-analysis/dat_boot_ord.Rdata")
ordinal_dat_mics |> 
  count(HH7A, stunting) |> 
  reframe(original_p = n/sum(n),
         .by = c(HH7A))

fn_err <- function(x){
  dat_boot |> 
    bind_cols(ordinal_dat_mics |> 
                count(HH7A, stunting) |> 
                reframe(original_p = n/sum(n),
                        .by = c(HH7A)) |> 
                select(original_p)) |> 
  mutate(
    ymin = original_p - 2*sd,
    ymax = original_p + 2*sd
  ) |> 
  pivot_wider(
    names_from = stunting,
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

fn_err(Normal)
fn_err(Stunted)
fn_err(`Severely Stunted`)
gridExtra::grid.arrange(fn_err(Normal),
                        fn_err(Stunted),
                        fn_err(`Severely Stunted`))
