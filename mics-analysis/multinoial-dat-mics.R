library(tidyverse)
library(haven)
library(foreign)

## reading the data

health_seek_data <- read_sav("micsdatasets/ch.sav") |>
  select(CA6A:CA6X, CA21A:CA21X, CAGE, HH7, HH7A, ethnicity,
         HL4, HH6, caretakerdis, windex5, melevel, BD3) |> 
  pivot_longer(cols = c(CA6A:CA6X, CA21A:CA21X),
               names_to = "health_seek",
               values_to = "health_value") |> 
  mutate(
    health_care = case_when(
      health_value %in% C(LETTERS[1:5], LETTERS[8]) ~ "govt",
      health_value %in% LETTERS[9:15] ~ "private", 
      health_value %in% LETTERS[16:18] ~ "Informal",
      health_value %in% c("W", "X") ~ NA_character_
    )) |>  drop_na() |> 
  mutate(
    caretakerdis = na_if(as.character(caretakerdis), "No information"),
    health_care = relevel(factor(health_care), ref = "Informal"),
    HH7 = as_factor(HH7),
    HH7A = as_factor(HH7A),
    HL4 = as_factor(HL4),
    windex5 = as_factor(windex5),
    melevel = na_if(as.character(melevel), "Missing/DK"),
    melevel = as_factor(melevel),
    HH6 = as_factor(HH6),
    caretakerdis = as_factor(caretakerdis),
    BD3 = na_if(as.character(BD3), "DK"),
    BD3 = na_if(as.character(BD3), "NO RESPONSE"),
    BD3 = as_factor(BD3),
    ethnicity = as_factor(ethnicity)
    ) |> 
  drop_na() |> 
  mutate(
    melevel = case_when(
      melevel == 0 ~ "None",
      melevel == 1 ~ "Primary",
      melevel == 2 ~ "Secondary",
      melevel == 3 ~ "Higher",
    .ptype = fct(levels = c("None", "Primary", "Secondary", "Higher"))),
    caretakerdis = case_when(
      caretakerdis == 1 ~ "Has Has functional difficulty",
      caretakerdis == 2 ~ "Has no functional difficulty",
      caretakerdis == 9 ~ NA
    ),
    BD3 = case_when(
      BD3 == 1 ~ "Breast Feeding",
      BD3 == 2 ~ "Not Breast Feeding"
    )
  ) |> drop_na()
  
#view(health_seek_data)
save(health_seek_data, file = "mics-analysis/health_seek_data.Rdata")
