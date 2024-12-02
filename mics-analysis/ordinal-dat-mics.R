library(tidyverse)
library(haven)
library(foreign)

ordinal_dat_mics <- read_sav("micsdatasets/ch.sav") |> 
  select(CAGE, HH7, HH7A, ethnicity, HAZ, WAZ, WHZ, FLAG, BMI, 
         HAZ2, WAZ2, WHZ2, ZBMI, HAZFLAG, WAZFLAG, WHZFLAG,
         BMIFLAG, HL4, HH6, caretakerdis, windex5, melevel) |> 
  mutate(
    # HAZ = if_else(FLAG == 1, NA, HAZ),
    # WAZ = if_else(FLAG == 1, NA, WAZ),
    # WHZ = if_else(FLAG == 1, NA, WHZ),
    HAZ2 = if_else(HAZFLAG == 1, NA, HAZ2),
    WAZ2 = if_else(WAZFLAG == 1, NA, WAZ2),
    WHZ2 = if_else(WHZFLAG == 1, NA, WHZ2),
    ZMBI = if_else(BMIFLAG == 1, NA, ZBMI),
    HH7 = as_factor(HH7),
    HH7A = as_factor(HH7A),
    HL4 = as_factor(HL4),
    windex5 = na_if(as.numeric(windex5), 0),
    melevel = na_if(as.character(melevel),
                    "Missing/DK"),
    melevel = case_when(
      melevel == 0 ~ "None",
      melevel == 1 ~ "Primary",
      melevel == 2 ~ "Secondary",
      melevel == 3 ~ "Higher",
      .ptype = fct(levels = c("None", "Primary", "Secondary", "Higher"))),
    HH6 = as_factor(HH6),
    caretakerdis = as_factor(caretakerdis),
    CAGE = case_when(
      CAGE < 12 ~ "<1 Year",
      CAGE >= 12 & CAGE <= 36 ~ "1-3 Years",
      CAGE > 36 & CAGE <= 60 ~ "3-5 Years",
      .ptype = fct(levels = c("<1 Year", "1-3 Years", 
                              "3-5 Years"))
    ),
    #ethnicity = as_factor(ethnicity),
    stunting = case_when(
      HAZ2 <= -3 ~ "Severely Stunted",
      HAZ2 > -3 & HAZ2 <= -2 ~ "Stunted",
      HAZ2 > -2 & HAZ2 <= 3 ~ "Normal",
      HAZ2 > 3 ~ NA,
      .ptype = fct(levels = c("Normal", "Stunted", "Severely Stunted"))
    ),
    wasting = case_when(
      WHZ2 <= -3 ~ "Severely Wasted",
      WHZ2 > -3 & WHZ2 <= -2 ~ "Wasted",
      WHZ2 > -2 & WHZ2 <= 3 ~ "Normal",
      WHZ2 > 3 ~ NA,
      .ptype = fct(levels = c("Normal", "Wasted", "Severely Wasted"))
    ),
    underweight = case_when(
      WAZ2 <= -3 ~ "Severely Underweight",
      WAZ2 > -3 & WAZ2 <= -2 ~ "Underweight",
      WAZ2 > -2 & WAZ2 <= 3 ~ "Normal",
      WAZ2 > 3 ~ NA,
      .ptype = fct(levels = c("Normal", "Underweight", "Severely Underweight"))
    )
  ) |> 
  drop_na()|>
  mutate(
    windex5 = case_when(
      windex5 == 1 ~ "Poorest",
      windex5 == 2 ~ "Second",
      windex5 == 3 ~ "Middle",
      windex5 == 4 ~ "Fourth",
      windex5 == 5 ~ "Richest",
      .ptype = fct(levels = c("Poorest", "Second", "Middle",
                              "Fourth", "Richest"))
    )
  )

# view(ordinal_dat)

save(ordinal_dat_mics, file = "mics-analysis/ordinal_dat_mics.Rdata")
