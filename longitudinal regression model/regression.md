treatment effect
================

## GitHub Documents

This file read in a prism exported text file and organize the table for
linear regression for the effect of the treatment (TGFB1 or MT-ND5
mutation effects).

It compares and selects the best fitted model.

### data 1: RL1, 121522-12a-soft-agar-d12-10000seeding

``` r
library(lme4)
```

    ## Loading required package: Matrix

``` r
library(lmerTest)
```

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
# Load the data
filename <- "/Users/yuanyuan/Downloads/github/treatment-effect-regression/data/12a_tgfb1_d12_soft_agar_both_count.txt"

df <- read.table(filename, header = TRUE, sep = "\t", quote = "\"", fill = TRUE)

library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ tidyr::expand() masks Matrix::expand()
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ tidyr::pack()   masks Matrix::pack()
    ## ✖ tidyr::unpack() masks Matrix::unpack()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
long_df <- df %>%
  pivot_longer(
    cols = everything(), # Selects all columns to transform
    names_to = "treatment", # This is the name of the new column for the variable names (corrected typo)
    values_to = "value" # This is the name of the new column for the values
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    tgfb1 = case_when(
      treatment == "c1"       ~ 0,
      treatment == "c2"       ~ 0,
      treatment == "c1.tgfb1" ~ 1,
      treatment == "c2.tgfb1" ~ 1,
      treatment == "n.tgfb1"  ~ 1,
      treatment == "n"        ~ 0,
      TRUE                    ~ NA_real_
    ),
    nd5 = case_when(
      treatment == "c1"       ~ 0,
      treatment == "c2"       ~ 0,
      treatment == "c1.tgfb1" ~ 0,
      treatment == "c2.tgfb1" ~ 0,
      treatment == "n.tgfb1"  ~ 1,
      treatment == "n"        ~ 1,
      TRUE                    ~ NA_real_
    ),
    bioreplicate = case_when(
      treatment == "c1"       ~ 1,
      treatment == "c2"       ~ 2,
      treatment == "c1.tgfb1" ~ 4,
      treatment == "c2.tgfb1" ~ 5,
      treatment == "n.tgfb1"  ~ 6,
      treatment == "n"        ~ 3,
      TRUE                    ~ NA_integer_
    ),
     nd5_3_level = case_when(
      treatment == "c1"       ~ 1,
      treatment == "c2"       ~ 2,
      treatment == "c1.tgfb1" ~ 1,
      treatment == "c2.tgfb1" ~ 2,
      treatment == "n.tgfb1"  ~ 3,
      treatment == "n"        ~ 3,
      TRUE                    ~ NA_integer_
    )
  )


model1 <- lm(value ~ tgfb1 * nd5, data = long_df)
summary(model1)
```

    ## 
    ## Call:
    ## lm(formula = value ~ tgfb1 * nd5, data = long_df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -81.467 -20.722   8.533  19.500  51.533 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  271.467      7.835  34.648   <2e-16 ***
    ## tgfb1       -164.133     10.608 -15.472   <2e-16 ***
    ## nd5           -3.356     12.794  -0.262    0.794    
    ## tgfb1:nd5     28.689     17.809   1.611    0.114    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 30.34 on 47 degrees of freedom
    ## Multiple R-squared:  0.8764, Adjusted R-squared:  0.8685 
    ## F-statistic: 111.1 on 3 and 47 DF,  p-value: < 2.2e-16

``` r
model2 <- lm(value ~ tgfb1 + nd5, data = long_df)
summary(model2)
```

    ## 
    ## Call:
    ## lm(formula = value ~ tgfb1 + nd5, data = long_df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -75.914 -23.961   7.634  21.860  57.086 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  265.914      7.152  37.181   <2e-16 ***
    ## tgfb1       -153.953      8.661 -17.775   <2e-16 ***
    ## nd5           11.452      9.046   1.266    0.212    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 30.84 on 48 degrees of freedom
    ## Multiple R-squared:  0.8696, Adjusted R-squared:  0.8642 
    ## F-statistic: 160.1 on 2 and 48 DF,  p-value: < 2.2e-16

``` r
anova(model2, model1)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: value ~ tgfb1 + nd5
    ## Model 2: value ~ tgfb1 * nd5
    ##   Res.Df   RSS Df Sum of Sq      F Pr(>F)
    ## 1     48 45666                           
    ## 2     47 43277  1    2389.5 2.5951 0.1139

``` r
#pick model 2 as model 1 is not significantly better

model3 <- lm(value ~ tgfb1 + nd5_3_level, data = long_df)
summary(model3)
```

    ## 
    ## Call:
    ## lm(formula = value ~ tgfb1 + nd5_3_level, data = long_df)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -66.506 -24.688   4.222  20.812  54.314 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  244.326     12.605  19.383   <2e-16 ***
    ## tgfb1       -152.908      8.358 -18.294   <2e-16 ***
    ## nd5_3_level   12.180      5.200   2.342   0.0234 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 29.7 on 48 degrees of freedom
    ## Multiple R-squared:  0.8791, Adjusted R-squared:  0.874 
    ## F-statistic: 174.5 on 2 and 48 DF,  p-value: < 2.2e-16

``` r
#####this result of model 3 is confusing##########
#to resolve this, c1 data is deleted#############

df_long2 <- long_df %>% filter(nd5_3_level!=1) #eliminate c1s
model4 <- lm(value ~ tgfb1 + nd5, data = df_long2)
summary(model4)
```

    ## 
    ## Call:
    ## lm(formula = value ~ tgfb1 + nd5, data = df_long2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -54.944 -22.528   1.389  17.278  44.389 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  278.6111     7.6227  36.550   <2e-16 ***
    ## tgfb1       -156.6667     8.8020 -17.799   <2e-16 ***
    ## nd5            0.1111     8.8020   0.013     0.99    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 26.41 on 33 degrees of freedom
    ## Multiple R-squared:  0.9057, Adjusted R-squared:  0.8999 
    ## F-statistic: 158.4 on 2 and 33 DF,  p-value: < 2.2e-16
