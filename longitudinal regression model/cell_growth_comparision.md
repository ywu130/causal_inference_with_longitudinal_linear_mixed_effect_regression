cell growth test
================
Yuanyuan Wu
2023-07-15

## Description

This analysis examines the longitudinal growth rate of cells and
identifies treatments that significantly affect cellular growth over
time. A semi-logarithmic regression model with interaction terms was
employed to capture the combined effects of time (Hours), the presence
of MT-ND5(nd5) mutations, and tgfb1 treatment. Interaction terms with
time were included to assess how these factors jointly influence growth
trajectories. Treatments with significant interaction terms with time
were considered to have a measurable impact on the cellular growth rate.

The script below contains data organization and transformation, fitting
of the model, model comparison, and model visualization on multiple
datasets

Cell growth data measured by Celigo, from 090622,293T, r8,d25

### Fitting with treatment time interation

$$\log_{10}(\text{cell count}) = \beta_0 + \beta_1 (\text{treatment}) + \beta_2 (\text{time}) + \beta_3 (\text{treatment} \times \text{time}) + \epsilon$$

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * time..hour., data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36941 -0.17756 -0.01379  0.13046  0.46721 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               6.257378   0.083814  74.658   <2e-16 ***
    ## treatmentnd5              0.287164   0.124316   2.310   0.0250 *  
    ## time..hour.               0.045000   0.001082  41.582   <2e-16 ***
    ## treatmentnd5:time..hour. -0.004016   0.001605  -2.502   0.0156 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1973 on 51 degrees of freedom
    ## Multiple R-squared:  0.9829, Adjusted R-squared:  0.9819 
    ## F-statistic: 974.8 on 3 and 51 DF,  p-value: < 2.2e-16

``` r
library(ggplot2)
# Create a new data frame for predictions
newdata <- expand.grid(time..hour. = seq(min(df_long$time..hour.), max(df_long$time..hour.), length.out = 100),
                       treatment = c("ctr", "nd5"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = time..hour., y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = time..hour., y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/plot-1.png)<!-- -->

``` r
# Extract coefficients and confidence intervals
coef <- summary(model)$coefficients
slopes <- c(coef["time..hour.", "Estimate"], coef["time..hour.", "Estimate"] + coef["treatmentnd5:time..hour.", "Estimate"])
CI <- confint(model)[c("time..hour.", "treatmentnd5:time..hour."), ]


# Extract coefficients, standard errors, and covariance
coef <- summary(model)$coefficients
cov <- vcov(model)

# Compute slopes and their standard errors
slope.ctr <- coef["time..hour.", "Estimate"]
slope.nd5 <- slope.ctr + coef["treatmentnd5:time..hour.", "Estimate"]
se.ctr <- coef["time..hour.", "Std. Error"]
se.nd5 <- sqrt(se.ctr^2 + coef["treatmentnd5:time..hour.", "Std. Error"]^2 + 2 * cov["time..hour.", "treatmentnd5:time..hour."])

# Compute 95% confidence intervals
CI.ctr <- c(slope.ctr - 1.96 * se.ctr, slope.ctr + 1.96 * se.ctr)
CI.nd5 <- c(slope.nd5 - 1.96 * se.nd5, slope.nd5 + 1.96 * se.nd5)

# Create a data frame for slopes
slope_data <- data.frame(treatment = c("ctr", "nd5"), slope = c(slope.ctr, slope.nd5), lower = c(CI.ctr[1], CI.nd5[1]), upper = c(CI.ctr[2], CI.nd5[2]))


# # Plot the observed data
# ggplot(df_long, aes(x = time..hour., y = log(cell_count), color = treatment)) +
#   geom_point() +
#   # Add the model predictions
#   geom_line(data = newdata, aes(x = time..hour., y = pred)) +
#   # Add the slopes
#   geom_text(data = slope_data, aes(x = Inf, y = slope, label = paste("Slope =", round(slope, 2))), hjust = "inward", vjust = "inward") +
#   labs(x = "Time (hours)", y = "Log-transformed cell count") +
#   theme_minimal()
# 
# # Plot the observed data
# ggplot(df_long, aes(x = time..hour., y = log(cell_count), color = treatment)) +
#   geom_point() +
#   # Add the model predictions
#   geom_line(data = newdata, aes(x = time..hour., y = pred)) +
#   # Add the slopes
#   geom_text(data = slope_data, aes(x = Inf, y = slope, label = paste("Slope =", round(slope, 2))), hjust = "inward", vjust = "inward") +
#   # Add error bars for the confidence intervals
#   geom_errorbarh(data = slope_data, aes(xmin = lower, xmax = upper, y = slope), height = 0.1, color = "black") +
#   labs(x = "Time (hours)", y = "Log-transformed cell count") +
#   theme_minimal()

# 
# # Plot the observed data
# ggplot(df_long, aes(x = time..hour., y = log(cell_count), color = treatment)) +
#   geom_point() +
#   # Add the model predictions
#   geom_line(data = newdata, aes(x = time..hour., y = pred)) +
#   # Add the slopes
#   geom_text(data = slope_data, aes(x = Inf, y = slope, label = paste("Slope =", round(slope, 2))), hjust = "inward", vjust = "inward") +
#   # Add error bars for the confidence intervals
#   geom_errorbarh(data = slope_data, aes(x = slope, xmin = lower, xmax = upper, y = treatment), height = 0.1, color = "black") +
#   labs(x = "Time (hours)", y = "Log-transformed cell count") +
#   theme_minimal()


# Convert treatment to a factor with numeric levels
df_long$treatment_num <- as.numeric(factor(df_long$treatment))
slope_data$treatment_num <- as.numeric(factor(slope_data$treatment))

# Plot the observed data
ggplot(df_long, aes(x = time..hour., y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = time..hour., y = pred)) +
  # Add the slopes
  geom_text(data = slope_data, aes(x = Inf, y = slope, label = paste("Slope =", round(slope, 2))), hjust = "inward", vjust = "inward") +
  # Add error bars for the confidence intervals
  geom_errorbarh(data = slope_data, aes(x = slope, xmin = lower, xmax = upper, y = treatment_num), height = 0.1, color = "black") +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_minimal()
```

    ## Warning in geom_errorbarh(data = slope_data, aes(x = slope, xmin = lower, :
    ## Ignoring unknown aesthetics: x

![](cell_growth_comparision_files/figure-gfm/plot-2.png)<!-- --> \###
Slope with se when fitting seperately

    ##          slope           se
    ## ctr 0.04499957 0.0009670167
    ## nd5 0.04098406 0.0013229113

### test on 090722, hek293_R8_d26

### Fitting hek293 101022_R3_d12

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("101022_R3_d12.txt", header=TRUE)

# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
df_long$treatment <- ifelse(startsWith(df_long$replicate, "ctr"), "ctr", "nd5")


# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary


model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.88918 -0.17946  0.03507  0.23956  0.67474 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         7.291698   0.108984  66.906   <2e-16 ***
    ## treatmentnd5       -0.095796   0.154127  -0.622    0.537    
    ## Hours               0.044332   0.001957  22.656   <2e-16 ***
    ## treatmentnd5:Hours -0.001046   0.002767  -0.378    0.707    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3523 on 56 degrees of freedom
    ## Multiple R-squared:  0.9472, Adjusted R-squared:  0.9444 
    ## F-statistic: 335.1 on 3 and 56 DF,  p-value: < 2.2e-16

``` r
# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("ctr", "nd5"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/101022_R3_d12%20data-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##          slope          se
    ## ctr 0.04433238 0.001987630
    ## nd5 0.04328651 0.001925311

### Fititng 101622-hek293-d18-R3

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("101622-hek293-d18-R3 for export.txt", header=TRUE)

# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
df_long$treatment <- ifelse(startsWith(df_long$replicate, "ctr"), "ctr", "nd5")


# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary


model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.37432 -0.13875  0.02593  0.13240  0.31755 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         7.165941   0.055423 129.296   <2e-16 ***
    ## treatmentnd5       -0.130441   0.078379  -1.664    0.102    
    ## Hours               0.050124   0.001273  39.361   <2e-16 ***
    ## treatmentnd5:Hours -0.002363   0.001801  -1.312    0.195    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1651 on 56 degrees of freedom
    ## Multiple R-squared:  0.9816, Adjusted R-squared:  0.9806 
    ## F-statistic: 993.9 on 3 and 56 DF,  p-value: < 2.2e-16

``` r
# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("ctr", "nd5"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/101622-hek293-d18-R3%20data-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##          slope          se
    ## ctr 0.05012363 0.001208241
    ## nd5 0.04776030 0.001335442

### Fitting 011723_hek293_rl2

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("011723_hek293_rl2_export.txt", header=TRUE)
```

    ## Warning in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,
    ## : number of items read is not a multiple of the number of columns

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")




# Create a treatment variable
df_long$treatment <- ifelse(startsWith(df_long$replicate, "c2"), "c2", 
                            ifelse(startsWith(df_long$replicate, "c1"), "c1", "n"))


#reorder to n, c1, c2, so that n is the reference level
df_long$treatment <- factor(df_long$treatment, levels = c("n", "c1", "c2"))


# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary


model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.22778 -0.09781  0.03166  0.06820  0.15469 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        7.8344326  0.0481319 162.770   <2e-16 ***
    ## treatmentc1       -0.0766264  0.0674945  -1.135    0.262    
    ## treatmentc2        0.0361337  0.0674945   0.535    0.595    
    ## Hours              0.0531816  0.0012313  43.192   <2e-16 ***
    ## treatmentc1:Hours -0.0003765  0.0016848  -0.223    0.824    
    ## treatmentc2:Hours -0.0007471  0.0016848  -0.443    0.659    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1066 on 47 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.9923, Adjusted R-squared:  0.9915 
    ## F-statistic:  1213 on 5 and 47 DF,  p-value: < 2.2e-16

``` r
# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("c1", "c2","n"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](cell_growth_comparision_files/figure-gfm/011723_hek293_rl2-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##         slope          se
    ## n  0.05318162 0.001249880
    ## c1 0.05280510 0.001249627
    ## c2 0.05243457 0.001022465

### Fitting 120921_12a_d11

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("120921_12a_d11.txt", header=TRUE)
```

    ## Warning in read.table("120921_12a_d11.txt", header = TRUE): incomplete final
    ## line found by readTableHeader on '120921_12a_d11.txt'

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")




# Create a treatment variable
df_long$treatment <- ifelse(startsWith(df_long$replicate, "c1"), "c1", "n")
                            #ifelse(startsWith(df_long$replicate, "c2-"), "c2-", "n"))





# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary


model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.190640 -0.063229  0.005006  0.045407  0.166651 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       6.4302389  0.0506520 126.949   <2e-16 ***
    ## treatmentn       -0.0183811  0.0716328  -0.257    0.799    
    ## Hours             0.0346659  0.0005964  58.128   <2e-16 ***
    ## treatmentn:Hours -0.0000419  0.0008434  -0.050    0.961    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08536 on 44 degrees of freedom
    ## Multiple R-squared:  0.9935, Adjusted R-squared:  0.9931 
    ## F-statistic:  2250 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("c1","n"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/120921_12a_d11-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##         slope           se
    ## c1 0.03466589 0.0005152990
    ## n  0.03462399 0.0006676694

### Fitting 111521 hek293 d30

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("111521_hek293_d30_export.txt", header=TRUE)
```

    ## Warning in read.table("111521_hek293_d30_export.txt", header = TRUE):
    ## incomplete final line found by readTableHeader on
    ## '111521_hek293_d30_export.txt'

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
df_long$treatment <- ifelse(startsWith(df_long$replicate, "ctr"), "ctr", "nd5")


# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary


model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.53919 -0.19402 -0.06024  0.17391  0.60052 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         5.2287367  0.2996123  17.452  < 2e-16 ***
    ## treatmentnd5        0.1573843  0.4237158   0.371    0.713    
    ## Hours               0.0626272  0.0042715  14.662 9.43e-16 ***
    ## treatmentnd5:Hours -0.0002077  0.0060408  -0.034    0.973    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3118 on 32 degrees of freedom
    ## Multiple R-squared:  0.9308, Adjusted R-squared:  0.9243 
    ## F-statistic: 143.5 on 3 and 32 DF,  p-value: < 2.2e-16

``` r
# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("ctr", "nd5"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/111521_hek293_d30-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##          slope          se
    ## ctr 0.06262717 0.004217258
    ## nd5 0.06241949 0.004324998

### Fitting 032322_12A_R4_d10

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("032322_12a_R4_d10.txt", header=TRUE)

# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 2)  # get the first two characters
  if (prefix %in% c('a1', 'a2', 'a3')) {
    return('ctr-')
  } else if (prefix %in% c('b1', 'b2', 'b3')) {
    return('ctr+')
  } else if (prefix %in% c('a4', 'a5', 'a6')) {
    return('nd5-')
  } else if (prefix %in% c('b4', 'b5', 'b6')) {
    return('nd5+')
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)


#export the long table
#write.csv(df_long, file = "032322_12a_R4_d10_organized.csv")



# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54145 -0.12989 -0.01651  0.13085  0.57675 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          6.9374397  0.0724589  95.743  < 2e-16 ***
    ## treatmentctr+       -0.2956057  0.1024723  -2.885  0.00442 ** 
    ## treatmentnd5-        0.0737007  0.1024723   0.719  0.47298    
    ## treatmentnd5+       -0.3182300  0.1024723  -3.106  0.00222 ** 
    ## Hours                0.0308182  0.0009420  32.717  < 2e-16 ***
    ## treatmentctr+:Hours -0.0062117  0.0013321  -4.663 6.23e-06 ***
    ## treatmentnd5-:Hours -0.0003959  0.0013321  -0.297  0.76668    
    ## treatmentnd5+:Hours -0.0039265  0.0013321  -2.948  0.00365 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2087 on 172 degrees of freedom
    ## Multiple R-squared:  0.9598, Adjusted R-squared:  0.9581 
    ## F-statistic:   586 on 7 and 172 DF,  p-value: < 2.2e-16

``` r
# ###testing emmeans
# library(emmeans)
# emmeans_result <- emmeans(model, pairwise ~ treatment | Hours)
# summary(emmeans_result)

#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["nd5+"]])
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54145 -0.12989 -0.01651  0.13085  0.57675 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          6.619210   0.072459  91.351  < 2e-16 ***
    ## treatmentnd5-        0.391931   0.102472   3.825 0.000183 ***
    ## treatmentctr+        0.022624   0.102472   0.221 0.825523    
    ## treatmentctr-        0.318230   0.102472   3.106 0.002222 ** 
    ## Hours                0.026892   0.000942  28.548  < 2e-16 ***
    ## treatmentnd5-:Hours  0.003531   0.001332   2.650 0.008792 ** 
    ## treatmentctr+:Hours -0.002285   0.001332  -1.715 0.088074 .  
    ## treatmentctr-:Hours  0.003927   0.001332   2.948 0.003648 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2087 on 172 degrees of freedom
    ## Multiple R-squared:  0.9598, Adjusted R-squared:  0.9581 
    ## F-statistic:   586 on 7 and 172 DF,  p-value: < 2.2e-16

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("ctr-", "ctr+","nd5-","nd5+"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/032322_12A_R4_d10-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##           slope           se
    ## nd5+ 0.02689169 0.0007589819
    ## nd5- 0.03042230 0.0006494532
    ## ctr+ 0.02460654 0.0013980311
    ## ctr- 0.03081820 0.0007726009

``` r
##################fitting with interaction term of the treatment##############################################################################
# Splitting the treatment column into two separate columns for nd5 and tgfb1
df_long$nd5 <- ifelse(grepl("nd5", df_long$treatment), "nd5", "ctr")
df_long$tgfb1 <- ifelse(grepl("\\+", df_long$treatment), "plus", "minus")

# Fit the linear model
model <- lm(log(cell_count) ~ Hours * nd5 * tgfb1, data=df_long)

# Get the summary of the model
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 * tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54145 -0.12989 -0.01651  0.13085  0.57675 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             6.9374397  0.0724589  95.743  < 2e-16 ***
    ## Hours                   0.0308182  0.0009420  32.717  < 2e-16 ***
    ## nd5nd5                  0.0737007  0.1024723   0.719  0.47298    
    ## tgfb1plus              -0.2956057  0.1024723  -2.885  0.00442 ** 
    ## Hours:nd5nd5           -0.0003959  0.0013321  -0.297  0.76668    
    ## Hours:tgfb1plus        -0.0062117  0.0013321  -4.663 6.23e-06 ***
    ## nd5nd5:tgfb1plus       -0.0963249  0.1449177  -0.665  0.50714    
    ## Hours:nd5nd5:tgfb1plus  0.0026811  0.0018839   1.423  0.15652    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2087 on 172 degrees of freedom
    ## Multiple R-squared:  0.9598, Adjusted R-squared:  0.9581 
    ## F-statistic:   586 on 7 and 172 DF,  p-value: < 2.2e-16

``` r
model2 <- lm(log(cell_count) ~ Hours*nd5 + Hours* tgfb1, data=df_long)
summary(model2)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 + Hours * tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.57793 -0.14128 -0.00704  0.12442  0.53884 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      6.9615209  0.0631311 110.271  < 2e-16 ***
    ## Hours            0.0301479  0.0008207  36.734  < 2e-16 ***
    ## nd5nd5           0.0255382  0.0728975   0.350    0.727    
    ## tgfb1plus       -0.3437682  0.0728975  -4.716 4.92e-06 ***
    ## Hours:nd5nd5     0.0009446  0.0009477   0.997    0.320    
    ## Hours:tgfb1plus -0.0048711  0.0009477  -5.140 7.33e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.21 on 174 degrees of freedom
    ## Multiple R-squared:  0.9588, Adjusted R-squared:  0.9576 
    ## F-statistic: 809.8 on 5 and 174 DF,  p-value: < 2.2e-16

``` r
model3 <- lm(log(cell_count) ~ Hours*nd5 + Hours* tgfb1 + Hours:nd5:tgfb1, data=df_long)
summary(model3)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + 
    ##     Hours:nd5:tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.53278 -0.12834 -0.02721  0.12354  0.57468 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             6.9615209  0.0626499 111.118  < 2e-16 ***
    ## Hours                   0.0305355  0.0008391  36.390  < 2e-16 ***
    ## nd5nd5                  0.0255382  0.0723419   0.353   0.7245    
    ## tgfb1plus              -0.3437682  0.0723419  -4.752 4.22e-06 ***
    ## Hours:nd5nd5            0.0001695  0.0010235   0.166   0.8686    
    ## Hours:tgfb1plus        -0.0056462  0.0010235  -5.517 1.24e-07 ***
    ## Hours:nd5nd5:tgfb1plus  0.0015502  0.0008078   1.919   0.0566 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2084 on 173 degrees of freedom
    ## Multiple R-squared:  0.9597, Adjusted R-squared:  0.9583 
    ## F-statistic: 685.8 on 6 and 173 DF,  p-value: < 2.2e-16

``` r
anova(model3, model)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + Hours:nd5:tgfb1
    ## Model 2: log(cell_count) ~ Hours * nd5 * tgfb1
    ##   Res.Df    RSS Df Sum of Sq      F Pr(>F)
    ## 1    173 7.5139                           
    ## 2    172 7.4946  1  0.019251 0.4418 0.5071

``` r
anova(model2,model3)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(cell_count) ~ Hours * nd5 + Hours * tgfb1
    ## Model 2: log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + Hours:nd5:tgfb1
    ##   Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
    ## 1    174 7.6738                              
    ## 2    173 7.5139  1   0.15996 3.6831 0.05661 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pick model 3 as the complex model is not better
```

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("032322_12a_R4_d10.txt", header=TRUE)

# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 2)  # get the first two characters
  if (prefix %in% c('a1', 'a2', 'a3')) {
    return('Control')
  } else if (prefix %in% c('b1', 'b2', 'b3')) {
    return('Control_TGFB1+')
  } else if (prefix %in% c('a4', 'a5', 'a6')) {
    return('MT-ND5')
  } else if (prefix %in% c('b4', 'b5', 'b6')) {
    return('MT-ND5_TGFB1+')
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)


#export the long table
#write.csv(df_long, file = "032322_12a_R4_d10_organized.csv")



# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54145 -0.12989 -0.01651  0.13085  0.57675 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    6.9374397  0.0724589  95.743  < 2e-16 ***
    ## treatmentControl_TGFB1+       -0.2956057  0.1024723  -2.885  0.00442 ** 
    ## treatmentMT-ND5                0.0737007  0.1024723   0.719  0.47298    
    ## treatmentMT-ND5_TGFB1+        -0.3182300  0.1024723  -3.106  0.00222 ** 
    ## Hours                          0.0308182  0.0009420  32.717  < 2e-16 ***
    ## treatmentControl_TGFB1+:Hours -0.0062117  0.0013321  -4.663 6.23e-06 ***
    ## treatmentMT-ND5:Hours         -0.0003959  0.0013321  -0.297  0.76668    
    ## treatmentMT-ND5_TGFB1+:Hours  -0.0039265  0.0013321  -2.948  0.00365 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2087 on 172 degrees of freedom
    ## Multiple R-squared:  0.9598, Adjusted R-squared:  0.9581 
    ## F-statistic:   586 on 7 and 172 DF,  p-value: < 2.2e-16

``` r
# ###testing emmeans
# library(emmeans)
# emmeans_result <- emmeans(model, pairwise ~ treatment | Hours)
# summary(emmeans_result)

#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["nd5+"]])
```

    ## NULL

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("Control", "Control_TGFB1+","MT-ND5","MT-ND5_TGFB1+"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
my_plot <- ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment,shape = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()

# Save the plot as an SVG file
#ggsave(filename = "12A_growth_curve.svg", plot = my_plot, device = "svg")
ggsave(filename = "12A_growth_curve.svg", plot = my_plot, device = "svg", width = 8, height = 4, units = "in")
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##                     slope           se
    ## MT-ND5_TGFB1+  0.02689169 0.0007589819
    ## MT-ND5         0.03042230 0.0006494532
    ## Control_TGFB1+ 0.02460654 0.0013980311
    ## Control        0.03081820 0.0007726009

``` r
##################fitting with interaction term of the treatment##############################################################################
# Splitting the treatment column into two separate columns for nd5 and tgfb1
df_long$nd5 <- ifelse(grepl("ND5", df_long$treatment), "nd5", "ctr")
df_long$tgfb1 <- ifelse(grepl("\\+", df_long$treatment), "plus", "minus")

# Fit the linear model
model <- lm(log(cell_count) ~ Hours * nd5 * tgfb1, data=df_long)

# Get the summary of the model
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 * tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54145 -0.12989 -0.01651  0.13085  0.57675 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             6.9374397  0.0724589  95.743  < 2e-16 ***
    ## Hours                   0.0308182  0.0009420  32.717  < 2e-16 ***
    ## nd5nd5                  0.0737007  0.1024723   0.719  0.47298    
    ## tgfb1plus              -0.2956057  0.1024723  -2.885  0.00442 ** 
    ## Hours:nd5nd5           -0.0003959  0.0013321  -0.297  0.76668    
    ## Hours:tgfb1plus        -0.0062117  0.0013321  -4.663 6.23e-06 ***
    ## nd5nd5:tgfb1plus       -0.0963249  0.1449177  -0.665  0.50714    
    ## Hours:nd5nd5:tgfb1plus  0.0026811  0.0018839   1.423  0.15652    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2087 on 172 degrees of freedom
    ## Multiple R-squared:  0.9598, Adjusted R-squared:  0.9581 
    ## F-statistic:   586 on 7 and 172 DF,  p-value: < 2.2e-16

``` r
model2 <- lm(log(cell_count) ~ Hours*nd5 + Hours* tgfb1, data=df_long)
summary(model2)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 + Hours * tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.57793 -0.14128 -0.00704  0.12442  0.53884 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      6.9615209  0.0631311 110.271  < 2e-16 ***
    ## Hours            0.0301479  0.0008207  36.734  < 2e-16 ***
    ## nd5nd5           0.0255382  0.0728975   0.350    0.727    
    ## tgfb1plus       -0.3437682  0.0728975  -4.716 4.92e-06 ***
    ## Hours:nd5nd5     0.0009446  0.0009477   0.997    0.320    
    ## Hours:tgfb1plus -0.0048711  0.0009477  -5.140 7.33e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.21 on 174 degrees of freedom
    ## Multiple R-squared:  0.9588, Adjusted R-squared:  0.9576 
    ## F-statistic: 809.8 on 5 and 174 DF,  p-value: < 2.2e-16

``` r
model3 <- lm(log(cell_count) ~ Hours*nd5 + Hours* tgfb1 + Hours:nd5:tgfb1, data=df_long)
summary(model3)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + 
    ##     Hours:nd5:tgfb1, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.53278 -0.12834 -0.02721  0.12354  0.57468 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             6.9615209  0.0626499 111.118  < 2e-16 ***
    ## Hours                   0.0305355  0.0008391  36.390  < 2e-16 ***
    ## nd5nd5                  0.0255382  0.0723419   0.353   0.7245    
    ## tgfb1plus              -0.3437682  0.0723419  -4.752 4.22e-06 ***
    ## Hours:nd5nd5            0.0001695  0.0010235   0.166   0.8686    
    ## Hours:tgfb1plus        -0.0056462  0.0010235  -5.517 1.24e-07 ***
    ## Hours:nd5nd5:tgfb1plus  0.0015502  0.0008078   1.919   0.0566 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2084 on 173 degrees of freedom
    ## Multiple R-squared:  0.9597, Adjusted R-squared:  0.9583 
    ## F-statistic: 685.8 on 6 and 173 DF,  p-value: < 2.2e-16

``` r
anova(model3, model)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + Hours:nd5:tgfb1
    ## Model 2: log(cell_count) ~ Hours * nd5 * tgfb1
    ##   Res.Df    RSS Df Sum of Sq      F Pr(>F)
    ## 1    173 7.5139                           
    ## 2    172 7.4946  1  0.019251 0.4418 0.5071

``` r
anova(model2,model3)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: log(cell_count) ~ Hours * nd5 + Hours * tgfb1
    ## Model 2: log(cell_count) ~ Hours * nd5 + Hours * tgfb1 + Hours:nd5:tgfb1
    ##   Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
    ## 1    174 7.6738                              
    ## 2    173 7.5139  1   0.15996 3.6831 0.05661 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#pick model 3 as the complex model is not better
```

### Fitting 031123_12a_d10

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
#df <- read.table("031123_12a_d10.txt", header=TRUE)

#There is NA in the row
df <- read.delim("031123_12a_d10.txt", header = TRUE, na.strings = "", fill = TRUE)
head(df)
```

    ##       Hours  c10 c10.1 c10.2  c20 c20.1 c20.2  nd0 nd0.1 nd0.2  c11 c11.1 c11.2
    ## 1  0.000000  752   855   785  698   722   861  752   770   718 1292  1130  1139
    ## 2  3.364722 1194   965   982  994   970   932  924  1019   835 1293  1537  1334
    ## 3 24.947222 1650  1689  1536 1585  1674  1634 1404  1530  1455 2099  2508  2431
    ## 4 53.379167 2824  3393  2809   NA    NA    NA 2775  3246  2555 4468  5641  5063
    ## 5 74.638056 5947  5981  5280 4229  4432  4630 3671  4161  4965 8804  7818  7875
    ##    c21 c21.1 c21.2  nd1 nd1.1 nd1.2
    ## 1 1444  1407  1298 1046  1090   934
    ## 2 1537  1604  1489 1392  1169  1178
    ## 3 2692  2555  2604 2029  2309  1975
    ## 4 4489  5447  5544 4071  4950  4022
    ## 5 9925  9472 10178 6723  6716  7872

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")




# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 3)  # get the first two characters
  if (prefix == 'c10') {
    return('c1-')
  } else if (prefix == 'c11') {
    return('c1+')
  } else if (prefix == 'nd0') {
    return('nd5-')
  } else if (prefix == 'nd1') {
    return('nd5+')
  }  else if (prefix == 'c20') {
    return('c2-')
  } else if (prefix == 'c21') {
    return('c2+')
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)

##split the treatment into 2 seperate columns
# Load the dplyr package
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
# Assuming your data frame is named df
# Update the name accordingly if it differs
df_long <- df_long %>%
  mutate(
    treatment = as.character(treatment),
    # Create the 'tgfb1' column based on whether 'treatment' ends with '+' or '-'
    tgfb1 = case_when(
      grepl("\\+$", treatment) ~ "plus",
      grepl("\\-$", treatment) ~ "minus",
      TRUE ~ NA_character_
    ),

    # Create the 'nd5' column based on whether 'treatment' starts with 'c'
    nd5 = case_when(
      startsWith(treatment, "c") ~ "ctr",
      TRUE ~ "nd5"
    )
  )

# View the updated dataframe
#head(df_long)

#change the variable columns back to factor format
# Convert 'treatment', 'nd5', and 'tgfb1' to factors
df_long <- df_long %>%
  mutate(
    treatment = as.factor(treatment),
    nd5 = as.factor(nd5),
    tgfb1 = as.factor(tgfb1)
  )

# Verify the changes
#str(df_long)







# 
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.195632 -0.060645  0.000602  0.054007  0.238155 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          6.7638221  0.0360046 187.860  < 2e-16 ***
    ## treatmentc1+         0.3525777  0.0509183   6.924 1.29e-09 ***
    ## treatmentc2-        -0.0219494  0.0509680  -0.431   0.6680    
    ## treatmentc2+         0.4684363  0.0509183   9.200 6.21e-14 ***
    ## treatmentnd5-       -0.0675585  0.0509183  -1.327   0.1886    
    ## treatmentnd5+        0.2295088  0.0509183   4.507 2.38e-05 ***
    ## Hours                0.0246935  0.0008461  29.185  < 2e-16 ***
    ## treatmentc1+:Hours   0.0009946  0.0011966   0.831   0.4085    
    ## treatmentc2-:Hours  -0.0020969  0.0012470  -1.682   0.0968 .  
    ## treatmentc2+:Hours   0.0010236  0.0011966   0.855   0.3950    
    ## treatmentnd5-:Hours -0.0019669  0.0011966  -1.644   0.1044    
    ## treatmentnd5+:Hours  0.0007299  0.0011966   0.610   0.5437    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09459 on 75 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.9868, Adjusted R-squared:  0.9848 
    ## F-statistic: 508.1 on 11 and 75 DF,  p-value: < 2.2e-16

``` r
# ###testing emmeans
# library(emmeans)
# emmeans_result <- emmeans(model, pairwise ~ treatment | Hours)
# summary(emmeans_result)


#######################effect of tgfb1 and nd5###################


#model <- lm(log(cell_count) ~ tgfb1*nd5 * Hours, data = df_long)
#summary(model)
# Fit the linear model with specified interaction terms
model <- lm(log(cell_count) ~ tgfb1 + nd5 + Hours + tgfb1:Hours + nd5:Hours, data = df_long)

# Print the summary of the model
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ tgfb1 + nd5 + Hours + tgfb1:Hours + 
    ##     nd5:Hours, data = df_long)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.226704 -0.077490  0.000982  0.078115  0.230237 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      6.7749234  0.0257888 262.708  < 2e-16 ***
    ## tgfb1plus        0.3789380  0.0326080  11.621  < 2e-16 ***
    ## nd5nd5          -0.1195951  0.0345816  -3.458 0.000869 ***
    ## Hours            0.0237475  0.0006225  38.150  < 2e-16 ***
    ## tgfb1plus:Hours  0.0021034  0.0007747   2.715 0.008092 ** 
    ## nd5nd5:Hours    -0.0007242  0.0008171  -0.886 0.378098    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1049 on 81 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.9824, Adjusted R-squared:  0.9813 
    ## F-statistic: 905.1 on 5 and 81 DF,  p-value: < 2.2e-16

``` r
#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["c1-"]])
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.195632 -0.060645  0.000602  0.054007  0.238155 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          6.7638221  0.0360046 187.860  < 2e-16 ***
    ## treatmentc1+         0.3525777  0.0509183   6.924 1.29e-09 ***
    ## treatmentc2-        -0.0219494  0.0509680  -0.431   0.6680    
    ## treatmentc2+         0.4684363  0.0509183   9.200 6.21e-14 ***
    ## treatmentnd5-       -0.0675585  0.0509183  -1.327   0.1886    
    ## treatmentnd5+        0.2295088  0.0509183   4.507 2.38e-05 ***
    ## Hours                0.0246935  0.0008461  29.185  < 2e-16 ***
    ## treatmentc1+:Hours   0.0009946  0.0011966   0.831   0.4085    
    ## treatmentc2-:Hours  -0.0020969  0.0012470  -1.682   0.0968 .  
    ## treatmentc2+:Hours   0.0010236  0.0011966   0.855   0.3950    
    ## treatmentnd5-:Hours -0.0019669  0.0011966  -1.644   0.1044    
    ## treatmentnd5+:Hours  0.0007299  0.0011966   0.610   0.5437    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09459 on 75 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.9868, Adjusted R-squared:  0.9848 
    ## F-statistic: 508.1 on 11 and 75 DF,  p-value: < 2.2e-16

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = c("c1-", "c1+","nd5-","nd5+","c2-","c2+"))

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](cell_growth_comparision_files/figure-gfm/031123_12a_d10-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##           slope           se
    ## nd5+ 0.02542342 0.0008367039
    ## nd5- 0.02272666 0.0009388385
    ## c2+  0.02571713 0.0006020527
    ## c2-  0.02259663 0.0010063236
    ## c1+  0.02568813 0.0007685567
    ## c1-  0.02469354 0.0009648141

### Fitting 052721_12A_sup

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
#df <- read.table("031123_12a_d10.txt", header=TRUE)

#There is NA in the row
df <- read.delim("052721_12A_sup.txt", header = TRUE, na.strings = "", fill = TRUE)
head(df)
```

    ##   Hours  ctrS ctrS.1 ctrS.2  ctrN ctrN.1 ctrN.2  nd5S nd5S.1 nd5S.2  nd5s
    ## 1    15  2870   2754   2915  2597   2145   2299  1784   1787   1641  2012
    ## 2    39  5683   5907   7048  6154   6485   6533  4655   3692   4342  4837
    ## 3    69 18614  17164  18843    NA  24128  20128 11219  10866   9870 11622
    ## 4   118 45642  55635  48772 50272  46718  47451 31481  33870  39539 40080
    ## 5   136 57573  62054  66285 60681  67081  61176 32097  35036  46219 54855
    ##   nd5s.1 nd5s.2  ndsS ndsS.1 ndsS.2  ndss ndss.1 ndss.2
    ## 1   3095   2935  1597   1617   1591  1622   1663   1737
    ## 2   5141   5498  3627   4020   3777  4406   4862   4710
    ## 3  11724  12820  9310  19471  10874 11567  11113  12147
    ## 4  42817  45711 29017  33048  33703 34044  37801  36009
    ## 5  42588  47550 34103  46421  44557 45800  55774  39656

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 4)  # get the first two characters
  if (prefix == 'ctrS') {
    return('ctr_s+')
  } else if (prefix == 'ctrN') {
    return('ctr_s-') #no supplements
  } else if (prefix == 'nd5S') {
    return('nd5_s+')
  } else if (prefix == 'nd5s') {
    return('nd5_s+U-') #S-U, no uridine
  }  else if (prefix == 'ndsS') {
    return('nd5_suspension_s+')
  } else if (prefix == 'ndss') {
    return('nd5_suspension_s+U-') #suspension, no uridine, supplemented
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)

##############################################################################
###split the treatment column into 4 columns
## nd5, other supplements, uridine, suspension


library(dplyr)

# Assuming your data frame is named df
# Update the name accordingly if it differs
df_long <- df_long %>%
  mutate(
    treatment = as.character(treatment),
    # Create the 'othersup' column based on whether 'treatment' ends with '+' or '-'
    othersup = case_when(
      grepl("\\+", treatment) ~ "plus",
      grepl("\\-", treatment) ~ "minus",
      TRUE ~ NA_character_
    ),

    # Create the 'nd5' column based on whether 'treatment' starts with 'c'
    nd5 = case_when(
      startsWith(treatment, "c") ~ "ctr",
      TRUE ~ "nd5"
    ),
    
    
    uridine = case_when(
      grepl("s\\+U-$", treatment) ~ "minus",
      grepl("s\\+", treatment) ~ "plus",
      
      TRUE ~ "minus"
    ),
    
    suspension = case_when(
      grepl("suspension", treatment) ~ "sus",
      TRUE ~ "no_sus"
    )
  )

# View the updated dataframe
#head(df_long)

#change the variable columns back to factor format
# Convert 'treatment', 'nd5', and 'tgfb1' to factors
df_long <- df_long %>%
  mutate(
    treatment = as.factor(treatment),
    nd5 = as.factor(nd5),
    othersup = as.factor(othersup),
    suspension = as.factor(suspension),
    uridine = as.factor(uridine)
  )

# Verify the changes
#str(df_long)






#write_csv(df_long, "052721_12A_sup_growth_rate.csv")
# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.3926 -0.1739  0.0030  0.1227  0.8357 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                         7.6347878  0.1187368  64.300   <2e-16 ***
    ## treatmentctr_s+                     0.1060735  0.1667246   0.636   0.5265    
    ## treatmentnd5_s+                    -0.3596574  0.1667246  -2.157   0.0341 *  
    ## treatmentnd5_s+U-                  -0.0531904  0.1667246  -0.319   0.7506    
    ## treatmentnd5_suspension_s+         -0.4315051  0.1667246  -2.588   0.0115 *  
    ## treatmentnd5_suspension_s+U-       -0.3708255  0.1667246  -2.224   0.0291 *  
    ## Hours                               0.0265226  0.0013274  19.981   <2e-16 ***
    ## treatmentctr_s+:Hours              -0.0009597  0.0018766  -0.511   0.6105    
    ## treatmentnd5_s+:Hours              -0.0008312  0.0018766  -0.443   0.6591    
    ## treatmentnd5_s+U-:Hours            -0.0017042  0.0018766  -0.908   0.3666    
    ## treatmentnd5_suspension_s+:Hours    0.0001109  0.0018766   0.059   0.9530    
    ## treatmentnd5_suspension_s+U-:Hours  0.0003902  0.0018766   0.208   0.8358    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2354 on 77 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.9685, Adjusted R-squared:  0.964 
    ## F-statistic: 215.3 on 11 and 77 DF,  p-value: < 2.2e-16

``` r
# ###testing emmeans
# library(emmeans)
# emmeans_result <- emmeans(model, pairwise ~ treatment | Hours)
# summary(emmeans_result)

##############################################################################
#####fitting each factor######################################################


model <- lm(log(cell_count) ~ othersup * Hours + nd5* Hours + uridine * Hours + suspension * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ othersup * Hours + nd5 * Hours + 
    ##     uridine * Hours + suspension * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41492 -0.16978  0.02297  0.11689  0.87726 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          7.6347878  0.1186593  64.342   <2e-16 ***
    ## othersupplus         0.2896468  0.2035722   1.423   0.1587    
    ## Hours                0.0265226  0.0013265  19.994   <2e-16 ***
    ## nd5nd5              -0.4042840  0.1547301  -2.613   0.0107 *  
    ## uridineplus         -0.1835733  0.1169649  -1.569   0.1205    
    ## suspensionsus       -0.1947414  0.1169649  -1.665   0.0999 .  
    ## othersupplus:Hours  -0.0012565  0.0022966  -0.547   0.5858    
    ## Hours:nd5nd5        -0.0001596  0.0017536  -0.091   0.9277    
    ## Hours:uridineplus    0.0002968  0.0013256   0.224   0.8234    
    ## Hours:suspensionsus  0.0015182  0.0013256   1.145   0.2555    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2353 on 79 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.9677, Adjusted R-squared:  0.9641 
    ## F-statistic: 263.3 on 9 and 79 DF,  p-value: < 2.2e-16

``` r
#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["nd5_sup_noUridine"]])
```

    ## NULL

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = treatments)

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](cell_growth_comparision_files/figure-gfm/052721_12A_sup-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##                          slope           se
    ## nd5_suspension_s+U- 0.02691281 0.0011865280
    ## nd5_suspension_s+   0.02663351 0.0016497731
    ## nd5_s+U-            0.02481845 0.0009626156
    ## nd5_s+              0.02569142 0.0012802895
    ## ctr_s+              0.02556295 0.0011496723
    ## ctr_s-              0.02652262 0.0016151681

``` r
###########################################################
##split treatments for linear regression#################
# Load the readxl package
library(readxl)

# File path - replace this with the path to your file
file_path <- "052721_12A_sup_growth_rate_split_treatments.xls"

# Read the Excel file
data <- read_excel(file_path)
data$cell_count <- as.numeric(data$cell_count)
```

    ## Warning: NAs introduced by coercion

``` r
# drop NAs
data <- na.omit(data, cols = "cell_count")

# View the first few rows of the data
model <- lm(log(cell_count) ~ nd5 * Hours + other_sup* Hours + Uridine*Hours + suspension*Hours , data = data)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ nd5 * Hours + other_sup * Hours + 
    ##     Uridine * Hours + suspension * Hours, data = data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41492 -0.16978  0.02297  0.11689  0.87726 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       7.6347878  0.1186593  64.342   <2e-16 ***
    ## nd5              -0.4042840  0.1547301  -2.613   0.0107 *  
    ## Hours             0.0265226  0.0013265  19.994   <2e-16 ***
    ## other_sup         0.2896468  0.2035722   1.423   0.1587    
    ## Uridine          -0.1835733  0.1169649  -1.569   0.1205    
    ## suspension       -0.1947414  0.1169649  -1.665   0.0999 .  
    ## nd5:Hours        -0.0001596  0.0017536  -0.091   0.9277    
    ## Hours:other_sup  -0.0012565  0.0022966  -0.547   0.5858    
    ## Hours:Uridine     0.0002968  0.0013256   0.224   0.8234    
    ## Hours:suspension  0.0015182  0.0013256   1.145   0.2555    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2353 on 79 degrees of freedom
    ## Multiple R-squared:  0.9677, Adjusted R-squared:  0.9641 
    ## F-statistic: 263.3 on 9 and 79 DF,  p-value: < 2.2e-16

### Fitting 052721_hela_sup

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
#df <- read.table("031123_12a_d10.txt", header=TRUE)

#There is NA in the row
df <- read.delim("052721_hela_sup.txt", header = TRUE, na.strings = "", fill = TRUE)
```

    ## Warning in read.table(file = file, header = header, sep = sep, quote = quote, :
    ## incomplete final line found by readTableHeader on '052721_hela_sup.txt'

``` r
head(df)
```

    ##   Hours  ctrS ctrS.1 ctrS.2  ctrN ctrN.1 ctrN.2  nd5S nd5S.1 nd5S.2 nd5s nd5s.1
    ## 1    24  2108   2048   2125  1988   1996   1919  1250   1489   1444 1401   1375
    ## 2    69  8089   4416   6111  5468   6580   5458  4663   5158   4770 5146   5147
    ## 3    93 14283  15376  14796 10180  11669  11071  9243  10220   9794 7521   8358
    ## 4   117 17813  22906  23158 14790  14790  14583 10727  12374  12737 8965  10011
    ##   nd5s.2
    ## 1   1554
    ## 2   5435
    ## 3   7933
    ## 4   9648

``` r
# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 4)  # get the first two characters
  if (prefix == 'ctrS') {
    return('ctr_s+')
  } else if (prefix == 'ctrN') {
    return('ctr_no_sup')
  } else if (prefix == 'nd5S') {
    return('nd5_supplemented')
  } else if (prefix == 'nd5s') {
    return('nd5_sup_noUridine')
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)


# Create a supplement variable
assign_sup <- function(replicate) {
  prefix <- substr(replicate, 1, 4)  # get the first two characters
  if (prefix == 'ctrS' | prefix == 'nd5S') {
    return('sup+')
  } else if (prefix == 'ctrN' | prefix == 'nd5s') {
    return('sup-')
  } else {
    return('unknown')  # or whatever default value you want
  }
}
df_long$sup <- sapply(df_long$replicate, assign_sup)

#write.table(df_long, "df_long.txt", sep = "\t", row.names = FALSE)


# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41056 -0.10828 -0.01042  0.14248  0.24701 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       7.098331   0.120250  59.030   <2e-16 ***
    ## treatmentctr_s+                  -0.077336   0.170059  -0.455   0.6517    
    ## treatmentnd5_sup_noUridine       -0.185652   0.170059  -1.092   0.2815    
    ## treatmentnd5_supplemented        -0.348869   0.170059  -2.051   0.0468 *  
    ## Hours                             0.022302   0.001446  15.427   <2e-16 ***
    ## treatmentctr_s+:Hours             0.003532   0.002044   1.728   0.0918 .  
    ## treatmentnd5_sup_noUridine:Hours -0.001419   0.002044  -0.694   0.4916    
    ## treatmentnd5_supplemented:Hours   0.001783   0.002044   0.872   0.3883    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1721 on 40 degrees of freedom
    ## Multiple R-squared:  0.965,  Adjusted R-squared:  0.9588 
    ## F-statistic: 157.3 on 7 and 40 DF,  p-value: < 2.2e-16

``` r
##test effects of supplements
model <- lm(log(cell_count) ~ sup * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ sup * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52500 -0.17595  0.04321  0.14758  0.43408 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    7.005504   0.122751  57.071   <2e-16 ***
    ## supsup+       -0.120276   0.173597  -0.693    0.492    
    ## Hours          0.021593   0.001476  14.632   <2e-16 ***
    ## supsup+:Hours  0.003367   0.002087   1.613    0.114    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2484 on 44 degrees of freedom
    ## Multiple R-squared:  0.9197, Adjusted R-squared:  0.9142 
    ## F-statistic: 167.9 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
# # ###testing emmeans
#  library(emmeans)
#  emmeans_result <- emmeans(model )
#  summary(emmeans_result)


## ANOCOVA on sup
model <- lm(log(cell_count) ~ Hours * sup, data = df_long)

summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * sup, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52500 -0.17595  0.04321  0.14758  0.43408 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    7.005504   0.122751  57.071   <2e-16 ***
    ## Hours          0.021593   0.001476  14.632   <2e-16 ***
    ## supsup+       -0.120276   0.173597  -0.693    0.492    
    ## Hours:supsup+  0.003367   0.002087   1.613    0.114    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2484 on 44 degrees of freedom
    ## Multiple R-squared:  0.9197, Adjusted R-squared:  0.9142 
    ## F-statistic: 167.9 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
# Get the ANCOVA table
anova(model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(cell_count)
    ##           Df  Sum Sq Mean Sq  F value Pr(>F)    
    ## Hours      1 30.7043 30.7043 497.5388 <2e-16 ***
    ## sup        1  0.2180  0.2180   3.5328 0.0668 .  
    ## Hours:sup  1  0.1606  0.1606   2.6031 0.1138    
    ## Residuals 44  2.7153  0.0617                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## add ctr/nd5 column

df_long$mutation  <- ifelse(startsWith(df_long$replicate, "ctr"), "ctr", "nd5")
## ANOCOVA on nd5 mutation only
model <- lm(log(cell_count) ~ Hours * mutation, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * mutation, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36063 -0.11415 -0.00598  0.16350  0.34256 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        7.059663   0.094237  74.914   <2e-16 ***
    ## Hours              0.024068   0.001133  21.244   <2e-16 ***
    ## mutationnd5       -0.228593   0.133272  -1.715   0.0933 .  
    ## Hours:mutationnd5 -0.001584   0.001602  -0.989   0.3283    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1907 on 44 degrees of freedom
    ## Multiple R-squared:  0.9526, Adjusted R-squared:  0.9494 
    ## F-statistic: 295.1 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
anova(model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(cell_count)
    ##                Df  Sum Sq Mean Sq  F value    Pr(>F)    
    ## Hours           1 30.7043 30.7043 844.1781 < 2.2e-16 ***
    ## mutation        1  1.4581  1.4581  40.0886 1.098e-07 ***
    ## Hours:mutation  1  0.0355  0.0355   0.9774    0.3283    
    ## Residuals      44  1.6004  0.0364                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
###fitted at the same time for mutatoin and sup########
model <- lm(log(cell_count) ~ Hours * mutation * sup, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * mutation * sup, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41056 -0.10828 -0.01042  0.14248  0.24701 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                7.0983305  0.1202499  59.030   <2e-16 ***
    ## Hours                      0.0223022  0.0014457  15.427   <2e-16 ***
    ## mutationnd5               -0.1856521  0.1700591  -1.092   0.2815    
    ## supsup+                   -0.0773358  0.1700591  -0.455   0.6517    
    ## Hours:mutationnd5         -0.0014192  0.0020445  -0.694   0.4916    
    ## Hours:supsup+              0.0035320  0.0020445   1.728   0.0918 .  
    ## mutationnd5:supsup+       -0.0858811  0.2404998  -0.357   0.7229    
    ## Hours:mutationnd5:supsup+ -0.0003296  0.0028914  -0.114   0.9098    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1721 on 40 degrees of freedom
    ## Multiple R-squared:  0.965,  Adjusted R-squared:  0.9588 
    ## F-statistic: 157.3 on 7 and 40 DF,  p-value: < 2.2e-16

``` r
#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["ctr_no_sup"]])
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41056 -0.10828 -0.01042  0.14248  0.24701 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                       7.098331   0.120250  59.030   <2e-16 ***
    ## treatmentctr_s+                  -0.077336   0.170059  -0.455   0.6517    
    ## treatmentnd5_sup_noUridine       -0.185652   0.170059  -1.092   0.2815    
    ## treatmentnd5_supplemented        -0.348869   0.170059  -2.051   0.0468 *  
    ## Hours                             0.022302   0.001446  15.427   <2e-16 ***
    ## treatmentctr_s+:Hours             0.003532   0.002044   1.728   0.0918 .  
    ## treatmentnd5_sup_noUridine:Hours -0.001419   0.002044  -0.694   0.4916    
    ## treatmentnd5_supplemented:Hours   0.001783   0.002044   0.872   0.3883    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1721 on 40 degrees of freedom
    ## Multiple R-squared:  0.965,  Adjusted R-squared:  0.9588 
    ## F-statistic: 157.3 on 7 and 40 DF,  p-value: < 2.2e-16

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = treatments)

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/052721_hela_sup-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##                        slope          se
    ## nd5_supplemented  0.02408537 0.001457344
    ## nd5_sup_noUridine 0.02088295 0.001623257
    ## ctr_s+            0.02583414 0.001637047
    ## ctr_no_sup        0.02230215 0.000959781

### Fitting 052721_hek293t_sup

``` r
##all libraries
library(ggplot2)
library(emmeans)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
# Load the data
df <- read.table("052721_hek293t_sup.txt", header=TRUE)
```

    ## Warning in read.table("052721_hek293t_sup.txt", header = TRUE): incomplete
    ## final line found by readTableHeader on '052721_hek293t_sup.txt'

``` r
# #There is NA in the row
# df <- read.delim("052721_hek293t_sup.txt", header = TRUE, na.strings = "", fill = TRUE)
# head(df)



# Reshape the data to long format
df_long <- reshape(df, varying=names(df)[-1], v.names="cell_count", 
                   timevar="replicate", times=names(df)[-1], direction="long")


# Create a treatment variable
assign_treatment <- function(replicate) {
  prefix <- substr(replicate, 1, 4)  # get the first two characters
  if (prefix == 'ctrS') {
    return('ctr_s+')
  } else if (prefix == 'ctrN') {
    return('ctr_no_sup')
  } else if (prefix == 'nd5S') {
    return('nd5_supplemented')
  } else if (prefix == 'nd5s') {
    return('nd5_no_sup')
  } else {
    return('unknown')  # or whatever default value you want
  }
}

df_long$treatment <- sapply(df_long$replicate, assign_treatment)


# Create a supplement variable
assign_sup <- function(replicate) {
  prefix <- substr(replicate, 1, 4)  # get the first two characters
  if (prefix == 'ctrS' | prefix == 'nd5S') {
    return('sup+')
  } else if (prefix == 'ctrN' | prefix == 'nd5s') {
    return('sup-')
  } else {
    return('unknown')  # or whatever default value you want
  }
}
df_long$sup <- sapply(df_long$replicate, assign_sup)

#write.table(df_long, "df_long.txt", sep = "\t", row.names = FALSE)


# # Fit the mixed model
# model <- lmer(log(cell_count) ~ treatment + (1|replicate), data=df_long)
# # Get the summary of the model
# summary(model)
# 
# ##results indicate random effect is not necessary



model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.62469 -0.14672 -0.00113  0.13828  0.34571 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                      6.1725378  0.1566205  39.411  < 2e-16 ***
    ## treatmentctr_s+                  0.0845382  0.2214948   0.382  0.70473    
    ## treatmentnd5_no_sup             -0.5655067  0.2214948  -2.553  0.01459 *  
    ## treatmentnd5_supplemented       -0.3807685  0.2214948  -1.719  0.09333 .  
    ## Hours                            0.0402959  0.0018829  21.401  < 2e-16 ***
    ## treatmentctr_s+:Hours            0.0003217  0.0026629   0.121  0.90446    
    ## treatmentnd5_no_sup:Hours        0.0080719  0.0026629   3.031  0.00426 ** 
    ## treatmentnd5_supplemented:Hours  0.0076389  0.0026629   2.869  0.00655 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2241 on 40 degrees of freedom
    ## Multiple R-squared:  0.9824, Adjusted R-squared:  0.9794 
    ## F-statistic: 319.5 on 7 and 40 DF,  p-value: < 2.2e-16

``` r
##test effects of supplements
model <- lm(log(cell_count) ~ sup * Hours, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ sup * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.62936 -0.14015  0.00314  0.13183  0.45478 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    5.890e+00  1.271e-01  46.346   <2e-16 ***
    ## supsup+        1.346e-01  1.797e-01   0.749    0.458    
    ## Hours          4.433e-02  1.528e-03  29.016   <2e-16 ***
    ## supsup+:Hours -5.565e-05  2.161e-03  -0.026    0.980    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2572 on 44 degrees of freedom
    ## Multiple R-squared:  0.9745, Adjusted R-squared:  0.9728 
    ## F-statistic: 561.6 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
# # ###testing emmeans
#  library(emmeans)
#  emmeans_result <- emmeans(model )
#  summary(emmeans_result)


## ANOCOVA on sup
model <- lm(log(cell_count) ~ Hours * sup, data = df_long)

summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * sup, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.62936 -0.14015  0.00314  0.13183  0.45478 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    5.890e+00  1.271e-01  46.346   <2e-16 ***
    ## Hours          4.433e-02  1.528e-03  29.016   <2e-16 ***
    ## supsup+        1.346e-01  1.797e-01   0.749    0.458    
    ## Hours:supsup+ -5.565e-05  2.161e-03  -0.026    0.980    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2572 on 44 degrees of freedom
    ## Multiple R-squared:  0.9745, Adjusted R-squared:  0.9728 
    ## F-statistic: 561.6 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
# Get the ANCOVA table
anova(model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(cell_count)
    ##           Df  Sum Sq Mean Sq   F value  Pr(>F)    
    ## Hours      1 111.241 111.241 1681.7924 < 2e-16 ***
    ## sup        1   0.204   0.204    3.0860 0.08593 .  
    ## Hours:sup  1   0.000   0.000    0.0007 0.97957    
    ## Residuals 44   2.910   0.066                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## add ctr/nd5 column

df_long$mutation  <- ifelse(startsWith(df_long$replicate, "ctr"), "ctr", "nd5")
## ANOCOVA on nd5 mutation only
model <- lm(log(cell_count) ~ Hours * mutation, data = df_long)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ Hours * mutation, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.70213 -0.12573 -0.00357  0.13622  0.39439 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        6.214807   0.111018  55.980  < 2e-16 ***
    ## Hours              0.040457   0.001335  30.312  < 2e-16 ***
    ## mutationnd5       -0.515407   0.157003  -3.283 0.002019 ** 
    ## Hours:mutationnd5  0.007695   0.001888   4.077 0.000188 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2247 on 44 degrees of freedom
    ## Multiple R-squared:  0.9806, Adjusted R-squared:  0.9793 
    ## F-statistic: 740.5 on 3 and 44 DF,  p-value: < 2.2e-16

``` r
anova(model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(cell_count)
    ##                Df  Sum Sq Mean Sq   F value    Pr(>F)    
    ## Hours           1 111.241 111.241 2203.7315 < 2.2e-16 ***
    ## mutation        1   0.055   0.055    1.0818 0.3039808    
    ## Hours:mutation  1   0.839   0.839   16.6181 0.0001885 ***
    ## Residuals      44   2.221   0.050                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#######################all pairs###################################
# Get the levels of the treatment factor
df_long$treatment <- as.factor(df_long$treatment)
treatments <- levels(df_long$treatment)

# Initialize a list to store the model summaries
model_summaries <- list()

# Loop over the treatments
for (treat in treatments) {
  # Relevel the treatment factor so that the current treatment is the reference category
  df_long$treatment <- relevel(df_long$treatment, ref = treat)
  
  # Fit the model
  model <- lm(log(cell_count) ~ treatment * Hours, data = df_long)
  
  # Store the summary of the model in the list
  model_summaries[[treat]] <- summary(model)
}

# Now model_summaries is a list of model summaries, one for each treatment. You can access 
# them using model_summaries[["ctr-"]], model_summaries[["ctr+"]], etc.

print(model_summaries[["nd5_no_sup"]])
```

    ## 
    ## Call:
    ## lm(formula = log(cell_count) ~ treatment * Hours, data = df_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.62469 -0.14672 -0.00113  0.13828  0.34571 
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                      5.607031   0.156620  35.800  < 2e-16 ***
    ## treatmentctr_s+                  0.650045   0.221495   2.935  0.00551 ** 
    ## treatmentctr_no_sup              0.565507   0.221495   2.553  0.01459 *  
    ## treatmentnd5_supplemented        0.184738   0.221495   0.834  0.40921    
    ## Hours                            0.048368   0.001883  25.687  < 2e-16 ***
    ## treatmentctr_s+:Hours           -0.007750   0.002663  -2.910  0.00587 ** 
    ## treatmentctr_no_sup:Hours       -0.008072   0.002663  -3.031  0.00426 ** 
    ## treatmentnd5_supplemented:Hours -0.000433   0.002663  -0.163  0.87166    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2241 on 40 degrees of freedom
    ## Multiple R-squared:  0.9824, Adjusted R-squared:  0.9794 
    ## F-statistic: 319.5 on 7 and 40 DF,  p-value: < 2.2e-16

``` r
###########################################################################################

# Create a new data frame for predictions
newdata <- expand.grid(Hours = seq(min(df_long$Hours), max(df_long$Hours), length.out = 100),
                       treatment = treatments)

# Generate predictions
newdata$pred <- predict(model, newdata)

# Convert treatment to factor for plotting
newdata$treatment <- as.factor(newdata$treatment)

# Plot the observed data
ggplot(df_long, aes(x = Hours, y = log(cell_count), color = treatment)) +
  geom_point() +
  # Add the model predictions
  geom_line(data = newdata, aes(x = Hours, y = pred)) +
  labs(x = "Time (hours)", y = "Log-transformed cell count") +
  theme_classic()
```

![](cell_growth_comparision_files/figure-gfm/052721_hek293t_sup-1.png)<!-- -->

``` r
##calculate se and slope when fitted seorately

# Fit a model to each treatment
models <- df_long %>%
  split(.$treatment) %>%
  map(~ lm(log(cell_count) ~ Hours, data = .))

# Tidy up the results
tidied <- map_df(models, tidy, .id = "treatment")

#library(emmeans)

# Compute the EMMs
#emm <- emmeans(tidied, "time..hour.", by = "treatment")


# Compute the slopes
slopes <- map(models, ~ coef(.)["Hours"])

# Convert the slopes to a data frame
slopes_df <- as.data.frame(do.call(rbind, slopes))
names(slopes_df) <- "slope"

# Compute the standard errors of the slopes
se <- map(models, ~ sqrt(vcov(.)["Hours", "Hours"]))

# Convert the standard errors to a data frame
se_df <- as.data.frame(do.call(rbind, se))
names(se_df) <- "se"

# Combine the slopes and standard errors into one data frame
results <- cbind(slopes_df, se_df)

# Print the results
print(results)
```

    ##                       slope          se
    ## nd5_supplemented 0.04793487 0.001671202
    ## nd5_no_sup       0.04836783 0.002429831
    ## ctr_s+           0.04061760 0.001816802
    ## ctr_no_sup       0.04029594 0.001477829
