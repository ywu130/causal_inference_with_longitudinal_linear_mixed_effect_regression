# causal_inference_with_longitudinal_linear_mixed_effect_regression

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
