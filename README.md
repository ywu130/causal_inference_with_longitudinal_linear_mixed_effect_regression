# causal_inference_with_longitudinal_linear_mixed_effect_regression

## Description

This analysis investigates the longitudinal growth rate of cells and identifies treatments that significantly influence cellular growth over time. A semi-logarithmic regression model with interaction terms was used to explore the combined effects of time (Hours), the presence of MT-ND5 (nd5) mutations, and TGFB1 treatment. Interaction terms with time were included to evaluate how these factors jointly affect growth trajectories. Treatments with significant interaction terms with time were identified as having a measurable impact on the cellular growth rate.

The script includes the following steps:
- Data organization and transformation
- Model fitting
- Model comparison
- Visualization of results across multiple datasets

**Dataset:**  
Cell growth data measured by Celigo, collected on 090622 from 293T, r8, and d25 samples.

## Regression Model

The regression model used to analyze the data is represented as follows:

$$\log_{10}(\text{cell count}) = \beta_0 + \beta_1 (\text{treatment}) + \beta_2 (\text{time}) + \beta_3 (\text{treatment} \times \text{time}) + \epsilon$$


Where:  
- **β0**: Intercept  
- **β1**: Effect of treatment  
- **β2**: Effect of time  
- **β3**: Interaction effect of treatment and time  
- **ε**: Random error term  

This model allows us to understand how treatments interact with time to influence the cell growth trajectory.

