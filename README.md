# SSvM (Simulation)

complete, self-contained R script for your Simulation Study section. It:

Implements correct SSvM sampling (rejection using von Mises as proposal).

Fits SSvM by MLE with proper constraints via reparameterization.

Computes observed-information standard errors (via numerical Hessian) and 95% Wald CIs.

Runs all 16 scenarios: 
μ∈{π/4,π/2}, 
𝜅∈{2,5},
η∈{0,0.4}, 
n∈{200,1000}, with nsim = 1000 

Produces bias, RMSE, and coverage summaries + clean, colorful ggplot2 figures.

Saves everything under paper_outputs/simulations/<timestamp>/ and also prints the key tables and plots in the console for instant visual checks.
