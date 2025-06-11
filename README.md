# Monte Carlo study on the bias of regression coefficients in LogReg model

## Assignment statement

Conduct a Monte Carlo (MC) study to examine the influence of sample size and the type of distribution of predictor variables on the bias of regression coefficients in a dichotomous logistic regression model.

Perform the experiment for all combinations of the following factors:

1. Sample size: $n \in \{ 50, 100, 200 \}$;
2. Type of distribution of predictor variables:
    1. Normal, $\mathcal{N} (0, 1)$,
    2. Centered gamma, $\Gamma(0.5, 1) - \mathbf{E} (\Gamma(0.5, 1))$.
    3. Contaminated normal (90% of data from $\mathcal{N} (0, 2)$ and 10% from $\mathcal{N} (20, 5)$ ), and
    4. Uniform, on the interval $[-5,+5]$;
3. Number of predictor variables: 1 and 2.

True/population regression coefficients should be:

- In the model with one predictor: $\beta_{10} = 3$, $\beta_{11} = -3$
- In the model with two predictors (same type of distribution): $\beta_{20} = 4$, $\beta_{21} = -2$, $\beta_{22} = -2$

For each combination of factors, perform $n_{rep} = 1000$ replications such that for each replication:

1. Generate $n$ or $2n$ random numbers according to the specified distributions.
2. Calculate the linear combination with the given coefficients ($\eta$).
3. Generate $y$ according to the Bernoulli distribution with parameter calculated as $\frac{\exp(\eta)}{1 + \exp(\eta)}$.

For each combination of factors, calculate the bias as the difference between the average of the $n_{rep}$ estimated $\beta_{ij}$ values and the true/population $\beta_{ij}$.

Express the bias in relative terms as $\frac{\mathbf{E} (\hat{\beta}_{ij})- \beta_{ij}}{\beta_{ij}} \cdot 100 $, and compare them.
