This repository contains an illustrations of dual systems system.



The contents of the both files illustrate:

1) the arithmetic to manually calculate the population estimates, coverage error, and approximate confidence intervals; and

2) the calculation of the population, coverage error, and confidence intervals from a generalized linear model.


The contents of the first file ('01_DSE.R') does so for a small known population and a post-enumeration survey. It uses a vector generalized linear model to estimate the coverage error and, in turn, the population.

The contents of the second file ('02_DSE.R') does so for a large population (Canada's general population accoridng to the 2021 Census of Population) and a post-enumeration survey. It uses a standard generalized linear model to estimate the coverage error and, in turn, the population.

The generalized linear models can also include attributes of the respondents of the census and pos-enumeration surveys to adjust for population and/or sample heterogeneity that make recapture more or less likely.
