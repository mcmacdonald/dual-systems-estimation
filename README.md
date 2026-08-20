This repository contains a demonstration of dual systems estimation (DSE) to estimate coverage errors in:
1) an incomplete census; and
2) an administrative database (police-reported incidents).

The contents of the both folders illustrate:

1) the arithmetic to manually calculate the population estimates, coverage error, and approximate confidence intervals; and

2) the calculation of the population, coverage error, and confidence intervals from regression models

Folder contents for demo that estimates coverage errors in an incomplete census:

* The contents of the first file ('demo_01.R') uses a census and a post-enumeration survey to estimate a small known population. It uses a vector generalized linear model (VGLM) to estimate the coverage error and, in turn, recover the population estimate.

* The contents of the second file ('demo_02.R') uses a census and a post-enumeration survey to estimate a large population. To demonstrate, the models try to recover Canada's general population, according to the 2021 Census of Population. It uses a logistic model to estimate the coverage error and, in turn, recover the population estimate.

* The models can also include attributes of the respondents of the census and post-enumeration surveys to adjust for population and/or sample heterogeneity that make recapture more or less likely.

Folder contents for demo that estimates coverage errors in incidents known to police: 

* The contents of the third folder contains both R ('demo.R') and Python ('demo.py') scripts to estimate the probability that an incident is known to the police. To do so, it uses 2019 police-reported crime incidents reported by Statistics Canada and self-reported information about victimization from the 2019 General Social Survey on Victimization. The scripts estimate and visualize the probability distribution by type of crime, adjusting for attributes of the respondents, offenders, and incidents that make reporting/non-reporting to the police more or less likely.
