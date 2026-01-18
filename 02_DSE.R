# ----------------------------------------------------------------------------------


# syntax for calculation of dual systems estimation

# ... this syntax calculates population estimates and confidence intervals by:
# ... 1) dual systems estimation procedure and Seber's approximate variance formula
# ... 2) log-linear model to calculate population estimates and confidence intervals


# ----------------------------------------------------------------------------------



# don't run
# install packages used in this script
# install.packages(
# c('dplyr', 
# 'tibble', 
# 'stats, 
# 'VGAM'
# )



# set seed for reproducibility
set.seed(123)



# define parameters for estimation -------------------------------------------------------------


# Canadian population, per 2021 Census of Population
P <- 36991981

# let's assume the actual census counts 99.8% of people
known <- 0.998 # let's assume that the census in this case samples 99,8% of the general population
census <- P * known

# a post-enumeration survey (PES) of the population
target <- 0.04 # let's assume that PES in this case samples 10% of the general population
pes <- P * target

# let's also assume 99.8% of people are part of both samples
recapture <- pes * 0.998



# 1) do the arithmetic to calculate the population estimates -------------------------------------------------------------

# dual systems estimation a.k.a. Lincoln-Petersen formula
P_hat <- (census * pes) / recapture

# Seber’s variance formula
# for reference: https://books.google.co.uk/books/about/Estimation_of_Animal_Abundance.html?id=iIAAPQAACAAJ&redir_esc=y
v <- (census^2 * pes * (pes - recapture)) / (recapture^3)

# standard error
se <- sqrt(v)

# critical value
z <- 1.96

# confidence intervals
ci_lower <- P_hat - z * se; ci_upper <- P_hat + z * se

# print results
results <- tibble::tibble(
  statistics = c(
    "'true' population count",
    "estimated population count",
    "standard error of the estimate",
    "95% CI lower bound",
    "95% CI upper bound",
    "coverage error, point estimate (%)",
    "coverage error, lower bound (%)",
    "coverage error, upper bound (%)"
    ),
  estimates = c(
    P,
    round(P_hat, 2),
    round(se, 2),
    round(ci_lower, 2),
    round(ci_upper, 2),
    round((1 - (census / P_hat)) * 100, 2),
    round((1 - (census / ci_lower)) * 100, 2),
    round((1 - (census / ci_upper)) * 100, 2)
    )
  )
print(results)



# 2) use the log-linear model to calculate the estimates -------------------------------------------------------------

# first, construct the dataset for the demonstration

# generate list of unique ids for each person in the population
ids <- 1:pes

# people named in both samples
ids_recapture <- sample(ids, recapture, replace = FALSE)

# people named in one sample or the other
ids_unique <- dplyr::setdiff(ids, ids_recapture)

# recaptures
df <- tibble::tibble(
  id = ids,
  recapture = as.integer(ids %in% ids_recapture)
  )

# second, estimate the model and use predicted probabilities to calculate population estimate

# the log-linear model
model <- stats::glm(recapture ~ 1, family = binomial(link = "logit"), data = df)

# ... manually calculate predicted probabilities and confidence intervals

# intercept on the log odds scale
b0 <- summary(model)$coefficients[1]

# standard errors
se <- summary(model)$coefficients[2]

# inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# point estimate and confidence intervals 
pp_hat <- inv_logit(b0)
pp_lo <- inv_logit(b0 - z * se)
pp_hi <- inv_logit(b0 + z * se)

# print results
results <- tibble::tibble(
  statistics = c(
    "'true' population count",
    "estimated population count",
    "95% CI lower bound",
    "95% CI upper bound",
    "coverage error, point estimate (%)",
    "coverage error, lower bound (%)",
    "coverage error, upper bound (%)"
    ),
  estimates = c(
    P,
    census / pp_hat,
    census / pp_lo,
    census / pp_hi,
    round(pp_hat, digits = 4),
    round(pp_lo, digits = 4),
    round(pp_hi, digits = 4)
    )
  )
print(results)



# ... close .R script



