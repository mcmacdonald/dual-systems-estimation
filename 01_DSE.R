# ----------------------------------------------------------------------------------


# syntax for calculation of dual systems estimation

# ... this syntax calculates population estimates and confidence intervals by:
# ... 1) dual systems estimation procedure and Seber's approximate variance formula
# ... 2) log-linear model to calculate population estimates and confidence intervals


# ----------------------------------------------------------------------------------



# don't run
# install packages used in this script
# install.packages('dplyr'); install.packages('tibble'); install.packages('stats); install.packages('VGAM')



# set seed for reproducibility
set.seed(123)



# define parameters for estimation -------------------------------------------------------------
  
  # let's assume that the population size is 10,100
  P <- 10100
  
  # the census counts 10,000 people
  census <- 10000
  
  # a post-enumeration survey (PES) also surveys 10,000 people
  pes <- 10000
  
  # let's assume 99% of people are part of both samples
  recapture <- 9900

  
  
# 1) do the arithmetic to calculate the population estimates -------------------------------------------------------------
  
  # dual systems estimation a.k.a. Lincoln-Petersen formula
  P_hat <- (census * pes) / recapture
  
  # Seber’s variance formula
  # https://books.google.co.uk/books/about/Estimation_of_Animal_Abundance.html?id=iIAAPQAACAAJ&redir_esc=y
  v <- (census^2 * pes * (pes - recapture)) / (recapture^3)

  # standard error
  s <- sqrt(v)
  
  # confidence intervals
  ci_lower <- P_hat - 1.96 * s; ci_upper <- P_hat + 1.96 * s

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
    round(s, 2),
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
  ids <- 1:P
  
  # people named in both samples
  ids_matched <- sample(ids, recapture, replace = FALSE)
  
  # people named in one sample or the other
  ids_loners <- dplyr::setdiff(ids, ids_matched)
  
  # number of people who only appear in the census
  census_unique <- sample(ids_loners, census - recapture, replace = FALSE)
  
  # number of people who only appear in the post-enumeration survey
  pes_unique <- sample(ids_loners, pes - recapture, replace = FALSE)
  
  # full census sample
  census_sample <- c(ids_matched, census_unique)
  
  # full post-numeration survey
  pes_sample <- c(ids_matched, pes_unique)

# construct data frame
df <- data.frame(
  id = 1:P,
  in_census = ifelse(1:P %in% census_sample, 1, 0),
  in_pes = ifelse(1:P %in% pes_sample, 1, 0)
  )

# second, estimate the model and use predicted probabilities to calculate population estimate

  # the dependent variable(s) is the columns of captures
  y <- cbind(df$in_census, df$in_pes)
  
  # the log-linear model
  model <- VGAM::vglm(y ~ 1, family = VGAM::posbernoulli.t(parallel = TRUE ~ 1), data = df)

  # ... manually calculate predicted probabilities and confidence intervals
  pp <- VGAM::predict(model, se.fit = TRUE)
  
  # inverse logit function
  inv_logit <- function(x) exp(x) / (1 + exp(x))

  # point estimate and confidence intervals 
  
    # ... for the census 
    pp_census_point <- inv_logit(pp$fitted.values[, 1])
    pp_census_lower <- inv_logit(pp$fitted.values[, 1] - 1.96 * pp$se.fit[, 1])
    pp_census_upper <- inv_logit(pp$fitted.values[, 1] + 1.96 * pp$se.fit[, 1])
    
    # ... for the post-enumeration survey 
    pp_pes_point <- inv_logit(pp$fitted.values[, 2])
    pp_pes_lower <- inv_logit(pp$fitted.values[, 2] - 1.96 * pp$se.fit[, 2])
    pp_pes_upper <- inv_logit(pp$fitted.values[, 2] + 1.96 * pp$se.fit[, 2])

# calculate the joint probability
joint_probability <- function(list1, list2) {
  pp <- 1 - (1 - list1) * (1 - list2)
  return(pp)
  }

# Horvitz-Thompson estimator
nhat <- function(pp){
  nhat <- sum(1/pp)
  return(nhat)
}
P_hat <- nhat(joint_probability(pp_census_point, pp_pes_point))
P_hi <- nhat(joint_probability(pp_census_lower, pp_pes_lower))
P_lo <- nhat(joint_probability(pp_census_upper, pp_pes_upper))

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
    round(model@extra$SE.N.hat, 2),
    round(P_lo, 2),
    round(P_hi, 2),
    round((1 - (census / P_hat)) * 100, 2),
    round((1 - (census / P_lo)) * 100, 2),
    round((1 - (census / P_hi)) * 100, 2)
  )
)
print(results)



# ... close .R script
