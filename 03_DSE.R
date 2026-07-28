# ----------------------------------------------------------------------------------------

# syntax for calculation of predicted probability that an incident is known to the police
# ... from self-reported victimization in the General Social Survey (GSS), 2019

# this syntax calculates the point estimates and confidence intervals by:
# ... 1) dual systems estimation procedure and Seber's approximate variance formula
# ... 2) logistic model to calculate point estimates and confidence intervals

# ----------------------------------------------------------------------------------------



# don't run
# install packages used in this script
# install.packages(
    # c('utils',
    #   'tibble',
    #   'stats',
    #   'broom',
    #   'magrittr',
    #   'ggplot2'
    #   )

# call pipe function
`%>%` <- magrittr::`%>%`

# download the file
utils::download.file(
  url = "https://borealisdata.ca/api/access/datafile/661489?format=RData", # define url
  destfile = "borealis_data.RData", # target file name
  mode = "wb"
  )

# load the data
temp <- new.env() # temp file
load("borealis_data.RData", envir = temp)
gss <- temp[[ls(temp)[1]]]

# recode values for variable that asks respondent whether police know about the incident?
gss$PFO_100[gss$PFO_100 > 1] <- 0 # recode reponse items = 0 as they're not known to the police

# define parameters for estimation -------------------------------------------------------------

# number of police-reported incidents according to Statistics Canada, 2019
# https://www150.statcan.gc.ca/n1/pub/85-002-x/2020001/article/00010-eng.htm
police <- 2.2E6

# GSS respondent responses
table(gss$PFO_100)

# sample size of the general social survey
victim <- nrow(gss)

# let's assume 99% of people are part of both samples
recapture <- table(gss$PFO_100)[2]



# 1) the arithmetic to calculate the probability any given incident is known to the police -------------------------------------------------------------

# dual systems estimation a.k.a. Lincoln-Petersen formula
N_hat <- (police * victim) / recapture

# Seber’s variance formula
# https://books.google.co.uk/books/about/Estimation_of_Animal_Abundance.html?id=iIAAPQAACAAJ&redir_esc=y
v <- (police^2 * victim * (victim - recapture)) / (recapture^3)

# standard error
s <- sqrt(v)

# critical value 
z <- stats::qnorm(0.975)

# confidence intervals
N_lo <- N_hat - z * s; N_hi <- N_hat + z * s

# print results
results <- tibble::tibble(
  statistics = c(
    "police-reported incidents, 2019",
    "estimated incidents",
    "standard error of the estimate",
    "95% CI lower bound",
    "95% CI upper bound",
    "incidents known to police, point estimate (%)",
    "incidents known to police, lower bound (%)",
    "incidents known to police, upper bound (%)"
    ),
  estimates = c(
    police,
    round(N_hat, 2),
    round(s, 2),
    round(N_lo, 2),
    round(N_hi, 2),
    round((police / N_hat) * 100, 2),
    round((police / N_lo) * 100, 2),
    round((police / N_hi) * 100, 2)
    )
  )
print(results)



# 2) logistic model to estimate the probability any given incident is known to police ------------------------------------------------------------- 
m_00 <- stats::glm(PFO_100 ~ 1, family = binomial(link = "logit"), data = gss)

# inverse logit function to calculate baseline probability that an incident is known to police
inv_logit <- function(x) { return(1 / (1 + exp(-x))) }

# calculate the baseline probability 
inv_logit(broom::tidy(m_00)$estimate)

# lower bound of the confidence interval
inv_logit(broom::tidy(m_00)$estimate - (stats::qnorm(0.975) * broom::tidy(m_00)$std.error))
  
# upper bound of the confidence interval
inv_logit(broom::tidy(m_00)$estimate + (stats::qnorm(0.975) * broom::tidy(m_00)$std.error))

# note: notice the arithmetic and the logistic regression result in near identical estimates



# 3) add covariate information into the logistic model to increase the precision of the estimate ---------------------------------------------------------
# note: labels for the numeric values for all variables can be found online at: https://borealisdata.ca/data-explorer/?siteUrl=https:%2F%2Fborealisdata.ca&dvLocale=en&fileId=661489

# type of crime/incident ---------------------------------------------

# first, group alike categories

  # attempted theft
  gss$CIR_D010[gss$CIR_D010 == 9] <- 2 
  
  # stolen property
  gss$CIR_D010[gss$CIR_D010 == 5] <- 4 
  gss$CIR_D010[gss$CIR_D010 == 6] <- 4
  
  # sexual assualt and related crimes
  gss$CIR_D010[gss$CIR_D010 == 13] <- 12 
  gss$CIR_D010[gss$CIR_D010 == 14] <- 12

# add labels to numeric values
gss$CIR_D010 <- factor(
  gss$CIR_D010,
  levels = c(1, 2, 3, 4, 7, 8 , 10, 11, 12, 15),
  labels = c(
    "household damage",
    "attemt to take or steal something",
    "attempt to break and enter",
    "stolen property",
    "theft or attempted theft of vehicle",
    "damange to vehicles",
    "assault",
    "threat of assault",
    "sexual assault",
    "another crime"
    )
  )
  
  # specify reference category
  gss$CIR_D010 <- relevel(gss$CIR_D010, ref = "stolen property")
  

# recode values for variable asking respondents about their relationship to the offender(s) ---------------------------------------------
gss$OFFENDRC[gss$OFFENDRC == 9] <- 6

# add labels to numeric values
gss$OFFENDRC <- factor(
  gss$OFFENDRC,
  levels = c(1, 2, 3, 4, 6),
  labels = c(
    "relative",
    "friend/neighbour/acquaintance",
    "stranger",
    "other type of relation",
    "valid skip i.e., no relation"
    )
  )

# specify reference category
gss$OFFENDRC <- relevel(gss$OFFENDRC, ref = "valid skip i.e., no relation")


# location where the incident happened ---------------------------------------------
gss$WHR_100[gss$WHR_100 == 9] <- 5

# add labels to numeric values
gss$WHR_100 <- factor(
  gss$WHR_100,
  levels = c(1, 2, 3, 4, 5),
  labels = c(
    "home/nearby",
    "other private residence",
    "commerical or institutional establishment",
    "street or other public place",
    "other place"
    )
  )

# specify the reference catefory
gss$WHR_100 <- relevel(gss$WHR_100, ref = "home/nearby")


# recode categorical variables that incidate if respondent / victom spoke with anyone about the incident ---------------------------------------------

  # spoke with family member
  gss$TTA_110[gss$TTA_110 == 1] <- 1
  gss$TTA_110[gss$TTA_110  > 1] <- 0
  
  # spoke with friend or neighbour
  gss$TTA_120[gss$TTA_120 == 1] <- 1
  gss$TTA_120[gss$TTA_120  > 1] <- 0
  
  # spoke with co-worker
  gss$TTA_130[gss$TTA_130 == 1] <- 1
  gss$TTA_130[gss$TTA_130  > 1] <- 0

  # spoke with doctor or nurse
  gss$TTA_140[gss$TTA_140 == 1] <- 1
  gss$TTA_140[gss$TTA_140  > 1] <- 0
  
  # spoke with lawyer
  gss$TTA_150[gss$TTA_150 == 1] <- 1
  gss$TTA_150[gss$TTA_150  > 1] <- 0
  
  # spoke with priest
  gss$TTA_160[gss$TTA_160 == 1] <- 1
  gss$TTA_160[gss$TTA_160  > 1] <- 0


# logistic model to estimate the probability is known to police , controlling for incident characteristics 
m_01 <- stats::glm(
  PFO_100 ~ 
    CIR_D010 + # type of incident 
    OFFENDRC + # victim-offender relationship
    WHR_100 + # location of incident
    # did the respondent/victim talk to anyone about the incident?
    TTA_110 + # spoke with family member
    TTA_120 + # spoke with friend or neighbour
    TTA_130 + # spoke with co-worker
    TTA_140 + # spoke with doctor or nurse
    TTA_150 + # spoke with lawyer
    TTA_160 ,  # spoke with priest
  family = binomial(link = "logit"), 
  data = gss
  )

# geometric average predicted probability, adjusts for skewness in probability distribution
p_y1 <- exp(mean(log(predict(m_01, type = "response")[m_01$y == 1])))
print(p_y1)

# function to calculate Tjur's R-squared ---------------------------------------
tjur <- function(model){
  
  # predicted probabilities
  pp <- predict(model, type = "response")
  
  # internal function to calculate geometric mean
  geomean <- function(x) {
    exp(mean(log(x[x > 0]))) # truncate to prevent log(0) errors
  }
  
  # calculate the modified Tjur's R-squared
  rsq <- geomean(pp[model$y == 1]) - geomean(pp[model$y == 0])
  
  message("R-squared ranges 0-1, where higher values indicate models that have more predictive power.")
  return(rsq)
}
tjur(m_01)

# predicted probabilities
pp <- predict(m_01, type = "response", se.fit = TRUE)

# mean or point estimates
mu <- pp$fit

# calculate confidence intervals
  
  # lower bound
  lo <- pp$fit - stats::qnorm(0.975) * pp$se.fit
  
  # upper bound
  hi <- pp$fit + stats::qnorm(0.975) * pp$se.fit

# join the predicted probabilities to dataset
gss <- cbind(gss, mu, lo, hi)



# plot the predicted probability distribution that incident was known to the police by crime type ----------------------------------------
fig01 <- ggplot2::ggplot(gss, ggplot2::aes(x = CIR_D010, y = mu)) +
  ggplot2::geom_hline(yintercept = p_y1, color = "red", linetype = "solid", linewidth = 1, alpha = 0.4) +
  ggplot2::geom_jitter(color = "darkgrey", alpha = 0.2, width = 0.2) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(
    x = "Type of incident",
    y = "Predicted probability",
    title = "The predicted probability that an incident is known to the police by type of crime",
    caption = stringr::str_wrap(
      "Figure notes: The boxplots present the predicted probability that an incident is known to the police calculated from self-reports to the 2019 General Social Survey (GSS) on victimization (N=5,522). The red line indicates the geometric average predicted probability, which accounts for skewness in the probability distribution.",
      width = 150
      )
    ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0, size = 12, face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0, size = 10, face = "plain"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
print(fig01)

# path to folder 
path <- "~/Desktop"

# function to output high resolution images
output <- function(filename, figure, path = path, width = 15, height = 10){
  ggplot2::ggsave(
    filename,
    figure,
    path = path, 
    width = width, 
    height = height, 
    device = 'png', 
    dpi = 250 # larger DPI increases the size of the plot aesthetics 
    )
}
output(
  filename = "fig01.png",
  figure = fig01,
  path = path,
  height = 10,
  width = 10
  )





# emd .R script


