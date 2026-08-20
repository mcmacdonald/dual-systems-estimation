# ----------------------------------------------------------------------------------------

# syntax for calculation of predicted probability that an incident is known to the police
# ... from self-reported victimization in the General Social Survey (GSS), 2019

# this syntax calculates the point estimates and confidence intervals by:
# ... 1) dual systems estimation procedure and Seber's approximate variance formula
# ... 2) logistic model to calculate point estimates and confidence intervals

# *** note: the input file (gss_34_incident_pumf.tab) must be located in the same directory as this script for the syntax to run


# ----------------------------------------------------------------------------------------



# don't run
# conda install pandas numpy scipy statsmodels matplotlib

import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm



# download the file
#  url = "https://borealisdata.ca/api/access/datafile/661489?format=RData", # define url
DATA_FILE = "gss_34_incident_pumf.tab" # target file name

if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(
        f"Cannot find {DATA_FILE}. Place the tab-delimited GSS file in the same "
        "folder as this Python script, or edit DATA_FILE to specify its location."
    )

# read the data
gss = pd.read_csv(DATA_FILE, sep="\t", low_memory=False) # low_memory=False prevents warnings

# select columns
selected_columns = [
    "PFO_100",
    "CIR_D010",
    "OFFENDRC",
    "WHR_100",
    "TTA_110",
    "TTA_120",
    "TTA_130",
    "TTA_140",
    "TTA_150",
    "TTA_160",
]

missing_columns = [column for column in selected_columns if column not in gss.columns]
if missing_columns:
    raise KeyError(
        "The tab-delimited file does not contain the selected variables: "
        f"{', '.join(missing_columns)}"
    )

# recode as numeric in case pandas read any of the columns as strings
for column in selected_columns:
    gss[column] = pd.to_numeric(gss[column], errors="coerce")

print(f"Loaded {len(gss):,} rows from {DATA_FILE}.")


# recode whether police knew about the incident ----------------------------------------
gss.loc[gss["PFO_100"] > 1, "PFO_100"] = 0 # = 1 if incident known or reported to police, = 0 if not

print("\nPFO_100 distribution after recoding:")
print(gss["PFO_100"].value_counts(dropna=False).sort_index())



# define parameters for estimation -------------------------------------------------------------

# number of police-reported incidents according to Statistics Canada, 2019: https://www150.statcan.gc.ca/n1/pub/85-002-x/2020001/article/00010-eng.htm
police = 2.2e6

# sample size of the General Social Survey incident file
victim = len(gss)

# number of GSS incidents reported as known to police
recapture = int((gss["PFO_100"] == 1).sum())

if recapture == 0:
    raise ValueError("No observations with PFO_100 == 1 were found after recoding.")

# drop missing observations
gss_baseline = gss.dropna(subset=["PFO_100"]).copy()



# 1) the arithmetic to calculate the probability any given incident is known to the police -------------------------------------------------------------

# dual systems estimation a.k.a. Lincoln-Petersen formula
N_hat = (police * victim) / recapture

# Seber's approximate variance formula
# https://books.google.co.uk/books/about/Estimation_of_Animal_Abundance.html?id=iIAAPQAACAAJ&redir_esc=y
variance = (police**2 * victim * (victim - recapture)) / (recapture**3)

# standard error
standard_error = np.sqrt(variance)

# critical value 
z_critical = norm.ppf(0.975) 

# 95% confidence intervals
N_lo = N_hat - z_critical * standard_error
N_hi = N_hat + z_critical * standard_error

# compile results into data frame
results = pd.DataFrame(
    {
        "statistics": [
            "police-reported incidents, 2019",
            "estimated incidents",
            "standard error of the estimate",
            "95% CI lower bound",
            "95% CI upper bound",
            "incidents known to police, point estimate (%)",
            "incidents known to police, lower bound (%)",
            "incidents known to police, upper bound (%)",
        ],
        "estimates": [
            police,
            round(N_hat, 2),
            round(standard_error, 2),
            round(N_lo, 2),
            round(N_hi, 2),
            round((police / N_hat) * 100, 2),
            round((police / N_lo) * 100, 2),
            round((police / N_hi) * 100, 2),
        ],
    }
)

# print results
print("\nDual-systems estimation results:")
print(results.to_string(index=False))



# 2) logistic model to estimate the probability any given incident is known to police ------------------------------------------------------------- 
m_00 = smf.glm(
    formula="PFO_100 ~ 1",
    data=gss_baseline,
    family=sm.families.Binomial(),
).fit()

# intercept
intercept = m_00.params["Intercept"]

# standard error
intercept_se = m_00.bse["Intercept"]

# inverse logit function to calculate baseline probability that an incident is known to police
def inv_logit(x):
    """Convert a log-odds value to its predicted probability."""
    return 1 / (1 + np.exp(-x))

# calculate the baseline probability 
baseline_probability = inv_logit(intercept)
baseline_lo = inv_logit(intercept - z_critical * intercept_se)
baseline_hi = inv_logit(intercept + z_critical * intercept_se)

# print results
print("\nIntercept-only logistic model:")
print(f"Probability incident was known to police: {baseline_probability:.4f}")
print(f"95% CI: [{baseline_lo:.4f}, {baseline_hi:.4f}]")

# note: notice the arithmetic and the logistic regression result in near identical estimates



# 3) add covariate information into the logistic model to increase the precision of the estimate ---------------------------------------------------------
# note: labels for the numeric values for all variables can be found online at: https://borealisdata.ca/data-explorer/?siteUrl=https:%2F%2Fborealisdata.ca&dvLocale=en&fileId=661489

# type of crime/incident ---------------------------------------------

# first, group alike categories

  # attempted theft
  gss.loc[gss["CIR_D010"] == 9, "CIR_D010"] = 2

   # stolen property
  gss.loc[gss["CIR_D010"].isin([5, 6]), "CIR_D010"] = 4

   # sexual assualt and related crimes
  gss.loc[gss["CIR_D010"].isin([13, 14]), "CIR_D010"] = 12

# add labels to numeric values
crime_labels = {
    1: "household damage",
    2: "attempt to take or steal something",
    3: "attempt to break and enter",
    4: "stolen property",
    7: "theft or attempted theft of vehicle",
    8: "damage to vehicles",
    10: "assault",
    11: "threat of assault",
    12: "sexual assault",
    15: "another crime",
}

# specify reference category
crime_order = [
    "stolen property", # first category is the reference category
    "household damage",
    "attempt to take or steal something",
    "attempt to break and enter",
    "theft or attempted theft of vehicle",
    "damage to vehicles",
    "assault",
    "threat of assault",
    "sexual assault",
    "another crime",
]

# map labels to numeric values
gss["CIR_D010"] = gss["CIR_D010"].map(crime_labels)
gss["CIR_D010"] = pd.Categorical(
    gss["CIR_D010"],
    categories=crime_order,
    ordered=False,
)



# recode values for variable asking respondents about their relationship to the offender(s) ---------------------------------------------
gss.loc[gss["OFFENDRC"] == 9, "OFFENDRC"] = 6

# add labels to numeric values
relationship_labels = {
    1: "relative",
    2: "friend/neighbour/acquaintance",
    3: "stranger",
    4: "other type of relation",
    6: "valid skip i.e., no relation",
}

# specify reference category
relationship_order = [
    "valid skip i.e., no relation", # first category is the reference category
    "relative",
    "friend/neighbour/acquaintance",
    "stranger",
    "other type of relation",
]

# map labels to numeric values
gss["OFFENDRC"] = gss["OFFENDRC"].map(relationship_labels)
gss["OFFENDRC"] = pd.Categorical(
    gss["OFFENDRC"],
    categories=relationship_order,
    ordered=False,
)



# location where the incident happened ---------------------------------------------
gss.loc[gss["WHR_100"] == 9, "WHR_100"] = 5

# add labels to numeric values
location_labels = {
    1: "home/nearby",
    2: "other private residence",
    3: "commercial or institutional establishment",
    4: "street or other public place",
    5: "other place",
}

# specify the reference catefory
location_order = [
    "home/nearby", # first category is the reference category
    "other private residence",
    "commercial or institutional establishment",
    "street or other public place",
    "other place",
]

# map labels to numeric values
gss["WHR_100"] = gss["WHR_100"].map(location_labels)
gss["WHR_100"] = pd.Categorical(
    gss["WHR_100"],
    categories=location_order,
    ordered=False,
)



# recode categorical variables that incidate if respondent / victom spoke with anyone about the incident ---------------------------------------------

# recode valid response values above 1 to 0, matching the R script.
talk_variables = [
    "TTA_110",  # spoke with family member
    "TTA_120",  # spoke with friend or neighbour
    "TTA_130",  # spoke with co-worker
    "TTA_140",  # spoke with doctor or nurse
    "TTA_150",  # spoke with lawyer
    "TTA_160",  # spoke with priest
]

for variable in talk_variables:
    gss.loc[gss[variable] > 1, variable] = 0


# logistic model to estimate the probability is known to police , controlling for incident characteristics 
m_01 = smf.glm(
    formula="""
     PFO_100 ~
     C(CIR_D010, Treatment(reference='stolen property')) +
     C(OFFENDRC, Treatment(reference='valid skip i.e., no relation')) +
     C(WHR_100, Treatment(reference='home/nearby')) +
     TTA_110 + 
     TTA_120 + TTA_130 + TTA_140 + TTA_150 + TTA_160
     """,
    data=gss,
    family=sm.families.Binomial(),
).fit()

# print model summary
print(m_01.summary())



# predicted probabilities and confidence intervals ------------------------------------

# for valid observations
estimation_data = gss.loc[m_01.model.data.row_labels].copy()

# predicted probability that each incident is known to police
estimation_data["mu"] = m_01.predict(estimation_data)

# confidence intervals for the predicted probability
prediction_frame = m_01.get_prediction(estimation_data).summary_frame(alpha=0.05)
estimation_data["lo"] = prediction_frame["mean_ci_lower"].to_numpy()
estimation_data["hi"] = prediction_frame["mean_ci_upper"].to_numpy()

# geometric average predicted probability, adjusts for skewness in probability distribution
positive_predictions = estimation_data.loc[
    estimation_data["PFO_100"] == 1,
    "mu",
]
positive_predictions = positive_predictions[positive_predictions > 0] # truncate to prevent log(0) errors

# internal function to calculate geometric mean
p_y1 = np.exp(np.mean(np.log(positive_predictions)))

# print result
print(f"\ngeometric average predicted probability that an incident is known to the police: {p_y1:.4f}")



# function to calculate Tjur's R-squared ---------------------------------------
def geometric_mean(values):
    """calculate a geometric mean after excluding zeros to avoid log(0)."""
    values = np.asarray(values, dtype=float)
    values = values[values > 0]

    if len(values) == 0:
        return np.nan

    return np.exp(np.mean(np.log(values)))


def tjur(model_data):
    """geometric mean of fitted probabilities where y=1 minus the geometric mean where y=0."""
    predicted_y1 = model_data.loc[model_data["PFO_100"] == 1, "mu"]
    predicted_y0 = model_data.loc[model_data["PFO_100"] == 0, "mu"]
    
  # calculate the modified Tjur's R-squared
    r_squared = geometric_mean(predicted_y1) - geometric_mean(predicted_y0)

    print(
        "R-squared ranges 0-1, where higher values indicate models with more predictive power."
    )
    return r_squared

tjur_r2 = tjur(estimation_data)
print(f"Tjur R-squared: {tjur_r2:.4f}")


# plot the predicted probability distribution that incident was known to the police by crime type ----------------------------------------
plot_data = estimation_data.dropna(subset=["CIR_D010", "mu"]).copy()
categories = list(plot_data["CIR_D010"].cat.categories)

# Matplotlib's boxplot accepts one numeric array per category
probability_by_crime = [
    plot_data.loc[plot_data["CIR_D010"] == category, "mu"].to_numpy()
    for category in categories
]

# dimensions
fig, ax = plt.subplots(figsize=(10, 10))

# geometric average predicted probability for incidents known to police
ax.axhline(
    y=p_y1,
    color="#F8766D",
    linestyle="-",
    linewidth=1,
    alpha=0.4,
    zorder=0,
)

# plot the data points
rng = np.random.default_rng(seed=1) # the fixed seed makes the position of jittered observations reproducible ach time the script is run

for position, probabilities in enumerate(probability_by_crime, start=1):
    x_jitter = rng.uniform(
        low=position - 0.20,
        high=position + 0.20,
        size=len(probabilities),
    )

    ax.scatter(
        x_jitter,
        probabilities,
        color="BDBDBD",
        alpha=0.2,
        s=12,
        edgecolors="none",
        zorder=1,
    )

# boxplots drawn over jittered points
boxprops = {
    "facecolor": "white",
    "edgecolor": "#333333",
    "linewidth": 1.2,
}
medianprops = {
    "color": "#333333",
    "linewidth": 1.8,
}
whiskerprops = {
    "color": "#333333",
    "linewidth": 1.0,
}
capprops = {
    "color": "#333333",
    "linewidth": 1.0,
}


# boxplots overlay the jittered points
ax.boxplot(
    probability_by_crime,
    positions=range(1, len(categories) + 1),
    widths=0.55,
    showfliers=True, # retain the outlier points
    patch_artist=True,
    boxprops=boxprops,
    medianprops=medianprops,
    whiskerprops=whiskerprops,
    capprops=capprops,
    flierprops={
        "marker": "o",
        "markerfacecolor": "#333333",
        "markeredgecolor": "#333333",
        "markersize": 3.5,
        "alpha": 0.95,
    },
    zorder=2,
)

# x-axis
ax.set_xlim(0.4, len(categories) + 0.6)
ax.set_xticks(range(1, len(categories) + 1))
ax.set_xticklabels(categories, rotation=45, ha="right", rotation_mode="anchor")
ax.set_xlabel("Type of incident", fontsize=12, color="black")

# y-axis
ax.set_ylim(0, 1)
ax.set_yticks([0.00, 0.25, 0.50, 0.75, 1.00])
ax.set_yticklabels(["0.00", "0.25", "0.50", "0.75", "1.00"])
ax.set_ylabel("Predicted probability", fontsize=12, color="black")

# y-axis label
ax.set_title(
    "The predicted probability that an incident is known to the police by type of crime",
    loc="left",
    fontsize=12,
    fontweight="bold",
    color="black",
    pad=8,
)

# approximate ggplot2::theme_classic().
ax.tick_params(axis="both", labelsize=10, colors="black", width=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#333333")
ax.spines["bottom"].set_color("#333333")
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

caption = (
    "Figure notes: The boxplots present the predicted probability that an incident "
    "is known to the police calculated from self-reports to the 2019 General Social "
    "Survey (GSS) on victimization (N=5,522). The red line indicates the geometric "
    "average predicted probability, which accounts for skewness in the probability "
    "distribution."
)

# wrap caption text
fig.text(
    0.08,
    0.02,
    "\n".join(textwrap.wrap(caption, width=125)),
    ha="left",
    va="bottom",
    fontsize=9,
    color="black",
)

# reserve space at the bottom for the caption
plt.tight_layout(rect=[0, 0.10, 1, 1])

# Output high-resolution figure to the Desktop.
output_path = os.path.expanduser("~/Desktop/fig01.png")
plt.savefig(output_path, dpi=250, facecolor="white", bbox_inches="tight")
plt.close()

print(f"\nFigure saved to: {output_path}")



# descriptive summary of probability distributions by crime type ----------------------

def safe_geometric_mean(values):
    """Geometric mean for grouped summaries, excluding missing and zero values."""
    values = values.dropna()
    values = values[values > 0]

    if len(values) == 0:
        return np.nan

    return np.exp(np.mean(np.log(values)))


descriptives = (
    plot_data.groupby("CIR_D010", observed=True)
    .agg(
        arithmetic_mean=("mu", "mean"),
        geometric_mean=("mu", safe_geometric_mean),
        median=("mu", "median"),
        sample_size=("mu", "size"),
    )
    .reset_index()
)

print("\nDescriptive summary by crime type:")
print(descriptives.to_string(index=False))





# end script


