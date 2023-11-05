setwd("C:/Users/Xinwen Zou/Documents/causal_forest_charls")
set.seed(123)
library(tableone)
## PS matching
library(Matching)
## Weighted analysis
library(survey)
library(reshape2)
library(ggplot2)
library(grf)
library(readxl)
library(policytree)
options(scipen = 999)
source("integrated_version.R")

df <- read_xlsx("summarydata.xlsx")

cols <- names(df)
vars <- c("sex","age","retire","marital_status","residence","live_near_to_child",
          "education","chronic_disease","sleep_hrs","smoke","alcohol","dep")
summary(df[,-(1:2)])
df <- df[df$age>=60 & df$cog_score!= -999 & df$orient != -999 & df$word_recall!=-999 &df$draw != -999 & df$serial != -999 & df$sp != -999,]
df <- na.omit(df)
df <- as.data.frame(df)

#perform psm and get the final matched data
matching_results <- perform_propensity_score_matching(data = df, 
                                                      treatment_var = "sp", 
                                                      covariates = vars)
dm <- matching_results$matched_data
outcome_variables <- c("orient", "word_recall", "draw", "serial", "cog_score")
X <- dm[,vars]


# using stats models for examining relationship
outcome_variables <- c("orient", "word_recall", "draw", "serial", "cog_score")




# Loop through outcome variables
for (outcome_variable in outcome_variables) {
  if (outcome_variable == "draw") {
    # Perform logistic regression for the binary outcome variable "draw"
    model <- glm(paste("draw", "~ sp + ", paste(vars, collapse = " + ")), family = binomial, data = dm)
  } else {
    # Perform linear regression for other outcome variables
    model <- lm(paste(outcome_variable, "~ sp + ", paste(vars, collapse = " + ")), data = dm)
  }
  cat("results of:", outcome_variable,"\n")
  print(summary(model))

}








forests <- list()
ates <- list()

for (outcome_variable in outcome_variables) {
  # Get the causal forest for the outcome variable
  cf <- get_cf(data = dm, outcome_var =  outcome_variable,treatment_var = "sp",covariates = vars)
  
  # Store the forest in the list
  forests[[outcome_variable]] <- cf
  
  # Compute the Average Treatment Effect (ATE)
  ates[[outcome_variable]]<- average_treatment_effect(cf)
  
  # Calculate the confidence interval for the ATE
  ci <- get_ci(ates[[outcome_variable]][[1]],ates[[outcome_variable]][[2]])
  
  # Print the results for the current outcome variable
  cat("Outcome Variable:", outcome_variable, "\n")
  cat("ATE:", ates[[outcome_variable]][[1]], "\n")
  cat("ATE std.err:",ates[[outcome_variable]][[2]],"\n")
  cat(ci, "\n")
  cat("\n")
  
}


# Split data into training and testing sets
n <- nrow(dm)
train_index <- sample(1:n, size = 0.5 * n)  # User-defined sample fraction


# Plot TOC for each outcome var
for (outcome_variable in outcome_variables) {
  # Calculate the TOC curve
  toc <- get_toc(data = dm, outcome_var = outcome_variable,treatment_var = "sp",covariates = vars, train_index=train_index)
  # Plot the TOC curve
  toc_plot <- plot(toc$TOC_Curve, main = paste("TOC for", outcome_variable))
}

trees <- list()

for (outcome_variable in outcome_variables){
  trees[[outcome_variable]] <- get_policy_tree(forests[[outcome_variable]],ates[[outcome_variable]],X,train_index,n=5)
  print(plot(trees[[outcome_variable]],leaf.labels = c("no treat", "treat")))
}

