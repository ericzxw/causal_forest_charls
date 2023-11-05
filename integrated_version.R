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

#df <- read_xlsx("summarydata.xlsx")

#cols <- names(df)
#vars <- c("sex","age","retire","marital_status","residence","live_near_to_child",
#          "education","chronic_disease","sleep_hrs","smoke","alcohol","dep")
#summary(df[,-(1:2)])
#df <- df[df$age>=60 & df$cog_score != -999,]
#df <- na.omit(df)

#df <- as.data.frame(df)

get_formula <-  function(x,covariates){
  Formula <- formula(paste(x, "~", paste(covariates, collapse = " + ")))
  return(Formula)
}


#倾向性得分匹配
perform_propensity_score_matching <- function(data, treatment_var, covariates, caliper = 0.2) {
  # Fit propensity score model
  psModel <- glm(formula = get_formula(treatment_var,covariates),
                 family  = binomial(link = "logit"),
                 data    = data)
  # Calculate propensity scores
  data$p_sp_1 <- predict(psModel, type = "response")
  data$p_sp_0 <- 1 - data$p_sp_1
  # Perform matching
  listMatch <- Match(Tr = (data[, treatment_var] == 1),
                     X = log(data$p_sp_1 / data$p_sp_0),
                     M = 1,
                     caliper = caliper,
                     replace = FALSE,
                     ties = TRUE,
                     version = "fast")
  # Check matching balance
  mb <- MatchBalance(psModel$formula, data = data, match.out = listMatch, nboots = 500)
  
  # Extract matched data
  matched_data <- data[unlist(listMatch[c("index.treated", "index.control")]), ]
  # Return results
  results <- list(
    psModel = psModel,
    listMatch = listMatch,
    mb = mb,
    matched_data = matched_data
    )
  return(results)
}


#matching_results <- perform_propensity_score_matching(data = df, 
#                                                      treatment_var = "sp", 
#                                                      covariates = vars)

#dm <- matching_results$matched_data

get_ci <- function(estimate_mean, standard_deviation) {
  z_value <- 1.96
  lower_limit <- estimate_mean - (z_value * standard_deviation)
  upper_limit <- estimate_mean + (z_value * standard_deviation)
  return(c(lower_limit, upper_limit))
}


#不划分训练集和测试集进行预测
get_cf <- function(data, outcome_var, treatment_var,covariates,min_node_size = 20){
  Y <- data[, outcome_var]
  W <- data[, treatment_var]
  X <- data[, covariates]
  cf <- causal_forest(X,Y,W,min.node.size = 20)
  return(cf)
}
#X <- dm[,vars]

#cf <- get_cf(data = dm, outcome_var = "cog_score", treatment_var = "sp",covariates = vars) 
#ate <- average_treatment_effect(cf,target.sample = "all")
#ate
#varimp <- variable_importance(cf)
#ranked.vars <- order(varimp, decreasing = TRUE)
#colnames(X)[ranked.vars[1:5]] #只选前5个最重要的影响因子
#best_linear_projection(cf,X[ranked.vars[1:5]])

#划分测试集和训练集得到TOC曲线
get_toc <- function(data, outcome_var, treatment_var, covariates, num_trees = 2000, min_node_size = 10, train_index) {
  # Extract Y, W, and X
  Y <- data[, outcome_var]
  W <- data[, treatment_var]
  X <- data[, covariates]
  

#  train_data <- data[train_index, ]
#  test_data <- data[-train_index, ]
  
  # Train a causal forest on the training data
  train_forest <- causal_forest(X = X[train_index, ], Y = Y[train_index], W = W[train_index], 
                                sample.fraction = 0.5, honesty = TRUE, num.trees = num_trees, min.node.size = min_node_size)
  
  # Predict CATE on the test data
  tau_hat_eval <- predict(train_forest, X[-train_index, ])$predictions
  
  # Train a separate forest for evaluation
  eval_forest <- causal_forest(X = X[-train_index, ], Y = Y[-train_index], W = W[-train_index],
                               sample.fraction = 0.5, honesty = TRUE, num.trees = num_trees, min.node.size = min_node_size)
  
  # Calculate the TOC curve
  rate_cate <- rank_average_treatment_effect(eval_forest, tau_hat_eval)
  
  
  # Return the trained forests if needed
  results <- list(
    Train_Forest = train_forest,
    Evaluation_Forest = eval_forest,
    TOC_Curve = rate_cate
  )
  
  return(results)
}


#results <- get_toc(data = dm, outcome_var = "cog_score", treatment_var = "sp",
#                                  covariates = vars,
#                                  sample_fraction = 0.5)

#train <- results$Train




get_policy_tree <- function(cf,ate, X, train, n, depth = 2, min_node_size = 100) {
  # Compute doubly robust scores
  dr.scores <- get_scores(cf)
  
  # Extract the ATE as a "cost" of program treatment
  cost <- ate[["estimate"]]
  
  # Compute rewards for control and treat groups
  dr.rewards <- data.frame(
    control = -dr.scores,
    treat = dr.scores - cost
  )
  
  # Calculate variable importance
  varimp <- variable_importance(cf)
  
  # Order variables by importance
  ranked_vars <- order(varimp, decreasing = TRUE)
  
  # Take the top 'n' variables from the ranked list
  selected_vars <- ranked_vars[1:n]
  
  # Fit a policy tree using the selected variables
  tree <- policy_tree(X[train, selected_vars], dr.rewards[train, ], depth = depth, min.node.size = min_node_size)
  
  return(tree)
}



#tree <- get_policy_tree(cf,ate,X,train,5)
#plot(tree,leaf.labels = c("no treat", "treat"))



# Compute doubly robust scores
#dr.scores <- get_scores(cf) #使用get scores得到的其实是Y1-Y0的结果，这里Y1和Y0都是doubly robust estimate,即样本实施政策与否的预测结果
# Use as the ATE as a "cost" of program treatment to find something non-trivial
#cost <- ate[["estimate"]] #把ate作为一个cost,这里ate是一个面向所有样本估计的ate
#dr.rewards<- cbind(control = -dr.scores, treat=dr.scores - cost) 
#此时dr.rewards即为：每一个样本，如果把它放在control组得到的回报值，和如果把它treat组得到的回报值
#如果control组，回报即为-dr.score = Y0-Y1，即如果不实施这个干预政策，这个样本可以得到一些reward
#如果treat组，回报为dr.scores -cost = Y1-Y0-cost,即如果实施了这个干预政策，这个样本得到reward，这个reward里去掉了ate，即如果每个人都有ate的效果，那么我们这些实施了干预政策的样本能够得到的额外的reward

# Fit depth 2 tree on training subset
#tree <- policy_tree(X[train,ranked.vars[1:5]], dr.rewards[train, ], depth = 2, min.node.size = 100)
#plot(tree, leaf.labels = c("no treat", "treat"))

#检验policy tree的效果
# Get the leaf node assigned to each test sample.找到每一个测试集样本的叶子节点
#

#node.id <- predict(tree, X[-train,ranked.vars[1:5]], type = "node.id")

# Doubly robust estimates of E[Y(control)] and E[Y(treated)] by leaf node.
#values <- aggregate(dr.scores[-train,], by = list(leaf.node = node.id),
#                    FUN = function(dr) c(mean = mean(dr), se = sd(dr) / sqrt(length(dr))))
#print(values, digits = 1)



