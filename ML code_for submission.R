
# This software is licensed under the MIT License.
# See the license.txt file for more information.

# Read parameters
param_df <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- setNames(as.list(as.numeric(param_df$value)), param_df$parameter)
ran_seed <- params$ran_seed   # 42 was on our works
cutoff <- params$cutoff   # 0 : no cutoff, 1 : cutoff at certain concentration of CA19-9
kind_of_system <- params$kind_of_system   # 0 : ELISA, 1 : AIQ-LFA
n_of_boot <- params$n_of_boot   # number of bootstraps
train_ratio <- params$train_ratio    # ratio for train set
FPR_rate <- 1 - (params$FPR_rate)    # FPR rate for set the cutoff, changed to specificity (1-FPR)
TPR_rate <- params$TPR_rate    # TPR rate for set the cutoff


# set the random seed
set.seed(ran_seed)

# loading the packages
library(readxl)
library(pROC)
library(openxlsx)
library(caret)
library(dplyr)

# Raw data reading
tdata <- read_excel("raw datas.xlsx", sheet = "Sheet1", col_names = TRUE)

# analysis mode settings

# p value variables
p_value_valid_SP_HE <- numeric(n_of_boot)
p_value_valid_SP_HL <- numeric(n_of_boot)
p_value_valid_SP_whole <- numeric(n_of_boot)
p_value_valid_Pc_HE <- numeric(n_of_boot)
p_value_valid_Pc_HL <- numeric(n_of_boot)
p_value_valid_Pc_whole <- numeric(n_of_boot)

dig_para <- data.frame(
  cutoff = c("nocutoff", "cutoff"),
  system = c("AI-ELISA", "AIQ-LFA")
)

folder_name_addon <- paste0(dig_para$cutoff[cutoff+1],"_",dig_para$system[kind_of_system+1])

# cutoff setting
if (cutoff == 1) {
  
  # cutoff value
  cutoff_ELISA <- 37
  cutoff_LFA <- 3.09
  
  tdata <- tdata %>%
    mutate(
      EC = ifelse(EC > cutoff_ELISA, 1, 0),
      QC = ifelse(QC > cutoff_LFA, 1, 0)
    )
}

# formular for model
temp_name <- c("EC", "QC")    # EC: ELISA CA19-9, QC: (Q)LFA CA19-9
formula <- reformulate(response = "stats", termlabels = temp_name[kind_of_system+1])
  
# size of each sample group
n_of_samples_H <- 50
n_of_samples_E <- 40
n_of_samples_L <- 60


# vector for stacking the AUC value
AUC_average_ROC_HE <- numeric(n_of_boot)
AUC_average_ROC_HL <- numeric(n_of_boot)
AUC_average_ROC_whole <- numeric(n_of_boot)

# modifying the data set and assign the train/test data sets
tdata$stats <- factor(tdata$stats, levels = c(0,1), labels = c("Class0", "Class1"))   # stats 0: healthy, 1: PDAC patient
tdata <- cbind(sample_no = seq(1:150), tdata)   # Sample no.1 to 50: healthy, no.51 to 90: early stage, no.91 to 150: Late stage

df_H <- tdata[1:n_of_samples_H,]
df_E <- tdata[(n_of_samples_H+1):(n_of_samples_H+n_of_samples_E),]
df_L <- tdata[(n_of_samples_H+n_of_samples_E+1):(n_of_samples_H+n_of_samples_E+n_of_samples_L),]

train_idx_H <- sample(1:n_of_samples_H, (n_of_samples_H*train_ratio))
train_idx_E <- sample(1:n_of_samples_E, (n_of_samples_E*train_ratio))
train_idx_L <- sample(1:n_of_samples_L, (n_of_samples_L*train_ratio))

df_H_train <- df_H[train_idx_H,]
df_E_train <- df_E[train_idx_E,]
df_L_train <- df_L[train_idx_L,]

df_HE_train_ori <- rbind (df_H_train, df_E_train)
df_HL_train_ori <- rbind (df_H_train, df_L_train)
df_whole_train_ori <- rbind (df_H_train, df_E_train, df_L_train)

df_H_test <- df_H[-train_idx_H,]
df_E_test <- df_E[-train_idx_E,]
df_L_test <- df_L[-train_idx_L,]

df_HE_test <- rbind (df_H_test, df_E_test)
df_HL_test <- rbind (df_H_test, df_L_test)
df_whole_test <- rbind (df_H_test, df_E_test, df_L_test)

# matrix for store the possibility results
p_HE <- matrix (NA, nrow = n_of_samples_H+n_of_samples_E, ncol = n_of_boot)
p_HL <- matrix (NA, nrow = 150, ncol = n_of_boot)
p_whole <- matrix (NA, nrow = n_of_samples_H+n_of_samples_E+n_of_samples_L, ncol = n_of_boot)

# matrix for store the count for pick up while bootstrapping
count_HE <- matrix (NA, nrow = n_of_samples_H+n_of_samples_E, ncol = n_of_boot)
count_HL <- matrix (NA, nrow = n_of_samples_H+n_of_samples_E+n_of_samples_L, ncol = n_of_boot)
count_whole <- matrix (NA, nrow = n_of_samples_H+n_of_samples_E+n_of_samples_L, ncol = n_of_boot)

# variable for AUC stacking
AUC_train_SP <- matrix (NA, nrow = n_of_boot, ncol = 3)
AUC_test_SP <- matrix (NA, nrow = n_of_boot, ncol = 3)
AUC_train_Pc <- matrix (NA, nrow = n_of_boot, ncol = 3)
AUC_test_Pc <- matrix (NA, nrow = n_of_boot, ncol = 3)


# Bootstrapping
for (i in 1:n_of_boot) {
  
  idx_H_train <- sample (1:(n_of_samples_H*train_ratio), replace = TRUE)
  idx_E_train <- sample (1:(n_of_samples_E*train_ratio), replace = TRUE)
  idx_L_train <- sample (1:(n_of_samples_L*train_ratio), replace = TRUE)
  
  df_HE_train <- rbind(df_H_train[idx_H_train,], df_E_train[idx_E_train,])
  df_HL_train <- rbind(df_H_train[idx_H_train,], df_L_train[idx_L_train,])
  df_whole_train <- rbind(df_H_train[idx_H_train,], df_E_train[idx_E_train,], df_L_train[idx_L_train,])
  
  # modeling
  
  model_HE <- glm (formula, data = df_HE_train, family = binomial)
  model_HL <- glm (formula, data = df_HL_train, family = binomial)
  model_whole <- glm (formula, data = df_whole_train, family = binomial)

  sample_HE_train <- unique(df_HE_train$sample_no)
  sample_HL_train <- unique(df_HL_train$sample_no)
  sample_whole_train <- unique(df_whole_train$sample_no)
  
  sample_HE_test <- unique(df_HE_test$sample_no)
  sample_HL_test <- unique(df_HL_test$sample_no)
  sample_whole_test <- unique(df_whole_test$sample_no)
  
  p_HE[df_HE_train$sample_no,i] <- predict(model_HE, type="response")
  p_HL[df_HL_train$sample_no,i] <- predict(model_HL, type="response")
  p_whole[df_whole_train$sample_no,i] <- predict(model_whole, type="response")
  
  p_HE[df_HE_test$sample_no,i] <- predict(model_HE, newdata = df_HE_test, type = "response")
  p_HL[df_HL_test$sample_no,i] <- predict(model_HL, newdata = df_HL_test, type = "response")
  p_whole[df_whole_test$sample_no,i] <- predict(model_whole, newdata = df_whole_test, type = "response")
  
  tab_H <- table(idx_H_train)
  tab_E <- table(idx_E_train)
  tab_L <- table(idx_L_train)
  
  count_HE[df_H_train$sample_no[as.numeric(names(tab_H))],i] <- as.numeric(tab_H)
  count_HE[df_E_train$sample_no[as.numeric(names(tab_E))],i] <- as.numeric(tab_E)
  count_HL[df_H_train$sample_no[as.numeric(names(tab_H))],i] <- as.numeric(tab_H)
  count_HL[df_L_train$sample_no[as.numeric(names(tab_L))],i] <- as.numeric(tab_L)
  count_whole[df_H_train$sample_no[as.numeric(names(tab_H))],i] <- as.numeric(tab_H)
  count_whole[df_E_train$sample_no[as.numeric(names(tab_E))],i] <- as.numeric(tab_E)
  count_whole[df_L_train$sample_no[as.numeric(names(tab_L))],i] <- as.numeric(tab_L)
  
  count_HE[df_HE_test$sample_no,i] <- 1
  count_HL[df_HL_test$sample_no,i] <- 1
  count_whole[df_whole_test$sample_no,i] <- 1

  
}


weights <- rowSums(count_HE, na.rm = TRUE)
tdata_HE <- cbind (rbind(df_H, df_E), Avg = rowSums(p_HE*count_HE, na.rm = TRUE)/weights)

p_HL <- p_HL[-(51:90),] # due to HL case, empty row (Early stage, sample 51 to 90) existed.
count_HL <- count_HL[-(51:90),]
weights <- rowSums(count_HL, na.rm = TRUE)
tdata_HL <- cbind (rbind(df_H, df_L), Avg = rowSums(p_HL*count_HL, na.rm = TRUE)/weights)

weights <- rowSums(count_whole, na.rm = TRUE)
tdata_whole <- cbind (rbind(df_H, df_E, df_L), Avg = rowSums(p_whole * count_whole, na.rm = TRUE)/weights)


HE_idx <- c(train_idx_H, train_idx_E+n_of_samples_H)
HL_idx <- c(train_idx_H, train_idx_L+n_of_samples_H)
whole_idx <- c(train_idx_H, train_idx_E+n_of_samples_H, train_idx_L+n_of_samples_H+n_of_samples_E)


# modeling of weighted average logit function

new_formula <- reformulate(response = "Avg", termlabels = temp_name[kind_of_system+1])

model_P_conc_HE <- glm (new_formula, data = tdata_HE[HE_idx,], family = binomial)
model_P_conc_HL <- glm (new_formula, data = tdata_HL[HL_idx,], family = binomial)
model_P_conc_whole <- glm (new_formula, data = tdata_whole[whole_idx,], family = binomial)

model_stats_P_HE <- glm (stats ~ Avg, data = tdata_HE[HE_idx,], family = binomial)
model_stats_P_HL <- glm (stats ~ Avg, data = tdata_HL[HL_idx,], family = binomial)
model_stats_P_whole <- glm (stats ~ Avg, data = tdata_whole[whole_idx,], family = binomial)

# Representative ROC curve

avg_ROC_train_HE <- roc (tdata_HE$stats[HE_idx], predict(model_stats_P_HE, type = "response"))
avg_ROC_train_HL <- roc (tdata_HL$stats[HL_idx], predict(model_stats_P_HL, type = "response"))
avg_ROC_train_whole <- roc (tdata_whole$stats[whole_idx], predict(model_stats_P_whole, type = "response"))

avg_ROC_test_HE <- roc (tdata_HE$stats[-HE_idx], predict(model_stats_P_HE, newdata = tdata_HE[-HE_idx,], type = "response"))
avg_ROC_test_HL <- roc (tdata_HL$stats[-HL_idx], predict(model_stats_P_HL, newdata = tdata_HL[-HL_idx,], type = "response"))
avg_ROC_test_whole <- roc (tdata_whole$stats[-whole_idx], predict(model_stats_P_whole, newdata = tdata_whole[-whole_idx,], type = "response"))


# bootstrapping for validation

FPR_sen <- matrix (NA, nrow = n_of_boot, ncol = 3)
TPR_spe <- matrix (NA, nrow = n_of_boot, ncol = 3)
cutoff_tpr <- matrix (NA, nrow = n_of_boot, ncol = 3)
cutoff_fpr <- matrix (NA, nrow = n_of_boot, ncol = 3)

# making temporary dataframe for HL

blank_rows <- tdata_HL[0,]
blank_rows <- blank_rows[rep(1, 40), ]
blank_rows[,] <- NA
temp_HL_test <- rbind(tdata_HL[1:50,],blank_rows,tdata_HL[51:110,])


set.seed(ran_seed)

for (i in 1:n_of_boot) {
  
  # resampling with bootstrapping
  
  valid_idx_HE <- sample(HE_idx, replace = TRUE)
  valid_idx_HL <- sample(HL_idx, replace = TRUE)
  valid_idx_whole <- sample(whole_idx, replace = TRUE)
  
  valid_test_HE <- sample(nrow(df_HE_test), replace = TRUE)
  valid_test_HL <- sample(nrow(df_HL_test), replace = TRUE)
  valid_test_whole <- sample(nrow(df_whole_test), replace = TRUE)
  
  
  # validation of stats <- P model (model_stats_P)
  
  temp_roc_train_SP_HE <- roc(tdata_HE$stats[valid_idx_HE], predict(model_stats_P_HE, newdata = tdata_HE[valid_idx_HE,], type = "response"))
  temp_roc_train_SP_HL <- roc(tdata_HL$stats[valid_idx_HL], predict(model_stats_P_HL, newdata = tdata_HL[valid_idx_HL,], type = "response"))
  temp_roc_train_SP_whole <- roc(tdata_whole$stats[valid_idx_whole], predict(model_stats_P_whole, newdata = tdata_whole[valid_idx_whole,], type = "response"))
  
  temp_roc_test_SP_HE <- roc(tdata_HE$stats[df_HE_test$sample_no[valid_test_HE]], predict(model_stats_P_HE, newdata = tdata_HE[df_HE_test$sample_no[valid_test_HE],], type = "response"))
  temp_roc_test_SP_HL <- roc(temp_HL_test$stats[df_HL_test$sample_no[valid_test_HL]], predict(model_stats_P_HL, newdata = temp_HL_test[df_HL_test$sample_no[valid_test_HL],], type = "response"))
  temp_roc_test_SP_whole <- roc(tdata_whole$stats[df_whole_test$sample_no[valid_test_whole]], predict(model_stats_P_whole, newdata = tdata_whole[df_whole_test$sample_no[valid_test_whole],], type = "response"))
  
  temp_test_SP_HE <- roc.test(temp_roc_train_SP_HE, temp_roc_test_SP_HE, method = "delong")
  temp_test_SP_HL <- roc.test(temp_roc_train_SP_HL, temp_roc_test_SP_HL, method = "delong")
  temp_test_SP_whole <- roc.test(temp_roc_train_SP_whole, temp_roc_test_SP_whole, method = "delong")
  
  p_value_valid_SP_HE[i] <- temp_test_SP_HE$p.value
  p_value_valid_SP_HL[i] <- temp_test_SP_HL$p.value
  p_value_valid_SP_whole[i] <- temp_test_SP_whole$p.value
  
  AUC_train_SP[i,1] <- auc(temp_roc_train_SP_HE)
  AUC_train_SP[i,2] <- auc(temp_roc_train_SP_HL)
  AUC_train_SP[i,3] <- auc(temp_roc_train_SP_whole)
  
  AUC_test_SP[i,1] <- auc(temp_roc_test_SP_HE)
  AUC_test_SP[i,2] <- auc(temp_roc_test_SP_HL)
  AUC_test_SP[i,3] <- auc(temp_roc_test_SP_whole)
  
  # validation of P <- conc model (model_P_conc)
  
  temp_roc_train_Pc_HE <- roc(tdata_HE$stats[valid_idx_HE], predict(model_P_conc_HE, newdata = tdata_HE[valid_idx_HE,], type = "response"))
  temp_roc_train_Pc_HL <- roc(tdata_HL$stats[valid_idx_HL], predict(model_P_conc_HL, newdata = tdata_HL[valid_idx_HL,], type = "response"))
  temp_roc_train_Pc_whole <- roc(tdata_whole$stats[valid_idx_whole], predict(model_P_conc_whole, newdata = tdata_whole[valid_idx_whole,], type = "response"))
  
  temp_roc_test_Pc_HE <- roc(tdata_HE$stats[df_HE_test$sample_no[valid_test_HE]], predict(model_P_conc_HE, newdata = tdata_HE[df_HE_test$sample_no[valid_test_HE],], type = "response"))
  temp_roc_test_Pc_HL <- roc(temp_HL_test$stats[df_HL_test$sample_no[valid_test_HL]], predict(model_P_conc_HL, newdata = temp_HL_test[df_HL_test$sample_no[valid_test_HL],], type = "response"))
  temp_roc_test_Pc_whole <- roc(tdata_whole$stats[df_whole_test$sample_no[valid_test_whole]], predict(model_P_conc_whole, newdata = tdata_whole[df_whole_test$sample_no[valid_test_whole],], type = "response"))
  
  temp_test_Pc_HE <- roc.test(temp_roc_train_Pc_HE, temp_roc_test_Pc_HE, method = "delong")
  temp_test_Pc_HL <- roc.test(temp_roc_train_Pc_HL, temp_roc_test_Pc_HL, method = "delong")
  temp_test_Pc_whole <- roc.test(temp_roc_train_Pc_whole, temp_roc_test_Pc_whole, method = "delong")
  
  p_value_valid_Pc_HE[i] <- temp_test_Pc_HE$p.value
  p_value_valid_Pc_HL[i] <- temp_test_Pc_HL$p.value
  p_value_valid_Pc_whole[i] <- temp_test_Pc_whole$p.value
  
  AUC_train_Pc[i,1] <- auc(temp_roc_train_Pc_HE)
  AUC_train_Pc[i,2] <- auc(temp_roc_train_Pc_HL)
  AUC_train_Pc[i,3] <- auc(temp_roc_train_Pc_whole)
  
  AUC_test_Pc[i,1] <- auc(temp_roc_test_Pc_HE)
  AUC_test_Pc[i,2] <- auc(temp_roc_test_Pc_HL)
  AUC_test_Pc[i,3] <- auc(temp_roc_test_Pc_whole)
  
  # obtaining the cutoff for FPR < 1% and TPR > 99%
  
  coords_fpr_HE <- coords(temp_roc_train_SP_HE, x = "all", ret = c("threshold", "specificity", "sensitivity"))
  coords_fpr_HL <- coords(temp_roc_train_SP_HL, x = "all", ret = c("threshold", "specificity", "sensitivity"))
  coords_fpr_whole <- coords(temp_roc_train_SP_whole, x = "all", ret = c("threshold", "specificity", "sensitivity"))
  
  cutoff_fpr_HE <- coords_fpr_HE[coords_fpr_HE$specificity > FPR_rate, "threshold"]
  cutoff_fpr_HL <- coords_fpr_HL[coords_fpr_HL$specificity > FPR_rate, "threshold"]
  cutoff_fpr_whole <- coords_fpr_whole[coords_fpr_whole$specificity > FPR_rate, "threshold"]
  
  cutoff_fpr_HE <- min(cutoff_fpr_HE)
  cutoff_fpr_HL <- min(cutoff_fpr_HL)
  cutoff_fpr_whole <- min(cutoff_fpr_whole)
  
  cutoff_tpr_HE <- coords_fpr_HE[coords_fpr_HE$sensitivity > TPR_rate, "threshold"]
  cutoff_tpr_HL <- coords_fpr_HL[coords_fpr_HL$sensitivity > TPR_rate, "threshold"]
  cutoff_tpr_whole <- coords_fpr_whole[coords_fpr_whole$sensitivity > TPR_rate, "threshold"]
  
  cutoff_tpr_HE <- max(cutoff_tpr_HE)
  cutoff_tpr_HL <- max(cutoff_tpr_HL)
  cutoff_tpr_whole <- max(cutoff_tpr_whole)
  
  temp_pred_HE <- tdata_HE$Avg[valid_idx_HE]
  temp_pred_HL <- tdata_HL$Avg[valid_idx_HL]
  temp_pred_whole <- tdata_whole$Avg[valid_idx_whole]
  
  temp_pred_fpr_HE <- ifelse(temp_pred_HE >= cutoff_fpr_HE, "Class1", "Class0")
  temp_pred_fpr_HL <- ifelse(temp_pred_HL >= cutoff_fpr_HL, "Class1", "Class0")
  temp_pred_fpr_whole <- ifelse(temp_pred_whole >= cutoff_fpr_whole, "Class1", "Class0")
  
  temp_pred_tpr_HE <- ifelse(temp_pred_HE >= cutoff_tpr_HE, "Class1", "Class0")
  temp_pred_tpr_HL <- ifelse(temp_pred_HL >= cutoff_tpr_HL, "Class1", "Class0")
  temp_pred_tpr_whole <- ifelse(temp_pred_whole >= cutoff_tpr_whole, "Class1", "Class0")
  
  
  FPR_sen[i,1] <- coords_fpr_HE$sensitivity[which(coords_fpr_HE$threshold == cutoff_fpr_HE)]
  FPR_sen[i,2] <- coords_fpr_HL$sensitivity[which(coords_fpr_HL$threshold == cutoff_fpr_HL)]
  FPR_sen[i,3] <- coords_fpr_whole$sensitivity[which(coords_fpr_whole$threshold == cutoff_fpr_whole)]
  
  TPR_spe[i,1] <- coords_fpr_HE$specificity[which(coords_fpr_HE$threshold == cutoff_tpr_HE)]
  TPR_spe[i,2] <- coords_fpr_HL$specificity[which(coords_fpr_HL$threshold == cutoff_tpr_HL)]
  TPR_spe[i,3] <- coords_fpr_whole$specificity[which(coords_fpr_whole$threshold == cutoff_tpr_whole)]
  
  cutoff_tpr[i,1] <- cutoff_tpr_HE
  cutoff_tpr[i,2] <- cutoff_tpr_HL
  cutoff_tpr[i,3] <- cutoff_tpr_whole
  
  cutoff_fpr[i,1] <- cutoff_fpr_HE
  cutoff_fpr[i,2] <- cutoff_fpr_HL
  cutoff_fpr[i,3] <- cutoff_fpr_whole
  
}

# Confusion matrix process (for validation)

test_HE <- tdata_HE[-HE_idx,]
test_HL <- tdata_HL[-HL_idx,]
test_whole <- tdata_whole[-whole_idx,]

pred_sample_HE <- predict(model_P_conc_HE, newdata = test_HE, type = "response")
pred_sample_HL <- predict(model_P_conc_HL, newdata = test_HL, type = "response")
pred_sample_whole <- predict(model_P_conc_whole, newdata = test_whole, type = "response")

prob_sample_tpr_HE <- ifelse(pred_sample_HE >= mean(na.omit(cutoff_tpr[,1])), "Class1", "Class0")
prob_sample_tpr_HL <- ifelse(pred_sample_HL >= mean(na.omit(cutoff_tpr[,2])), "Class1", "Class0")
prob_sample_tpr_whole <- ifelse(pred_sample_whole >= mean(na.omit(cutoff_tpr[,3])), "Class1", "Class0")

prob_sample_fpr_HE <- ifelse(pred_sample_HE >= mean(na.omit(cutoff_fpr[,1])), "Class1", "Class0")
prob_sample_fpr_HL <- ifelse(pred_sample_HL >= mean(na.omit(cutoff_fpr[,2])), "Class1", "Class0")
prob_sample_fpr_whole <- ifelse(pred_sample_whole >= mean(na.omit(cutoff_fpr[,3])), "Class1", "Class0")

#conf_matrix_fpr_HE <- confusionMatrix(factor(temp_pred_fpr_HE), factor(tdata_HE$stats[valid_idx_HE]), positive = "Class1")
#conf_matrix_fpr_HL <- confusionMatrix(factor(temp_pred_fpr_HL), factor(tdata_HL$stats[valid_idx_HL]), positive = "Class1")
#conf_matrix_fpr_whole <- confusionMatrix(factor(temp_pred_fpr_whole), factor(tdata_whole$stats[valid_idx_whole]), positive = "Class1")

conf_matrix_tpr_HE <- confusionMatrix(factor(prob_sample_tpr_HE), factor(test_HE$stats), positive = "Class1")
conf_matrix_tpr_HL <- confusionMatrix(factor(prob_sample_tpr_HL), factor(test_HL$stats), positive = "Class1")
conf_matrix_tpr_whole <- confusionMatrix(factor(prob_sample_tpr_whole), factor(test_whole$stats), positive = "Class1")

conf_matrix_fpr_HE <- confusionMatrix(factor(prob_sample_fpr_HE), factor(test_HE$stats), positive = "Class1")
conf_matrix_fpr_HL <- confusionMatrix(factor(prob_sample_fpr_HL), factor(test_HL$stats), positive = "Class1")
conf_matrix_fpr_whole <- confusionMatrix(factor(prob_sample_fpr_whole), factor(test_whole$stats), positive = "Class1")


# xlsx files generation

# dataframe for Bland-Altman plot

df_BA_HE <- data.frame(
  mean_auc = (AUC_train_SP[,1] + AUC_test_SP[,1]) / 2,
  diff_auc = AUC_train_SP[,1] - AUC_test_SP[,1],
  mean_diff = mean(AUC_train_SP[,1] - AUC_test_SP[,1]),
  loa_upper = mean(AUC_train_SP[,1] - AUC_test_SP[,1]) + 1.96 * sd(AUC_train_SP[,1] - AUC_test_SP[,1]),
  loa_lower = mean(AUC_train_SP[,1] - AUC_test_SP[,1]) - 1.96 * sd(AUC_train_SP[,1] - AUC_test_SP[,1])
)

df_BA_HL <- data.frame(
  mean_auc = (AUC_train_SP[,2] + AUC_test_SP[,2]) / 2,
  diff_auc = AUC_train_SP[,2] - AUC_test_SP[,2],
  mean_diff = mean(AUC_train_SP[,2] - AUC_test_SP[,2]),
  loa_upper = mean(AUC_train_SP[,2] - AUC_test_SP[,2]) + 1.96 * sd(AUC_train_SP[,2] - AUC_test_SP[,2]),
  loa_lower = mean(AUC_train_SP[,2] - AUC_test_SP[,2]) - 1.96 * sd(AUC_train_SP[,2] - AUC_test_SP[,2])
)

df_BA_whole <- data.frame(
  mean_auc = (AUC_train_SP[,3] + AUC_test_SP[,3]) / 2,
  diff_auc = AUC_train_SP[,3] - AUC_test_SP[,3],
  mean_diff = mean(AUC_train_SP[,3] - AUC_test_SP[,3]),
  loa_upper = mean(AUC_train_SP[,3] - AUC_test_SP[,3]) + 1.96 * sd(AUC_train_SP[,3] - AUC_test_SP[,3]),
  loa_lower = mean(AUC_train_SP[,3] - AUC_test_SP[,3]) - 1.96 * sd(AUC_train_SP[,3] - AUC_test_SP[,3])
)

# dataframe for DeLong

df_Delong_HE <- data.frame(
  AUC_train = mean(AUC_train_SP[,1]),
  AUC_test = mean(AUC_test_SP[,1]),
  p_value = mean(na.omit(p_value_valid_SP_HE))
)

df_Delong_HL <- data.frame(
  AUC_train = mean(AUC_train_SP[,2]),
  AUC_test = mean(AUC_test_SP[,2]),
  p_value = mean(na.omit(p_value_valid_SP_HL))
)

df_Delong_whole <- data.frame(
  AUC_train = mean(AUC_train_SP[,3]),
  AUC_test = mean(AUC_test_SP[,3]),
  p_value = mean(na.omit(p_value_valid_SP_whole))
)

df_Delong <- rbind(df_Delong_HE, df_Delong_HL, df_Delong_whole)
rownames(df_Delong) <- c("HE", "HL", "whole")

df_Pvalues <- data.frame(
  HE = p_value_valid_SP_HE,
  HL = p_value_valid_SP_HL,
  whole = p_value_valid_SP_whole
)

# just for test (for me)

df_Delong_Pc_HE <- data.frame(
  AUC_train = mean(AUC_train_Pc[,1]),
  AUC_test = mean(AUC_test_Pc[,1]),
  p_value = mean(na.omit(p_value_valid_Pc_HE))
)

df_Delong_Pc_HL <- data.frame(
  AUC_train = mean(AUC_train_Pc[,2]),
  AUC_test = mean(AUC_test_Pc[,2]),
  p_value = mean(na.omit(p_value_valid_Pc_HL))
)

df_Delong_Pc_whole <- data.frame(
  AUC_train = mean(AUC_train_Pc[,3]),
  AUC_test = mean(AUC_test_Pc[,3]),
  p_value = mean(na.omit(p_value_valid_Pc_whole))
)

df_Delong_Pc <- rbind(df_Delong_Pc_HE, df_Delong_Pc_HL, df_Delong_Pc_whole)
rownames(df_Delong_Pc) <- c("HE", "HL", "whole")


# dataframe for ROC cuve

df_roc_train_HE <- data.frame(
  thresholds = avg_ROC_train_HE$thresholds,
  Specificity = avg_ROC_train_HE$specificities,
  FPR = 1 - avg_ROC_train_HE$specificities,
  Sensitivity = avg_ROC_train_HE$sensitivities
)

df_roc_test_HE <- data.frame(
  thresholds = avg_ROC_test_HE$thresholds,
  Specificity = avg_ROC_test_HE$specificities,
  FPR = 1 - avg_ROC_test_HE$specificities,
  Sensitivity = avg_ROC_test_HE$sensitivities
)

df_roc_train_HL <- data.frame(
  thresholds = avg_ROC_train_HL$thresholds,
  Specificity = avg_ROC_train_HL$specificities,
  FPR = 1 - avg_ROC_train_HL$specificities,
  Sensitivity = avg_ROC_train_HL$sensitivities
)

df_roc_test_HL <- data.frame(
  thresholds = avg_ROC_test_HL$thresholds,
  Specificity = avg_ROC_test_HL$specificities,
  FPR = 1 - avg_ROC_test_HL$specificities,
  Sensitivity = avg_ROC_test_HL$sensitivities
)

df_roc_train_whole <- data.frame(
  thresholds = avg_ROC_train_whole$thresholds,
  Specificity = avg_ROC_train_whole$specificities,
  FPR = 1 - avg_ROC_train_whole$specificities,
  Sensitivity = avg_ROC_train_whole$sensitivities
)

df_roc_test_whole <- data.frame(
  thresholds = avg_ROC_test_whole$thresholds,
  Specificity = avg_ROC_test_whole$specificities,
  FPR = 1 - avg_ROC_test_whole$specificities,
  Sensitivity = avg_ROC_test_whole$sensitivities
)

# dataframe for 95% CI

df_CI_HE <- data.frame(
  TrainLowerRange = quantile(AUC_train_SP[,1], probs = c(0.025, 0.975))[1],
  TrainUpperRange = quantile(AUC_train_SP[,1], probs = c(0.025, 0.975))[2],
  TestLowerRange = quantile(AUC_test_SP[,1], probs = c(0.025, 0.975))[1],
  TestUpperRange = quantile(AUC_test_SP[,1], probs = c(0.025, 0.975))[2]
)

df_CI_HL <- data.frame(
  TrainLowerRange = quantile(AUC_train_SP[,2], probs = c(0.025, 0.975))[1],
  TrainUpperRange = quantile(AUC_train_SP[,2], probs = c(0.025, 0.975))[2],
  TestLowerRange = quantile(AUC_test_SP[,2], probs = c(0.025, 0.975))[1],
  TestUpperRange = quantile(AUC_test_SP[,2], probs = c(0.025, 0.975))[2]
)

df_CI_whole <- data.frame(
  TrainLowerRange = quantile(AUC_train_SP[,3], probs = c(0.025, 0.975))[1],
  TrainUpperRange = quantile(AUC_train_SP[,3], probs = c(0.025, 0.975))[2],
  TestLowerRange = quantile(AUC_test_SP[,3], probs = c(0.025, 0.975))[1],
  TestUpperRange = quantile(AUC_test_SP[,3], probs = c(0.025, 0.975))[2]
)

df_CI <- rbind(df_CI_HE, df_CI_HL, df_CI_whole)
rownames(df_CI) <- c("HE", "HL", "whole")

df_FPR <- as.data.frame(FPR_sen)
colnames(df_FPR) <- c("HE", "HL", "whole")

df_TPR <- as.data.frame(TPR_spe)
colnames(df_TPR) <- c("HE", "HL", "whole")

df_cutoff_fpr <- as.data.frame(cutoff_fpr)
colnames(df_cutoff_fpr) <- c("HE", "HL", "whole")

df_cutoff_tpr <- as.data.frame(cutoff_tpr)
colnames(df_cutoff_tpr) <- c("HE", "HL", "whole")

# parameter dataframe

df_para <- data.frame (
  random_seed = ran_seed,
  cutoff = cutoff,   # 0 : no cutoff, 1 : cutoff at certain concentration of CA19-9
  kind_of_system = kind_of_system,   # 0 : ELISA, 1 : LFA
  n_of_boot = n_of_boot,   # number of bootstraps
  train_ratio = train_ratio,    # ratio for train set
  FPR_rate = FPR_rate,
  TPR_rate = TPR_rate
)

# get the current fold name and generate new folder

current_dir <- getwd()
folder_name <- format(Sys.time(), "%y%m%d_%H%M_")
folder_name <- paste0(folder_name, folder_name_addon)

new_dir <- file.path(current_dir, folder_name)
if (!dir.exists(new_dir)) {
  dir.create(new_dir)
}

setwd(new_dir)

write.xlsx(tdata_HE[HE_idx,], "train_for_HE.xlsx", rowNames = FALSE)
write.xlsx(tdata_HL[HL_idx,], "train_for_HL.xlsx", rowNames = FALSE)
write.xlsx(tdata_whole[whole_idx,], "train_for_all.xlsx", rowNames = FALSE)

write.xlsx(tdata_HE[-HE_idx,], "test_for_HE.xlsx", rowNames = FALSE)
write.xlsx(tdata_HL[-HL_idx,], "test_for_HL.xlsx", rowNames = FALSE)
write.xlsx(tdata_whole[-whole_idx,], "test_for_all.xlsx", rowNames = FALSE)

write.xlsx(AUC_train_SP[ , 1, drop=FALSE], "AUC_train_HE.xlsx", rowNames = FALSE)
write.xlsx(AUC_train_SP[ , 2, drop=FALSE], "AUC_train_HL.xlsx", rowNames = FALSE)
write.xlsx(AUC_train_SP[ , 3, drop=FALSE], "AUC_train_all.xlsx", rowNames = FALSE)

write.xlsx(AUC_test_SP[ , 1, drop=FALSE], "AUC_test_HE.xlsx", rowNames = FALSE)
write.xlsx(AUC_test_SP[ , 2, drop=FALSE], "AUC_test_HL.xlsx", rowNames = FALSE)
write.xlsx(AUC_test_SP[ , 3, drop=FALSE], "AUC_test_all.xlsx", rowNames = FALSE)

write.xlsx(df_Pvalues, "P values (Status from P).xlsx", rowNames = FALSE)

write.xlsx(df_BA_HE, "Bland-Altman_HE.xlsx", rowNames = FALSE)
write.xlsx(df_BA_HL, "Bland-Altman_HL.xlsx", rowNames = FALSE)
write.xlsx(df_BA_whole, "Bland-Altman_all.xlsx", rowNames = FALSE)

write.xlsx(df_roc_train_HE, "ROC_train_HE.xlsx", rowNames = FALSE)
write.xlsx(df_roc_train_HL, "ROC_train_HL.xlsx", rowNames = FALSE)
write.xlsx(df_roc_train_whole, "ROC_train_all.xlsx", rowNames = FALSE)

write.xlsx(df_roc_test_HE, "ROC_test_HE.xlsx", rowNames = FALSE)
write.xlsx(df_roc_test_HL, "ROC_test_HL.xlsx", rowNames = FALSE)
write.xlsx(df_roc_test_whole, "ROC_test_all.xlsx", rowNames = FALSE)

write.xlsx(df_Delong, "Delong.xlsx", rowNames = TRUE)
write.xlsx(df_Delong_Pc, "Delong_only for JKim (P to conc).xlsx", rowNames = TRUE)
write.xlsx(df_CI, "CI range.xlsx", rowNames = TRUE)

write.xlsx(df_FPR, "Sensitivities at FPR less than FPR_rate.xlsx", rowNames = TRUE)
write.xlsx(df_TPR, "Specificities at TPR more than TPR_rate.xlsx", rowNames = TRUE)
write.xlsx(df_cutoff_fpr,"Cutoffs (TPR) at FPR less than FPR_rate.xlsx", rowNames = TRUE)
write.xlsx(df_cutoff_tpr,"Cutoffs (FPR) at TPR more than TPR_rate.xlsx", rowNames = TRUE)

write.xlsx(df_para, "Parameters.xlsx", rowNames = TRUE)

# write the confusion matrix as xlsx file

cm_table_fpr_HE <- as.data.frame(conf_matrix_fpr_HE$table)
cm_table_fpr_HL <- as.data.frame(conf_matrix_fpr_HL$table)
cm_table_fpr_whole <- as.data.frame(conf_matrix_fpr_whole$table)
cm_table_tpr_HE <- as.data.frame(conf_matrix_tpr_HE$table)
cm_table_tpr_HL <- as.data.frame(conf_matrix_tpr_HL$table)
cm_table_tpr_whole <- as.data.frame(conf_matrix_tpr_whole$table)

cm_stats_fpr_HE <- as.data.frame(as.list(conf_matrix_fpr_HE$overall))
cm_stats_fpr_HL <- as.data.frame(as.list(conf_matrix_fpr_HL$overall))
cm_stats_fpr_whole <- as.data.frame(as.list(conf_matrix_fpr_whole$overall))
cm_stats_tpr_HE <- as.data.frame(as.list(conf_matrix_tpr_HE$overall))
cm_stats_tpr_HL <- as.data.frame(as.list(conf_matrix_tpr_HL$overall))
cm_stats_tpr_whole <- as.data.frame(as.list(conf_matrix_tpr_whole$overall))

cm_class_fpr_HE <- as.data.frame(as.list(conf_matrix_fpr_HE$byClass))
cm_class_fpr_HL <- as.data.frame(as.list(conf_matrix_fpr_HL$byClass))
cm_class_fpr_whole <- as.data.frame(as.list(conf_matrix_fpr_whole$byClass))
cm_class_tpr_HE <- as.data.frame(as.list(conf_matrix_tpr_HE$byClass))
cm_class_tpr_HL <- as.data.frame(as.list(conf_matrix_tpr_HL$byClass))
cm_class_tpr_whole <- as.data.frame(as.list(conf_matrix_tpr_whole$byClass))

wb <- createWorkbook()
addWorksheet(wb, "Confusion Matrix")
addWorksheet(wb, "Overall Stats")
addWorksheet(wb, "Class Stats")

writeData(wb, "Confusion Matrix", cm_table_fpr_HE)
writeData(wb, "Overall Stats", cm_stats_fpr_HE)
writeData(wb, "Class Stats", cm_class_fpr_HE)
saveWorkbook(wb, "confusion_matrix_fpr_HE.xlsx", overwrite = TRUE)

writeData(wb, "Confusion Matrix", cm_table_fpr_HL)
writeData(wb, "Overall Stats", cm_stats_fpr_HL)
writeData(wb, "Class Stats", cm_class_fpr_HL)
saveWorkbook(wb, "confusion_matrix_fpr_HL.xlsx", overwrite = TRUE)

writeData(wb, "Confusion Matrix", cm_table_fpr_whole)
writeData(wb, "Overall Stats", cm_stats_fpr_whole)
writeData(wb, "Class Stats", cm_class_fpr_whole)
saveWorkbook(wb, "confusion_matrix_fpr_all.xlsx", overwrite = TRUE)

writeData(wb, "Confusion Matrix", cm_table_tpr_HE)
writeData(wb, "Overall Stats", cm_stats_tpr_HE)
writeData(wb, "Class Stats", cm_class_tpr_HE)
saveWorkbook(wb, "confusion_matrix_tpr_HE.xlsx", overwrite = TRUE)

writeData(wb, "Confusion Matrix", cm_table_tpr_HL)
writeData(wb, "Overall Stats", cm_stats_tpr_HL)
writeData(wb, "Class Stats", cm_class_tpr_HL)
saveWorkbook(wb, "confusion_matrix_tpr_HL.xlsx", overwrite = TRUE)

writeData(wb, "Confusion Matrix", cm_table_tpr_whole)
writeData(wb, "Overall Stats", cm_stats_tpr_whole)
writeData(wb, "Class Stats", cm_class_tpr_whole)
saveWorkbook(wb, "confusion_matrix_tpr_all.xlsx", overwrite = TRUE)

# save the model for P ~ conc as RData file

save(model_P_conc_HE, file = "model for conc to P_HE.RData")
save(model_P_conc_HL, file = "model for conc to P_HL.RData")
save(model_P_conc_whole, file = "model for conc to P_all.RData")




setwd(current_dir)


