library("MASS")


##############################3(a)##############################
# Assign values to parameters
a = 2
b = 0
n = 500
sigma = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = T)
mju = c(0,0,0)
set.seed(1)

# Construct the complete dataset
z = mvrnorm(n = n, mu = mju, Sigma = sigma)
z_1 = z[, 1]; z_2 = z[, 2]; z_3 = z[, 3]
y_1 = 1 + z_1; y_2 = 5 + 2*z_1 + z_2
dataset = data.frame(Y_1 = y_1, Y_2 = y_2)

# Find the missingness
missingness = a*(y_1-1) + b*(y_2-5) + z_3 < 0
y_2_obs = y_2[which(missingness)]

# Construct the observed dataset
miss_value = replace(missingness, missingness == TRUE, NA)
dataset_missing = dataset
dataset_missing[, 2] = y_2*(1 - miss_value)

# Draw the plot for marginal distribution of Y_2 for the complete and observed data
plot(density(dataset$Y_2), lwd = 2, col = "blue", xlab = "Y_2", main = "MAR", ylim = c(0, 0.3))
lines(density(y_2_obs), lwd = 2, col = "red")
legend(7, 0.3, legend = c("Complete data", "Observed data"),col = c("blue", "red"),
       lty = c(1,1), lwd = c(2,2), bty ="n")

# T-test with null hypothesis of equality of mean
mean(dataset$Y_2)
sd(dataset$Y_2)
mean(y_2_obs)
sd(y_2_obs)
t.test(dataset$Y_2, y_2_obs, paired = FALSE, var.equal = FALSE)


##############################3(b)##############################
# Fit a simple linear regression model
fity2 = lm(Y_2 ~ Y_1, data = dataset_missing)
summary(fity2)

# Complete the dataset using SRI
predy2 = predict(fity2, newdata = dataset_missing) + rnorm(n, 0, sigma(fity2))
y_2_completed = ifelse(is.na(dataset_missing$Y_2) == TRUE, predy2, dataset_missing$Y_2)
dataset_completed = dataset_missing
dataset_completed[, 2] = y_2_completed

plot(density(dataset$Y_2), lwd = 2, col = "blue", xlab = "Y_2", main = "MAR", ylim = c(0, 0.3))
lines(density(y_2_obs), lwd = 2, col = "red")
lines(density(dataset_completed$Y_2), lwd = 2, col = "green")
legend(7, 0.3, legend = c("Observed data", "Complete data", "Completed data"),col = c("red", "blue", "green"),
       lty = c(1,1), lwd = c(2,2), bty ="n")

mean(y_2_completed)
sd(y_2_completed)
t.test(dataset$Y_2, y_2_completed, paired = FALSE, var.equal = FALSE)


##############################3(c)##############################
# Change values of parameters
a = 0
b = 2
n = 500
set.seed(1)

z = mvrnorm(n = n, mu = mju, Sigma = sigma)
z_1 = z[, 1]; z_2 = z[, 2]; z_3 = z[, 3]
y_1 = 1 + z_1; y_2 = 5 + 2*z_1 + z_2
dataset = data.frame(Y_1 = y_1, Y_2 = y_2)

missingness = a*(y_1-1) + b*(y_2-5) + z_3 < 0
y_2_obs = y_2[which(missingness)]

miss_value = replace(missingness, missingness == TRUE, NA)
dataset_missing = dataset
dataset_missing[, 2] = y_2*(1 - miss_value)

plot(density(dataset$Y_2), lwd = 2, col = "blue", xlab = "Y_2", main = "MNAR", ylim = c(0, 0.3))
lines(density(y_2_obs), lwd = 2, col = "red")
legend(7, 0.3, legend = c("Complete data", "Observed data"),col = c("blue", "red"),
       lty = c(1,1), lwd = c(2,2), bty ="n")

mean(dataset$Y_2)
sd(dataset$Y_2)
mean(y_2_obs)
sd(y_2_obs)
t.test(dataset$Y_2, y_2_obs, paired = FALSE, var.equal = FALSE)


##############################3(d)##############################
fity2 = lm(Y_2 ~ Y_1, data = dataset_missing)
summary(fity2)

predy2 = predict(fity2, newdata = dataset_missing) + rnorm(n, 0, sigma(fity2))
y_2_completed = ifelse(is.na(dataset_missing$Y_2) == TRUE, predy2, dataset_missing$Y_2)
dataset_completed = dataset_missing
dataset_completed[, 2] = y_2_completed

plot(density(dataset$Y_2), lwd = 2, col = "blue", xlab = "Y_2", main = "MNAR", ylim = c(0, 0.3))
lines(density(y_2_obs), lwd = 2, col = "red")
lines(density(dataset_completed$Y_2), lwd = 2, col = "green")
legend(7, 0.3, legend = c("Observed data", "Complete data", "Completed data"),col = c("red", "blue", "green"),
       lty = c(1,1), lwd = c(2,2), bty ="n")

mean(y_2_completed)
sd(y_2_completed)
t.test(dataset$Y_2, y_2_completed, paired = FALSE, var.equal = FALSE)


##############################4(a)##############################
load("D:/databp.Rdata")
summary(databp)

# Complete case analysis
ind = which(is.na(databp$recovtime) == FALSE)
# Calculate the mean of recovery time under complete case analysis
mean_recov = mean(databp$recovtime, na.rm = TRUE)
# Calculate the standard error of recovery time under complete case analysis
se_recov = sd(databp$recovtime, na.rm = TRUE)/sqrt(length(ind))
# Calculate the correlation between the recovery time and the dose ......
cov_recov_dose = cor(databp$recovtime, databp$logdose, use = "complete")
# Calculate the correlation between the recovery time and the blood pressure ......
cov_recov_blood = cor(databp$recovtime, databp$bloodp, use = "complete")
mean_recov; se_recov; cov_recov_dose; cov_recov_blood


##############################4(b)##############################
# Mean imputation
bp_mi = ifelse(is.na(databp$recovtime) == TRUE, mean_recov, databp$recovtime)
n = length(bp_mi)
mi_mean_recov = mean(bp_mi)
mi_se_recov = sd(bp_mi)/sqrt(n)
mi_cov_recov_dose = cor(bp_mi, databp$logdose)
mi_cov_recov_blood = cor(bp_mi, databp$bloodp)
mi_mean_recov; mi_se_recov; mi_cov_recov_dose; mi_cov_recov_blood


##############################4(c)##############################
# Mean regression imputation
fit_bp = lm(recovtime ~ logdose + bloodp, data = databp)
summary(fit_bp)
pred_ri = predict(fit_bp, newdata = databp)
bp_ri = ifelse(is.na(databp$recovtime) == TRUE, pred_ri, databp$recovtime)
n = length(bp_ri)
ri_mean_recov = mean(bp_ri)
ri_se_recov = sd(bp_ri)/sqrt(n)
ri_cov_recov_dose = cor(bp_ri, databp$logdose)
ri_cov_recov_blood = cor(bp_ri, databp$bloodp)
ri_mean_recov; ri_se_recov; ri_cov_recov_dose; ri_cov_recov_blood


##############################4(d)##############################
# Stochastic regression imputation
set.seed(1)
pred_sri = predict(fit_bp, newdata = databp) + rnorm(n, 0, sigma(fit_bp))
bp_sri = ifelse(is.na(databp$recovtime) == TRUE, pred_sri, databp$recovtime)
n = length(bp_sri)
sri_mean_recov = mean(bp_sri)
sri_se_recov = sd(bp_sri)/sqrt(n)
sri_cov_recov_dose = cor(bp_sri, databp$logdose)
sri_cov_recov_blood = cor(bp_sri, databp$bloodp)
sri_mean_recov; sri_se_recov; sri_cov_recov_dose; sri_cov_recov_blood

# Draw a simple diagnostic plot
plot(databp, lwd = 2)


##############################4(e)##############################
# Index for missing values
mis = which(is.na(databp$recovtime) == TRUE)

# Construct the completed data
databp_pmm = databp
for (i in mis) {
  # Find the index of donor
  donor = as.integer(names(which.min((pred_ri[ind] - pred_ri[i])^2)))
  databp_pmm$recovtime[i] = databp$recovtime[donor]
  # Print out the imputation information
  print(paste0("The ", i,"'th recovery time is imputed by the ", donor, "'th recovery time"))
}

# Print out the completed dataset
databp_pmm

n = length(databp_pmm$recovtime)
pmm_mean_recov = mean(databp_pmm$recovtime)
pmm_se_recov = sd(databp_pmm$recovtime)/sqrt(n)
pmm_cov_recov_dose = cor(databp_pmm$recovtime, databp_pmm$logdose)
pmm_cov_recov_blood = cor(databp_pmm$recovtime, databp_pmm$bloodp)
pmm_mean_recov; pmm_se_recov; pmm_cov_recov_dose; pmm_cov_recov_blood

