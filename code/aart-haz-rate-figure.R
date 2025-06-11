dir.create("./results/") #to store results

require('ggplot2')
require('extrafont') #may need to load fonts
require('cowplot')

################################################################################
# AART-2017 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-25mo-trapezoid-dim.csv')

Delta = trap_param$delta #+ 1 #four for 2017-50mo data ---> need to check (also 2017-38mo)
#M = trap_param$m
M = max(obs_data$Y) - Delta #for 36mo loans bc of max cap
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 #should smallest X in data? ---> need to check
maxU = xi
minV = Delta + 1
maxV = Delta + M

obs_data$D = 1 - obs_data$C
names(obs_data) = c("Zi", "Yi", "Ci", "Di")

#censoring check for omega
obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0),
                     1,
                     obs_data$Di)

source("./code/ime-formulas.R")
source("./code/rt_geom_sim_studies_RC_formulas.R")

z = sort(unique(obs_data$Zi))
n = nrow(obs_data)

lam_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  lam_hat = append(lam_hat,lnx(i))
}

Var_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))

theta_est = thm_formulas(obs_data)
p_hat = theta_est[1]

p + geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                     y = p_hat, yend = p_hat), color = "red",
                 linetype = "dashed")

#Vn estimate
Vi = c()
for(i in c(1:n)){
  
  Zi = obs_data$Zi[i]
  Yi = obs_data$Yi[i]
  Di = obs_data$Di[i]
  
  Vi = append(Vi, (psi(Yi, Zi, Di, theta_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = theta_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = theta_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)
#(0.02262015, 0.03991475)

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")

df$year_lt = "AART-2017-25M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1725.csv")

################################################################################
# AART-2019 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2019-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2019-25mo-trapezoid-dim.csv')

Delta = trap_param$delta #+ 1 #four for 2017-50mo data ---> need to check (also 2017-38mo)
#M = trap_param$m
M = max(obs_data$Y) - Delta #for 36mo loans bc of max cap
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 #should smallest X in data? ---> need to check
maxU = xi
minV = Delta + 1
maxV = Delta + M

obs_data$D = 1 - obs_data$C
names(obs_data) = c("Zi", "Yi", "Ci", "Di")

#censoring check for omega
obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0),
                     1,
                     obs_data$Di)

source("./code/ime-formulas.R")
source("./code/rt_geom_sim_studies_RC_formulas.R")

z = sort(unique(obs_data$Zi))
n = nrow(obs_data)

lam_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  lam_hat = append(lam_hat,lnx(i))
}

Var_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))

theta_est = thm_formulas(obs_data)
p_hat = theta_est[1]

p + geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                     y = p_hat, yend = p_hat), color = "red",
                 linetype = "dashed")

#Vn estimate
Vi = c()
for(i in c(1:n)){
  
  Zi = obs_data$Zi[i]
  Yi = obs_data$Yi[i]
  Di = obs_data$Di[i]
  
  Vi = append(Vi, (psi(Yi, Zi, Di, theta_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = theta_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = theta_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)
#(0.03374927, 0.05242667)

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")

df$year_lt = "AART-2019-25M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1925.csv")

################################################################################
# AART-2017 ~ 50MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-50mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-50mo-trapezoid-dim.csv')

Delta = trap_param$delta + 1 #four for 2017-50mo data bc min Yi = 5
#M = trap_param$m
M = max(obs_data$Y) - Delta #for 36mo loans bc of max cap
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 #should smallest X in data? ---> need to check
maxU = xi
minV = Delta + 1
maxV = Delta + M

obs_data$D = 1 - obs_data$C
names(obs_data) = c("Zi", "Yi", "Ci", "Di")

#censoring check for omega
obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0),
                     1,
                     obs_data$Di)

source("./code/ime-formulas.R")
source("./code/rt_geom_sim_studies_RC_formulas.R")

z = sort(unique(obs_data$Zi))
n = nrow(obs_data)

lam_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  lam_hat = append(lam_hat,lnx(i))
}

Var_hat = vector()
for (i in c((Delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((Delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))

theta_est = thm_formulas(obs_data)
p_hat = theta_est[1]

p + geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                     y = p_hat, yend = p_hat), color = "red",
                 linetype = "dashed")

#Vn estimate
Vi = c()
for(i in c(1:n)){
  
  Zi = obs_data$Zi[i]
  Yi = obs_data$Yi[i]
  Di = obs_data$Di[i]
  
  Vi = append(Vi, (psi(Yi, Zi, Di, theta_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = theta_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = theta_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)
#(0.02274845, 0.02794531)

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")

df$year_lt = "AART-2017-50M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1750.csv")

################################################################################
# construct the figure
################################################################################
rm(list=ls())

df1 = read.csv('./results/df1725.csv')
df1 = df1[,-1]
df2 = read.csv('./results/df1925.csv')
df2 = df2[,-1]
df3 = read.csv('./results/df1750.csv')
df3 = df3[,-1]

plot_df = rbind(df1, df2)
plot_df = plot_df[plot_df$age >= 8,]

p1 <-
  ggplot() +
  geom_line(data=plot_df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  geom_line(data=plot_df, aes(x=age, y=phat), color="red", linetype = "dashed") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=p_low, ymax=p_upp),
              fill="red", alpha=0.35) +
  facet_grid(cols=vars(year_lt)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))# +

plot_df = df3
#plot_df = plot_df[-1,] #remove NA row
#plot_df = plot_df[plot_df$age >= 7,]

p2 <-
  ggplot() +
  geom_line(data=plot_df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  geom_line(data=plot_df, aes(x=age, y=phat), color="red", linetype = "dashed") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=p_low, ymax=p_upp),
              fill="red", alpha=0.35) +
  facet_grid(cols=vars(year_lt)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))# +

set_null_device(cairo_pdf)
plot_grid(p1, p2, nrow=2)
ggsave("./results/aart_comp.pdf",height=4,width=6,device=cairo_pdf)
file.remove('./Rplot001.pdf')










