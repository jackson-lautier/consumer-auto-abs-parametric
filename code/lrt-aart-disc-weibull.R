#2017 - 50 month loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-38mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-38mo-trapezoid-dim.csv')

Delta = trap_param$delta #+ 1 #four for 50mo data ---> need to check (and 2017-38mo loans)
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
source("./code/rt_discrete_weibull_sim_studies_RC_formulas.R")

#find parameters
#init = c(0.98388, 1.198) #from visual analysis; need increasing hazard (2017-25mo)
#init = c(0.98912, 1.219) #from visual analysis; need increasing hazard (2017-50mo-2)
#init = c(0.9874, 1.168) #from visual analysis; need increasing hazard (2017-50mo)
#init = c(0.9895, 1.271) #from visual analysis; need increasing hazard (2017-37mo)
#init = c(0.98258,1.2792) #from visual analysis; need increasing hazard (2019-25mo)
#init = c(0.98428, 1.15) #from visual analysis; need increasing hazard (2017-37mo) - super-prime only
init = c(0.992, 1.336) #from visual analysis; need increasing hazard (2017-38mo) - super-prime only

p_hat =
  optim(init,P_constraint,method="L-BFGS-B",
      lower=c(0.990,1.330),
      upper=c(0.995,1.134))$par

G_hat =
  sapply(c((Delta + 1):(M+Delta)),
         g_tau_hat,
         p_input = p_hat)


#likelihood using discrete weibull
theta_hat = c(p_hat, G_hat)
l.0 = -log_like_fn(theta_hat)

#unrestricted likelihood
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)

#f + #g – 2 – (#g – 1 + #p).
deg_free =
  length(U_est) + length(G_est) - 2 - (length(G_hat) - 1 + length(p_hat))

#simplified version ---> check formula
#deg_free = xi - (Delta + 1) - length(p_hat)

#critical value
qchisq(0.95, deg_free)

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)













