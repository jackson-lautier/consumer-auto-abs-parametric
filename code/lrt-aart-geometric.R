#2017 - 50 month loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-38mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-38mo-trapezoid-dim.csv')

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

#find parameters
#init = c(0.98912, 1.219) #from visual analysis; need increasing hazard (50mo-2)
# init = c(0.9874, 1.168) #from visual analysis; need increasing hazard (50mo-2)
# 
# p_hat =
#   optim(init,P_constraint,method="L-BFGS-B",
#       lower=c(0.987,1.15),
#       upper=c(0.990,1.225))$par
# 
# G_hat =
#   sapply(c((Delta + 1):(M+Delta)),
#          g_tau_hat,
#          p_input = p_hat)


#likelihood using geometric
theta_hat = thm_formulas(obs_data)
l.0 = -log_like_fn(theta_hat)

#unrestricted likelihood
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)

#f + #g – 2 – (#g – 1 + #p).
deg_free =
  length(U_est) + length(G_est) - 2 - (length(theta_hat) - 1)

#simplified version ---> check formula
#deg_free = xi - (Delta + 1) - length(p_hat)

#critical value
qchisq(0.95, deg_free)

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)













