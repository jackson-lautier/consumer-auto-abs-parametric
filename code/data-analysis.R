#packages
require('ggplot2')
require('extrafont') #may need to load fonts
require('cowplot')
require('latex2exp')
require('microbenchmark') #used for timing purposes


#dependent files

# supporting files:
# ".\raw-data\aart173_compiledr.csv"
# ".\raw-data\aart193_compiledr.csv"
#
# "./code/default_time.R"
# "./code/rt_geom_sim_studies_LT_formulas.R"
# "./code/rt_geom_sim_studies_RC_formulas.R"
# "./code/hsim-formulas.R"
# "./code/ime-formulas.R"
# "./code/rt_discrete_weibull_sim_studies_RC_formulas.R"

#supplemental files available as reference

# "discrete-weibull-optim-aart-2017-25mo.R"
# "discrete-weibull-optim-aart-2017-50mo.R"

#where processed data will be stored
dir.create('./results/')



################################################################################
################################################################################
# SECTION 2, TABLE 2
################################################################################
################################################################################
rm(list=ls())

source("./code/default_time.R")

#2017

path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

nrow(aart)

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)

mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

summary(aart$reportingPeriodBeginningLoanBalanceAmount)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

#2019

path = "./raw-data/"
aart <- read.csv(paste(path,'aart193_compiledr.csv',sep=""))

nrow(aart)

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)

mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

summary(aart$reportingPeriodBeginningLoanBalanceAmount)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

################################################################################
################################################################################
# SECTION 4, NUMERIC OPTIMIZAION FAIL
#
# NOTE: RUN-TIME APPROX 25 MINS FOR NUMERIC OPTIMIZATION PROCEDURE
################################################################################
################################################################################

rm(list=ls())

#problem set-up
m = 20
Delta = 0
omega = 24

p = 0.05

G = c(sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.4,
      sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_LT_formulas.R')

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_3.1", "thm_3.4")
colnames(disp_mat) = c("p", paste("g", c((Delta + 1):(Delta + m)), sep=""))

#simulate data from h_star
n = 1000
set.seed(1)
sample = sapply(runif(n), h_sim)

#observed data
Yi = sample[2,]
Xi = sample[1,]

obs_data = data.frame(
  "Yi" = Yi,
  "Xi" = Xi
)

#direct numeric optimization

#set up constraints
num_v = (m + Delta) - (Delta + 1) + 1

ui = matrix(NA, nrow = (2 + 2 * num_v + 2), ncol = num_v + 1)

ui[1,] = c(1, rep(0,num_v))
ui[2,] = c(-1, rep(0,num_v))

for(k in c(1:num_v)){
  row = rep(0,(num_v + 1))
  row[k + 1] = 1
  ui[3 + 2 * (k - 1), ] = row
  ui[3 + 2 * (k - 1) + 1, ] = -row
}

ui[2 + 2 * num_v + 1, ] = c(0, rep(1,num_v))
ui[2 + 2 * num_v + 2, ] = c(0, rep(-1,num_v))

ci = c(rep(c(0,-1), (num_v + 1)), c(0.999,-1.001))

init = c(0.5, rep(1/(num_v), num_v))
all(ui %*% init - ci > 0)

start.time <- Sys.time()
cpu_optim = constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#Time difference of 25.0694 mins

disp_mat["cpu_optim",] = cpu_optim$par

#stationary point theorem
start.time <- Sys.time()
p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#Time difference of 0.4090872 secs

disp_mat["thm_3.1", ] = c(p_hat, G_hat)

#closed form MLE solutions

#for calculating G-MLE
start.time <- Sys.time()
theta_est = thm_formulas(obs_data)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#Time difference of 0.1050742 secs

disp_mat["thm_3.4",] = theta_est

#summary
disp_mat

################################################################################
################################################################################
# SECTION 4, TABLE 3
################################################################################
################################################################################

################################################################################
#NO CENSORING
################################################################################

rm(list=ls())

#problem set-up
m = 20; M = 20
Delta = 0
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
#epsilon = 38
#tau = epsilon - (m + Delta + 1)

p = 0.05

G = c(sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.4,
      sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_LT_formulas.R')
source('./code/hsim-formulas.R')

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur = append(cur,
                 (psi(u,v,THETA)^2) * h_star(u,v,THETA))
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur2 = append(cur2,
                  dpsi_dp(u,v,THETA) * h_star(u,v,THETA))
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)


for(n in samp_sizes){
  
  e = "none"
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Xi = sample[1,]
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Xi" = Xi
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Xi = obs_data$Xi[i]
      Yi = obs_data$Yi[i]
      
      Vi = append(Vi, (psi(Xi, Yi, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(e, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tab3-cens-none.csv")
  
}


################################################################################
# e = 26
################################################################################

rm(list=ls())
#problem set-up

m = 20; M = 20
Delta = 0
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 26
tau = epsilon - (m + Delta + 1)

p = 0.05

G = c(sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.4,
      sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tab3-cens-e26.csv")
  
}

################################################################################
# e = 32
################################################################################

rm(list=ls())
#problem set-up

m = 20; M = 20
Delta = 0
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 32
tau = epsilon - (m + Delta + 1)

p = 0.05

G = c(sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.4,
      sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tab3-cens-e32.csv")
  
}

################################################################################
# e = 32
################################################################################

rm(list=ls())
#problem set-up

m = 20; M = 20
Delta = 0
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 38
tau = epsilon - (m + Delta + 1)

p = 0.05

G = c(sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.4,
      sapply(c(0:9), dbinom, prob = 0.35, size = 9) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tab3-cens-e38.csv")
  
}

################################################################################
################################################################################
# SECTION 5, FIGURE 2
################################################################################
################################################################################

################################################################################
# AART-2017 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-25mo-trapezoid-dim.csv')

Delta = trap_param$delta
M = max(obs_data$Y) - Delta
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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

Delta = trap_param$delta 
#M = trap_param$m
M = max(obs_data$Y) - Delta 
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 
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
M = max(obs_data$Y) - Delta 
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 
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

################################################################################
################################################################################
# SECTION 5, TABLE 4
################################################################################
################################################################################

################################################################################
# GEOMETRIC
################################################################################

################################################################################
# AART-2017 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-25mo-trapezoid-dim.csv')

Delta = trap_param$delta
M = max(obs_data$Y) - Delta 
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1 
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

#likelihood using geometric
theta_hat = thm_formulas(obs_data)
l.0 = -log_like_fn(theta_hat)

#unrestricted likelihood
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)

deg_free =
  length(U_est) + length(G_est) - 2 - (length(theta_hat) - 1)

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)

################################################################################
# AART-2019 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2019-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2019-25mo-trapezoid-dim.csv')

Delta = trap_param$delta
M = max(obs_data$Y) - Delta 
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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

#likelihood using geometric
theta_hat = thm_formulas(obs_data)
l.0 = -log_like_fn(theta_hat)

#unrestricted likelihood
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)

deg_free =
  length(U_est) + length(G_est) - 2 - (length(theta_hat) - 1)

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)


################################################################################
# AART-2019 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-50mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-50mo-trapezoid-dim.csv')

Delta = trap_param$delta + 1
M = max(obs_data$Y) - Delta 
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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

#likelihood using geometric
theta_hat = thm_formulas(obs_data)
l.0 = -log_like_fn(theta_hat)

#unrestricted likelihood
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)

deg_free =
  length(U_est) + length(G_est) - 2 - (length(theta_hat) - 1)

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)

################################################################################
# DISCRETE WEIBULL
################################################################################
# AART-2017 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-25mo-trapezoid-dim.csv')

Delta = trap_param$delta
M = max(obs_data$Y) - Delta
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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
#initial values from visual analysis
init = c(0.9899, 1.34) #see discrete-weibull-optim-aart-2017-25mo.R (2017-25mo)

p_hat =
  optim(init,P_constraint,method="L-BFGS-B",
        lower=c(0.985,1.32),
        upper=c(0.995,1.36))$par

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

deg_free =
  length(U_est) + length(G_est) - 2 - (length(G_hat) - 1 + length(p_hat))

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)

################################################################################
# AART-2019 ~ 25MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2019-25mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2019-25mo-trapezoid-dim.csv')

Delta = trap_param$delta
M = max(obs_data$Y) - Delta
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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
#initial values from visual analysis
init = c(0.98774, 1.3888) #visual analysis; details below

p_hat =
  optim(init,P_constraint,method="L-BFGS-B",
        lower=c(0.985,1.37),
        upper=c(0.990,1.40))$par

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

deg_free =
  length(U_est) + length(G_est) - 2 - (length(G_hat) - 1 + length(p_hat))

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)

################################################################################
# AART-2017 ~ 50MO Loans
rm(list=ls())

obs_data = read.csv('./data-clean/aart-2017-50mo.csv')
obs_data = obs_data[,-1]

trap_param = read.csv('./data-clean/aart-2017-50mo-trapezoid-dim.csv')

Delta = trap_param$delta + 1 #for 50mo data
M = max(obs_data$Y) - Delta
epsilon = trap_param$e
tau = epsilon - (M + Delta + 1)
omega = trap_param$omega
xi = min(omega, epsilon - 1)

minU = Delta + 1
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
init = c(0.98964, 1.2316) #discrete-weibull-optim-aart-2017-50mo.R

p_hat =
  optim(init,P_constraint,method="L-BFGS-B",
        lower=c(0.985,1.21),
        upper=c(0.995,1.25))$par

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

deg_free =
  length(U_est) + length(G_est) - 2 - (length(G_hat) - 1 + length(p_hat))

#p.value
pchisq(Q, deg_free, lower.tail = FALSE)


################################################################################
################################################################################
# SECTION F, TABLE F1
################################################################################
################################################################################

#see 'disp_mat' in SECTION 4, NUMERIC OPTIMIZAION FAIL

################################################################################
################################################################################
# SECTION F, TABLE F2
################################################################################
################################################################################

################################################################################
#LEFT-TRUNCATION (S3) & PL Geometric
################################################################################

rm(list=ls())

m = 3
Delta = 0
omega = 4

p = 0.3
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_LT_formulas.R')

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_3.1", "thm_3.4")
colnames(disp_mat) = c("p", "g1", "g2", "g3")

#simulate data from h_star
n = 1000
set.seed(1)
sample = sapply(runif(n), h_sim)

#observed data
Yi = sample[2,]
Xi = sample[1,]

obs_data = data.frame(
  "Yi" = Yi,
  "Xi" = Xi
)

#direct numeric optimization

#set up constraints
num_v = (m + Delta) - (Delta + 1) + 1

ui = matrix(NA, nrow = (2 + 2 * num_v + 2), ncol = num_v + 1)

ui[1,] = c(1, rep(0,num_v))
ui[2,] = c(-1, rep(0,num_v))

for(k in c(1:num_v)){
  row = rep(0,(num_v + 1))
  row[k + 1] = 1
  ui[3 + 2 * (k - 1), ] = row
  ui[3 + 2 * (k - 1) + 1, ] = -row
}

ui[2 + 2 * num_v + 1, ] = c(0, rep(1,num_v))
ui[2 + 2 * num_v + 2, ] = c(0, rep(-1,num_v))

ci = c(rep(c(0,-1), (num_v + 1)), c(0.999,-1.001))

init = c(0.5, rep(1/(num_v), num_v))
all(ui %*% init - ci > 0)

cpu_optim = constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci)
disp_mat["cpu_optim",] = cpu_optim$par

#stationary point theorem
p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
disp_mat["thm_3.1", ] = c(p_hat, G_hat)

#closed form MLE solutions

#for calculating G-MLE
disp_mat["thm_3.4",] = thm_formulas(obs_data)

#summary
disp_mat

#timing
microbenchmark(constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci),
               c(optimize(P_constraint, c(0,1), tol = 1e-10)$minimum, 
                 mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)),
               thm_formulas(obs_data),
               times = 100)

################################################################################
#LEFT-TRUNCATION & RIGHT-CENSORING (S4) & PL Geometric
################################################################################

rm(list=ls())

m = 3; M = 3
Delta = 0
omega = 4

epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.3
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_4.1", "cor_4.2.2")
colnames(disp_mat) = c("p", "g1", "g2", "g3")

#simulate data from h_star
n = 1000
set.seed(9999)
sample = sapply(runif(n), h_sim)

#observed data
Yi = sample[2,]
Zi = pmin(sample[1,], sample[2,] + tau)
Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)

obs_data = data.frame(
  "Yi" = Yi,
  "Zi" = Zi,
  "Di" = Di
)

#direct numeric optimization

#set up constraints
num_v = (m + Delta) - (Delta + 1) + 1

ui = matrix(NA, nrow = (2 + 2 * num_v + 2), ncol = num_v + 1)

ui[1,] = c(1, rep(0,num_v))
ui[2,] = c(-1, rep(0,num_v))

for(k in c(1:num_v)){
  row = rep(0,(num_v + 1))
  row[k + 1] = 1
  ui[3 + 2 * (k - 1), ] = row
  ui[3 + 2 * (k - 1) + 1, ] = -row
}

ui[2 + 2 * num_v + 1, ] = c(0, rep(1,num_v))
ui[2 + 2 * num_v + 2, ] = c(0, rep(-1,num_v))

ci = c(rep(c(0,-1), (num_v + 1)), c(0.999,-1.001))

init = c(0.5, rep(1/(num_v), num_v))
all(ui %*% init - ci > 0)

cpu_optim = constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci)
disp_mat["cpu_optim",] = cpu_optim$par

#stationary point theorem
p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
disp_mat["thm_4.1", ] = c(p_hat, G_hat)

#closed form MLE solutions

#for calculating G-MLE
disp_mat["cor_4.2.2",] = thm_formulas(obs_data)

#summary
disp_mat

#timing
microbenchmark(constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci),
               c(optimize(P_constraint, c(0,1), tol = 1e-10)$minimum, 
                 mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)),
               thm_formulas(obs_data),
               times = 100,
               unit = "microseconds")

################################################################################
#LEFT-TRUNCATION (S3) & Shifted-Binomial
################################################################################
rm(list=ls())

m = 3
Delta = 0
omega = 4

p = 0.75
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/binomial_sim_studies_LT_formulas.R')

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_3.1", "")
colnames(disp_mat) = c("theta", "g1", "g2", "g3")

#simulate data from h_star
n = 1000
set.seed(1)
sample = sapply(runif(n), h_sim)

#observed data
Yi = sample[2,]
Xi = sample[1,]

obs_data = data.frame(
  "Yi" = Yi,
  "Xi" = Xi
)

#direct numeric optimization

#set up constraints
num_v = (m + Delta) - (Delta + 1) + 1

ui = matrix(NA, nrow = (2 + 2 * num_v + 2), ncol = num_v + 1)

ui[1,] = c(1, rep(0,num_v))
ui[2,] = c(-1, rep(0,num_v))

for(k in c(1:num_v)){
  row = rep(0,(num_v + 1))
  row[k + 1] = 1
  ui[3 + 2 * (k - 1), ] = row
  ui[3 + 2 * (k - 1) + 1, ] = -row
}

ui[2 + 2 * num_v + 1, ] = c(0, rep(1,num_v))
ui[2 + 2 * num_v + 2, ] = c(0, rep(-1,num_v))

ci = c(rep(c(0,-1), (num_v + 1)), c(0.999,-1.001))

init = c(0.5, rep(1/(num_v), num_v))
all(ui %*% init - ci > 0)

cpu_optim = constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci)
disp_mat["cpu_optim",] = cpu_optim$par

#stationary point theorem
p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
disp_mat["thm_3.1", ] = c(p_hat, G_hat)

#closed form MLE solutions

#for calculating G-MLE
#disp_mat["mle",] = thm_formulas(obs_data)

#summary
disp_mat

#timing
microbenchmark(constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci),
               c(optimize(P_constraint, c(0,1), tol = 1e-10)$minimum, 
                 mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)),
               #thm_formulas(obs_data),
               times = 100,
               unit = "microseconds")

################################################################################
#LEFT-TRUNCATION & RIGHT-CENSORING (S4) & Shifted-Binomial
################################################################################
rm(list=ls())

m = 3
Delta = 0
omega = 4

epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.75
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/binomial_sim_studies_RC_formulas.R')

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_4.1", "")
colnames(disp_mat) = c("theta", "g1", "g2", "g3")

#simulate data from h_star
n = 1000
set.seed(9999)
sample = sapply(runif(n), h_sim)

#observed data
Yi = sample[2,]
Zi = pmin(sample[1,], sample[2,] + tau)
Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)

obs_data = data.frame(
  "Yi" = Yi,
  "Zi" = Zi,
  "Di" = Di
)

#direct numeric optimization

#set up constraints
num_v = (m + Delta) - (Delta + 1) + 1

ui = matrix(NA, nrow = (2 + 2 * num_v + 2), ncol = num_v + 1)

ui[1,] = c(1, rep(0,num_v))
ui[2,] = c(-1, rep(0,num_v))

for(k in c(1:num_v)){
  row = rep(0,(num_v + 1))
  row[k + 1] = 1
  ui[3 + 2 * (k - 1), ] = row
  ui[3 + 2 * (k - 1) + 1, ] = -row
}

ui[2 + 2 * num_v + 1, ] = c(0, rep(1,num_v))
ui[2 + 2 * num_v + 2, ] = c(0, rep(-1,num_v))

ci = c(rep(c(0,-1), (num_v + 1)), c(0.999,-1.001))

init = c(0.5, rep(1/(num_v), num_v))
all(ui %*% init - ci > 0)

cpu_optim = constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci)
disp_mat["cpu_optim",] = cpu_optim$par

#stationary point theorem
p_hat = optimize(P_constraint, c(0,1), tol = 1e-10)$minimum
G_hat = mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)
disp_mat["thm_4.1", ] = c(p_hat, G_hat)

#closed form MLE solutions

#for calculating G-MLE
#disp_mat["mle",] = thm_formulas(obs_data)

#summary
disp_mat

#timing
microbenchmark(constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci),
               c(optimize(P_constraint, c(0,1), tol = 1e-10)$minimum, 
                 mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)),
               #thm_formulas(obs_data),
               times = 100,
               unit = "microseconds")

################################################################################
################################################################################
# SECTION F, TABLE F3
################################################################################
################################################################################

################################################################################
#NO CENSORING
################################################################################

rm(list=ls())
#problem set-up

m = 15; M = 15
Delta = 5
omega = 24

p = 0.45

minU = Delta + 1 
maxU = omega
minV = Delta + 1
maxV = Delta + M

range.1 = c(1:(floor((maxV - minV + 1)/2))) - 1
range.2 = c(1:(floor((maxV - minV + 1)/2) + (maxV - minV + 1) %% 2)) - 1

G = c(sapply(range.1, dbinom, prob = 0.35, size = length(range.1)) * 0.4,
      sapply(range.2, dbinom, prob = 0.35, size = length(range.2)) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_LT_formulas.R')
source('./code/hsim-formulas.R')

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur = append(cur,
                 (psi(u,v,THETA)^2) * h_star(u,v,THETA))
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur2 = append(cur2,
                  dpsi_dp(u,v,THETA) * h_star(u,v,THETA))
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)


for(n in samp_sizes){
  
  e = "none"
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Xi = sample[1,]
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Xi" = Xi
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Xi = obs_data$Xi[i]
      Yi = obs_data$Yi[i]
      
      Vi = append(Vi, (psi(Xi, Yi, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(e, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tabf3-cens-none.csv")
  
}


################################################################################
# e = 25
################################################################################

rm(list=ls())
#problem set-up

m = 15; M = 15
Delta = 5
omega = 24

epsilon = 25
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.45

minU = Delta + 1 
maxU = omega
minV = Delta + 1
maxV = Delta + M

range.1 = c(1:(floor((maxV - minV + 1)/2))) - 1
range.2 = c(1:(floor((maxV - minV + 1)/2) + (maxV - minV + 1) %% 2)) - 1

G = c(sapply(range.1, dbinom, prob = 0.35, size = length(range.1)) * 0.4,
      sapply(range.2, dbinom, prob = 0.35, size = length(range.2)) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tabf3-cens-e25.csv")
  
}

################################################################################
# e = 30
################################################################################

rm(list=ls())
#problem set-up

m = 15; M = 15
Delta = 5
omega = 24

epsilon = 30
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.45

minU = Delta + 1 
maxU = omega
minV = Delta + 1
maxV = Delta + M

range.1 = c(1:(floor((maxV - minV + 1)/2))) - 1
range.2 = c(1:(floor((maxV - minV + 1)/2) + (maxV - minV + 1) %% 2)) - 1

G = c(sapply(range.1, dbinom, prob = 0.35, size = length(range.1)) * 0.4,
      sapply(range.2, dbinom, prob = 0.35, size = length(range.2)) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tabf3-cens-e30.csv")
  
}

################################################################################
# e = 35
################################################################################

rm(list=ls())
#problem set-up

m = 15; M = 15
Delta = 5
omega = 24

epsilon = 35
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.45

minU = Delta + 1 
maxU = omega
minV = Delta + 1
maxV = Delta + M

range.1 = c(1:(floor((maxV - minV + 1)/2))) - 1
range.2 = c(1:(floor((maxV - minV + 1)/2) + (maxV - minV + 1) %% 2)) - 1

G = c(sapply(range.1, dbinom, prob = 0.35, size = length(range.1)) * 0.4,
      sapply(range.2, dbinom, prob = 0.35, size = length(range.2)) * 0.6)

THETA = c(p, G)
source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')


#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

results = matrix(NA, )

samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#####################################
for(n in samp_sizes){
  
  print(n)
  
  #run thru 1,000 replicates
  cov_ind = c()
  p_est_vec = c()
  for(r in c(1:replicates)){
    
    #generate random sample
    #k = (which(samp_sizes == n) - 1) * replicates
    #set.seed(k + r)
    sample = sapply(runif(n), h_sim)
    
    #observed data
    Yi = sample[2,]
    Zi = pmin(sample[1,], sample[2,] + tau)
    Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
    
    obs_data = data.frame(
      "Yi" = Yi,
      "Zi" = Zi,
      "Di" = Di
    )
    
    #estimate parameters
    THETA_est = thm_formulas(obs_data)
    p_est_vec = append(p_est_vec, THETA_est[1])
    
    #Vn estimate
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
      
    }
    
    Vn = (1/n) * sum(Vi)
    sigma = sqrt(1 / Vn)
    
    #CI estimates
    CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
    CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind = append(cov_ind, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  mat_row =  c(epsilon, n, p,
               mean(p_est_vec),
               sd(p_est_vec),
               sqrt(true_var / n),
               sum(cov_ind) / replicates)
  
  results = rbind(results, mat_row)  
  write.csv(results, "./results/tabf3-cens-e35.csv")
  
}

################################################################################
################################################################################
# SECTION F, FIGURE F1
################################################################################
################################################################################

################################################################################
#LEFT-TRUNCATION (S3) & PL Geometric
################################################################################

rm(list=ls())

m = 3
Delta = 0
omega = 4

p = 0.3
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_LT_formulas.R')

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur = append(cur,
                 (psi(u,v,THETA)^2) * h_star(u,v,THETA))
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur2 = append(cur2,
                  dpsi_dp(u,v,THETA) * h_star(u,v,THETA))
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#verify with simulation
replicates = 1000
rep_results = c()
n = 1000

for(r in c(1:replicates)){
  
  set.seed(r)
  sample = sapply(runif(n), h_sim)
  
  #observed data
  Yi = sample[2,]
  Xi = sample[1,]
  
  obs_data = data.frame(
    "Yi" = Yi,
    "Xi" = Xi
  )
  
  rep_results = append(rep_results,
                       thm_formulas(obs_data)[1])
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

mean( sqrt(n) * (rep_results - THETA[1]) )
var( sqrt(n) * (rep_results - THETA[1]) ); true_var
sd( sqrt(n) * (rep_results - THETA[1]) ); sqrt(true_var)

#plot results
df = data.frame("sim_result" = sqrt(n) * (rep_results - THETA[1]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var))) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(expression( sqrt( italic(n) ) * ( hat( italic(p) )[italic(n)] - italic(p)[0] ) ) ) +
  ylab("Density Height") +
  #guides(linetype=guide_legend("")) +
  #scale_linetype(labels = c("gamma-kernel", "nonparametric loess")) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "bottom")

x = seq(from = min(df$sim_result), to = max(df$sim_result),
        by =
          (max(df$sim_result) - min(df$sim_result))/ (length(df$sim_result) - 1))
y = dnorm(x, mean = 0, sd = sqrt(true_var))

facet_data = data.frame("sim_result" = df$sim_result,
                        "true_density" = y,
                        "x_value" = x,
                        "scenario" = "truncated-geometric",
                        "setting" = "left-truncation")

write.csv(facet_data, "./results/fd_lt_rt_geom.csv")


################################################################################
#LEFT-TRUNCATION & RIGHT-CENSORING (S4) & PL Geometric
################################################################################
rm(list=ls())

m = 3; M = 3
Delta = 0
omega = 4

epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.3
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#verify with simulation
replicates = 1000
rep_results = c()
n = 1000

for(r in c(1:replicates)){
  
  set.seed(n + r)
  sample = sapply(runif(n), h_sim)
  
  #observed data
  Yi = sample[2,]
  Zi = pmin(sample[1,], sample[2,] + tau)
  Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
  
  obs_data = data.frame(
    "Yi" = Yi,
    "Zi" = Zi,
    "Di" = Di
  )
  
  rep_results = append(rep_results,
                       thm_formulas(obs_data)[1])
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

mean( sqrt(n) * (rep_results - THETA[1]) )
var( sqrt(n) * (rep_results - THETA[1]) ); true_var
sd( sqrt(n) * (rep_results - THETA[1]) ); sqrt(true_var)

#plot results
df = data.frame("sim_result" = sqrt(n) * (rep_results - THETA[1]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var))) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(expression( sqrt( italic(n) ) * ( hat( italic(p) )[italic(n)] - italic(p)[0] ) ) ) +
  ylab("Density Height") +
  #guides(linetype=guide_legend("")) +
  #scale_linetype(labels = c("gamma-kernel", "nonparametric loess")) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "bottom")

#ggsave("LT_geom_hist.pdf",height=4,width=6,device = cairo_pdf)

x = seq(from = min(df$sim_result), to = max(df$sim_result),
        by =
          (max(df$sim_result) - min(df$sim_result))/ (length(df$sim_result) - 1))
y = dnorm(x, mean = 0, sd = sqrt(true_var))

facet_data = data.frame("sim_result" = df$sim_result,
                        "true_density" = y,
                        "x_value" = x,
                        "scenario" = "truncated-geometric",
                        "setting" = "right-censoring")

write.csv(facet_data, "./results/fd_rc_rt_geom.csv")


################################################################################
#LEFT-TRUNCATION (S3) & Shifted-Binomial
################################################################################
rm(list=ls())

m = 3
Delta = 0
omega = 4

p = 0.75
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/binomial_sim_studies_LT_formulas.R')

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur = append(cur,
                 (psi(u,v,THETA)^2) * h_star(u,v,THETA))
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    cur2 = append(cur2,
                  dpsi_dp(u,v,THETA) * h_star(u,v,THETA))
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#verify with simulation
replicates = 1000
rep_results = c()
n = 1000

for(r in c(1:replicates)){
  
  set.seed(r)
  sample = sapply(runif(n), h_sim)
  
  #observed data
  Yi = sample[2,]
  Xi = sample[1,]
  
  obs_data = data.frame(
    "Yi" = Yi,
    "Xi" = Xi
  )
  
  rep_results = append(rep_results,
                       optimize(P_constraint, c(0,1), tol = 1e-10)$minimum)
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

mean( sqrt(n) * (rep_results - THETA[1]) )
var( sqrt(n) * (rep_results - THETA[1]) ); true_var
sd( sqrt(n) * (rep_results - THETA[1]) ); sqrt(true_var)

#plot results
df = data.frame("sim_result" = sqrt(n) * (rep_results - THETA[1]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var))) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(expression( sqrt( italic(n) ) * ( hat( italic(p) )[italic(n)] - italic(p)[0] ) ) ) +
  ylab("Density Height") +
  #guides(linetype=guide_legend("")) +
  #scale_linetype(labels = c("gamma-kernel", "nonparametric loess")) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "bottom")

#ggsave("LT_geom_hist.pdf",height=4,width=6,device = cairo_pdf)

x = seq(from = min(df$sim_result), to = max(df$sim_result),
        by =
          (max(df$sim_result) - min(df$sim_result))/ (length(df$sim_result) - 1))
y = dnorm(x, mean = 0, sd = sqrt(true_var))

facet_data = data.frame("sim_result" = df$sim_result,
                        "true_density" = y,
                        "x_value" = x,
                        "scenario" = "shifted-binomial",
                        "setting" = "left-truncation")

write.csv(facet_data, "./results/fd_lt_binom.csv")

################################################################################
#LEFT-TRUNCATION & RIGHT-CENSORING (S4) & Shifted-Binomial
################################################################################
rm(list=ls())

m = 3
Delta = 0
omega = 4

epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.75
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/binomial_sim_studies_RC_formulas.R')

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur = append(cur, 0)
      }
      if(prob > 0){
        cur = append(cur, (psi(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

sum(cur)

cur2 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur2 = append(cur2, 0)
      }
      if(prob > 0){
        cur2 = append(cur2, dpsi_dp(v, u, d, THETA) * prob)
      }
    }
  }
}

sum(cur2)

true_var = sum(cur) / (sum(cur2)^2)

#verify with simulation
replicates = 1000
rep_results = c()
n = 1000

for(r in c(1:replicates)){
  
  set.seed(n + r)
  sample = sapply(runif(n), h_sim)
  
  #observed data
  Yi = sample[2,]
  Zi = pmin(sample[1,], sample[2,] + tau)
  Di = ifelse(sample[1,] <= sample[2,] + tau, 1, 0)
  
  obs_data = data.frame(
    "Yi" = Yi,
    "Zi" = Zi,
    "Di" = Di
  )
  
  rep_results = append(rep_results,
                       optimize(P_constraint, c(0,1), tol = 1e-10)$minimum)
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

mean( sqrt(n) * (rep_results - THETA[1]) )
var( sqrt(n) * (rep_results - THETA[1]) ); true_var
sd( sqrt(n) * (rep_results - THETA[1]) ); sqrt(true_var)

#plot results
df = data.frame("sim_result" = sqrt(n) * (rep_results - THETA[1]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var))) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(expression( sqrt( italic(n) ) * ( hat( italic(p) )[italic(n)] - italic(p)[0] ) ) ) +
  ylab("Density Height") +
  #guides(linetype=guide_legend("")) +
  #scale_linetype(labels = c("gamma-kernel", "nonparametric loess")) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "bottom")

#ggsave("LT_geom_hist.pdf",height=4,width=6,device = cairo_pdf)

x = seq(from = min(df$sim_result), to = max(df$sim_result),
        by =
          (max(df$sim_result) - min(df$sim_result))/ (length(df$sim_result) - 1))
y = dnorm(x, mean = 0, sd = sqrt(true_var))

facet_data = data.frame("sim_result" = df$sim_result,
                        "true_density" = y,
                        "x_value" = x,
                        "scenario" = "shifted-binomial",
                        "setting" = "right-censoring")

write.csv(facet_data, "./results/fd_rc_binom.csv")

################################################################################
#Final figure plot for the manuscript
################################################################################

rm(list=ls())

df1 = read.csv('./results/fd_lt_binom.csv')
df1 = df1[,-1]
df2 = read.csv('./results/fd_lt_rt_geom.csv')
df2 = df2[,-1]
df3 = read.csv('./results/fd_rc_binom.csv')
df3 = df3[,-1]
df4 = read.csv('./results/fd_rc_rt_geom.csv')
df4 = df4[,-1]

plot_df = rbind(df1, df2, df3, df4)

ggplot() +
  geom_density(data = plot_df, aes(x = sim_result),
               color = "blue", linetype = "dashed") +
  geom_line(data = plot_df, aes(x = x_value, y = true_density)) +
  facet_grid(setting ~ scenario) +
  xlab(
    expression( 'shifted-binomial:'~
                  sqrt( italic(n) ) * ( hat( theta )[italic(n)] - theta[0] )
                ~ ' & ' ~
                  'PL-geometric:'~
                  sqrt( italic(n) ) * ( hat( italic(p) )[italic(n)] - italic(p)[0] )) )+
  ylab("Density Height") +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        strip.text = element_text(size = 8, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "none")

ggsave("./results/sim_comps.pdf",height=4,width=6,device = cairo_pdf)

################################################################################
################################################################################
# SECTION G, FIGURE G1
################################################################################
################################################################################


################################################################################
################################################################################
# SECTION G, TABLE G2
################################################################################
################################################################################

#see 'theta_hat; p_hat' values in SECTION 5, TABLE 4














