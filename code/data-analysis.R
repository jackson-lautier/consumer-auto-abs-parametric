#packages



#dependent files

# supporting files:
# ".\raw-data\aart173_compiledr.csv"
# ".\raw-data\aart193_compiledr.csv"
#
# "./code/default_time.R"
# "./code/rt_geom_sim_studies_LT_formulas.R"
# "./code/rt_geom_sim_studies_RC_formulas.R"
# './code/hsim-formulas.R'


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
# SECTION F, TABLE F1
################################################################################
################################################################################

#see 'disp_mat' in SECTION 4, NUMERIC OPTIMIZAION FAIL















