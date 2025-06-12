################################################################################
################################################################################
# Table 1: Numeric Validation and Performance Summary
################################################################################
################################################################################

require('ggplot2')
require('extrafont') #may need to load fonts
require('latex2exp')
require('microbenchmark') #used for timing purposes

#please run sequentially downwards
#includes all numeric, simulations studies from section 5

dir.create("./results/") #to store results


################################################################################
################################################################################
# Robustness Simulation Study (non-zero Delta)
################################################################################
################################################################################


################################################################################
#NO CENSORING
################################################################################

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

M = 15
Delta = 5
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
#epsilon = 25
#tau = epsilon - (M + Delta + 1)
#xi = min(omega, epsilon - 1)

p = 0.05

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
    write.csv(results, "./results/robust_none_supp.csv")
    
}


################################################################################
#CENSORING
################################################################################

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

M = 15
Delta = 5
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 25
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.05

minU = Delta + 1 
maxU = xi
minV = Delta + 1
maxV = Delta + M

range.1 = c(1:(floor((maxV - minV + 1)/2))) - 1
range.2 = c(1:(floor((maxV - minV + 1)/2) + (maxV - minV + 1) %% 2)) - 1

G = c(sapply(range.1, dbinom, prob = 0.35, size = length(range.1)) * 0.4,
      sapply(range.2, dbinom, prob = 0.35, size = length(range.2)) * 0.6)

THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')

results = matrix(NA, )

cens_time = c(25, 30, 35)
samp_sizes = c(50, 100, 250, 500)
replicates = 1000

results = matrix(NA, nrow = 1, ncol = 7)
colnames(results) = c("e", "n", "p0", "emp_mean",
                      "emp_sd", "thm_sd", "cov_prob")

#calculate true variance

#true variance
cur = c()
for(v in c((Delta + 1):(Delta + M))){
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
for(v in c((Delta + 1):(Delta + M))){
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


for(e in cens_time){
  for(n in samp_sizes){
    
    print(c(e,n))
    
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
    
    mat_row =  c(e, n, p,
                 mean(p_est_vec),
                 sd(p_est_vec),
                 sqrt(true_var / n),
                 sum(cov_ind) / replicates)
    
    results = rbind(results, mat_row)  
    write.csv(results, "./results/robust_cens_supp.csv")
    
  }
}



