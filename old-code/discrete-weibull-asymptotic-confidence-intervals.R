################################################################################
# LEFT-TRUNCATION
################################################################################
rm(list=ls())

m = 3; M = 3
Delta = 0
omega = 4

p = 0.3
b = 0.5
P = c(p,b)
G = c(0.5, 0.3, 0.2)
THETA = c(P, G)

source('./code/rt_discrete_weibull_sim_studies_LT_formulas.R')
source('./code/hsim-formulas.R')

n = 500
replicates = 1000

results = matrix(NA, nrow = replicates, ncol = length(c((Delta + 1):(omega))))

for(r in c(1:replicates)){
  
  #sim data
  set.seed(r)
  sample = sapply(runif(n), h_sim)
  
  #observed data
  Yi = sample[2,]
  Xi = sample[1,]
  
  obs_data = data.frame(
    "Yi" = Yi,
    "Xi" = Xi
  )
  
  init = c(0.5, 1)
  #point estimates
  THETA_est = optim(init,P_constraint,method="L-BFGS-B",
                    lower=c(0.001,0.001),
                    upper=c(0.999,999))
  
  #variance matrix
  V = matrix(NA, 2, 2)
  
  Vi = c()
  for(i in c(1:n)){
    
    Xi = obs_data$Xi[i]
    Yi = obs_data$Yi[i]
    
    Vi = append(Vi, (psi1(Xi, Yi, THETA_est$par))^2 )
    
  }
  
  V[1,1] = (1/n) * sum(Vi)
  
  Vi = c()
  for(i in c(1:n)){
    
    Xi = obs_data$Xi[i]
    Yi = obs_data$Yi[i]
    
    Vi = append(Vi,
                psi1(Xi, Yi, THETA_est$par) * psi2(Xi, Yi, THETA_est$par) )
    
  }
  
  V[2,1] = (1/n) * sum(Vi)
  V[1,2] = V[2,1]
  
  #Vn estimate b
  Vi = c()
  for(i in c(1:n)){
    
    Xi = obs_data$Xi[i]
    Yi = obs_data$Yi[i]
    
    Vi = append(Vi, (psi2(Xi, Yi, THETA_est$par))^2 )
    
  }
  
  V[2,2] = (1/n) * sum(Vi)
  
  
  V_est = solve(V)
  
  dw_lam_hat = sapply(c((Delta + 1):(omega)), lambda, THETA = THETA_est$par)
  
  sig_est = c()
  for(k in c((Delta + 1):(omega))){
    
    G.dat = c( dlambda_dp(k, THETA_est$par), dlambda_db(k, THETA_est$par) )
    G_est = matrix(G.dat, nrow = 1, ncol = 2)
    sig_est = append(sig_est, G_est %*% V_est %*% t(G_est))
    
  }
  
  CI_low = dw_lam_hat - qnorm(0.975) * sqrt(sig_est / n)
  CI_upp = dw_lam_hat + qnorm(0.975) * sqrt(sig_est / n)
  
  res.row = 1 *
    ( (CI_low <= sapply(c(1:4), lambda, THETA = THETA)) &
        (CI_upp >= sapply(c(1:4), lambda, THETA = THETA)) )
  
  results[r, ] = res.row
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

colSums(results) / 1000

################################################################################
# LEFT-TRUNCATION + RIGHT-CENSORING
################################################################################
rm(list=ls())

m = 3; M = 3
Delta = 0
omega = 4

epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.3
b = 0.5
P = c(p,b)
G = c(0.5, 0.3, 0.2)
THETA = c(P, G)

source('./code/rt_discrete_weibull_sim_studies_RC_formulas.R')
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

n = 500
replicates = 1000

results = matrix(NA, nrow = replicates, ncol = length(c((Delta + 1):(omega))))

for(r in c(1:replicates)){
  
  #sim data
  set.seed(r + replicates)
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
  
  init = c(0.5, 1)
  #point estimates
  THETA_est = optim(init,P_constraint,method="L-BFGS-B",
                    lower=c(0.001,0.001),
                    upper=c(0.999,999))
  
  #Vn estimate
  V = matrix(NA, 2, 2)
  
  Vi = c()
  for(i in c(1:n)){
    
    Zi = obs_data$Zi[i]
    Yi = obs_data$Yi[i]
    Di = obs_data$Di[i]
    
    Vi = append(Vi, (psi1(Yi, Zi, Di, THETA_est$par))^2 )
    
  }
  
  V[1,1] = (1/n) * sum(Vi)
  
  Vi = c()
  for(i in c(1:n)){
    
    Zi = obs_data$Zi[i]
    Yi = obs_data$Yi[i]
    Di = obs_data$Di[i]
    
    Vi = append(Vi,
                psi1(Yi, Zi, Di, THETA_est$par) * psi2(Yi, Zi, Di, THETA_est$par) )
    
  }
  
  V[2,1] = (1/n) * sum(Vi)
  V[1,2] = V[2,1]
  
  #Vn estimate b
  Vi = c()
  for(i in c(1:n)){
    
    Zi = obs_data$Zi[i]
    Yi = obs_data$Yi[i]
    Di = obs_data$Di[i]
    
    Vi = append(Vi, (psi2(Yi, Zi, Di, THETA_est$par))^2 )
    
  }
  
  V[2,2] = (1/n) * sum(Vi)
  
  
  V_est = solve(V)
  
  dw_lam_hat = sapply(c((Delta + 1):(omega)), lambda, THETA = THETA_est$par)
  
  sig_est = c()
  for(k in c((Delta + 1):(omega))){
    
    G.dat = c( dlambda_dp(k, THETA_est$par), dlambda_db(k, THETA_est$par) )
    G_est = matrix(G.dat, nrow = 1, ncol = 2)
    sig_est = append(sig_est, G_est %*% V_est %*% t(G_est))
    
  }
  
  CI_low = dw_lam_hat - qnorm(0.975) * sqrt(sig_est / n)
  CI_upp = dw_lam_hat + qnorm(0.975) * sqrt(sig_est / n)
  
  res.row = 1 *
    ( (CI_low <= sapply(c(1:4), lambda, THETA = THETA)) &
        (CI_upp >= sapply(c(1:4), lambda, THETA = THETA)) )
  
  results[r, ] = res.row
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

colSums(results) / 1000
















samp_sizes = c(100, 250, 1000)
replicates = 1000

init = c(0.5, 1)

cov_prob = matrix(NA, nrow = 2, ncol = length(samp_sizes))
colnames(cov_prob) = paste("n", samp_sizes, sep="")
rownames(cov_prob) = c("p", "beta")

for(n in samp_sizes){
  
  print(paste("n", samp_sizes, sep=""))
  
  cov_ind_p = c()
  cov_ind_b = c()
  for(r in c(1:replicates)){
    
    #sim data
    k = (which(samp_sizes == n) - 1) * replicates
    set.seed(k + r)
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
    
    #point estimates
    THETA_est = optim(init,P_constraint,method="L-BFGS-B",
                      lower=c(0.001,0.001),
                      upper=c(0.999,999))
    
    #Vn estimate
    V = matrix(NA, 2, 2)
    
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi1(Yi, Zi, Di, THETA_est$par))^2 )
      
    }
    
    V[1,1] = (1/n) * sum(Vi)
    
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi,
                  psi1(Yi, Zi, Di, THETA_est$par) * psi2(Yi, Zi, Di, THETA_est$par) )
      
    }
    
    V[2,1] = (1/n) * sum(Vi)
    V[1,2] = V[2,1]
    
    #Vn estimate b
    Vi = c()
    for(i in c(1:n)){
      
      Zi = obs_data$Zi[i]
      Yi = obs_data$Yi[i]
      Di = obs_data$Di[i]
      
      Vi = append(Vi, (psi2(Yi, Zi, Di, THETA_est$par))^2 )
      
    }
    
    V[2,2] = (1/n) * sum(Vi)
    
    
    V_est = solve(V)
    
    sigma1 = sqrt( V_est[1,1] )
    sigma2 = sqrt( V_est[2,2] )
    
    #CI estimates
    CI_low = THETA_est$par[1] - qnorm(0.975) * sigma1 / sqrt(n)
    CI_upp = THETA_est$par[1] + qnorm(0.975) * sigma1 / sqrt(n)
    if( CI_low < 0 ){print("check")}
    if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[1]) & (THETA[1] <= CI_upp) )
    
    cov_ind_p = append(cov_ind_p, res)
    
    #CI estimates
    CI_low = THETA_est$par[2] - qnorm(0.975) * sigma2 / sqrt(n)
    CI_upp = THETA_est$par[2] + qnorm(0.975) * sigma2 / sqrt(n)
    if( CI_low < 0 ){print("check")}
    #if( CI_upp > 1){print("check")}
    res = 1 * ( (CI_low <= THETA[2]) & (THETA[2] <= CI_upp) )
    
    cov_ind_b = append(cov_ind_b, res)
    
    if( (r/100) %in% c(1:10)){
      print(r)
    }
    
  }
  
  cov_prob[1,paste("n", n, sep="")] = sum(cov_ind_p) / replicates
  cov_prob[2,paste("n", n, sep="")] = sum(cov_ind_b) / replicates
  
  
}

cov_prob