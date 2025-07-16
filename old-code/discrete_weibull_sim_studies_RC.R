require('ggplot2')
require('extrafont') #may need to load fonts
require(latex2exp)
require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

m = 3
Delta = 0
omega = 4

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 6
tau = epsilon - (m + Delta + 1)

p = 0.3
b = 0.5
P = c(p,b)
G = c(0.5, 0.3, 0.2)
THETA = c(P, G)

source('./code/rt_discrete_weibull_sim_studies_RC_formulas.R')

################################################################################
################################################################################
# STATIONARY POINT THEOREM & MLE THEOREM
################################################################################
################################################################################

disp_mat = matrix(NA, nrow = 3, ncol = length(THETA))
rownames(disp_mat) = c("cpu_optim", "thm_3.1", "mle")

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

init = c(0.5, 1)

optim(init,P_constraint,method="L-BFGS-B",
      lower=c(0.001,0.001),
      upper=c(0.999,999))




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
################################################################################
# M-estimator variance theorem
################################################################################
################################################################################

#calculate true variance
V = matrix(NA, 2, 2)

#true variance
cur11 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur11 = append(cur11, 0)
      }
      if(prob > 0){
        cur11 = append(cur11, (psi1(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

V[1,1] = sum(cur11)

cur21 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur21 = append(cur21, 0)
      }
      if(prob > 0){
        cur21 = append(cur21,
                       (psi1(v, u, d, THETA) * psi2(v, u, d, THETA)) * prob)
      }
    }
  }
}

V[2,1] = sum(cur21)
V[1,2] = sum(cur21)

cur22 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur22 = append(cur22, 0)
      }
      if(prob > 0){
        cur22 = append(cur22,
                       (psi2(v, u, d, THETA)^2) * prob)
      }
    }
  }
}

V[2,2] = sum(cur22)

U = matrix(NA, 2, 2)

#true variance
cur11 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur11 = append(cur11, 0)
      }
      if(prob > 0){
        cur11 = append(cur11, (dpsi1_dp(v, u, d, THETA)) * prob)
      }
    }
  }
}

U[1,1] = sum(cur11)

cur21 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur21 = append(cur21, 0)
      }
      if(prob > 0){
        cur21 = append(cur21,
                       (dpsi1_db(v, u, d, THETA)) * prob)
      }
    }
  }
}

U[2,1] = sum(cur21)
U[1,2] = sum(cur21)

cur22 = c()
for(v in c((Delta + 1):(Delta + m))){
  for(u in c(v:omega)){
    for(d in c(0:1)){
      prob1 = d * h_star(u,v,THETA) * 1 * (u <= v + tau)
      prob2 = (1 - d) * h_bar_star(u, v, THETA) * 1 * (v + tau == u)
      prob = prob1 + prob2
      if(prob == 0){
        cur22 = append(cur22, 0)
      }
      if(prob > 0){
        cur22 = append(cur22,
                       (dpsi2_db(v, u, d, THETA)) * prob)
      }
    }
  }
}

U[2,2] = sum(cur22)

true_var = solve(U) %*% V %*% solve(t(U))

#verify with simulation
replicates = 1000
rep_results_p = c()
rep_results_b = c()
n = 250

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
  
  ans = optim(init,P_constraint,method="L-BFGS-B",
              lower=c(0.001,0.001),
              upper=c(0.999,999))
  
  rep_results_p = append(rep_results_p,
                         ans$par[1])
  
  rep_results_b = append(rep_results_b,
                         ans$par[2])
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

mean( sqrt(n) * (rep_results_p - THETA[1]) )
var( sqrt(n) * (rep_results_p - THETA[1]) ); true_var[1,1]
sd( sqrt(n) * (rep_results_p - THETA[1]) ); sqrt(true_var[1,1])

mean( sqrt(n) * (rep_results_b - THETA[2]) )
var( sqrt(n) * (rep_results_b - THETA[2]) ); true_var[2,2]
sd( sqrt(n) * (rep_results_b - THETA[2]) ); sqrt(true_var[2,2])

cov( sqrt(n) * (rep_results_p - THETA[1]),
     sqrt(n) * (rep_results_b - THETA[2]) )

mean(rep_results_p); sd(rep_results_p); sqrt(true_var[1,1] / n)
mean(rep_results_b); sd(rep_results_b); sqrt(true_var[2,2] / n)

#plot results
#plot results
df = data.frame("sim_result" = sqrt(n) * (rep_results_p - THETA[1]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var[1,1]))) +
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

df = data.frame("sim_result" = sqrt(n) * (rep_results_b - THETA[2]))

ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(true_var[2,2]))) +
  #xlab(TeX(("$\\sqrt{ n }(\\hat{p}_n - p_0)$"))) +
  xlab(expression( sqrt( italic(n) ) * ( hat( italic(b) )[italic(n)] - italic(b)[0] ) ) ) +
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

write.csv(facet_data, "fd_rc_binom.csv")
  
################################################################################
################################################################################
# Coverage probabilities
################################################################################
################################################################################

samp_sizes = c(100, 250, 1000)
replicates = 1000

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





