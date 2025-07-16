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
#LEFT-TRUNCATION (S3) & PL Geometric
################################################################################

#require('ggplot2')
#require('extrafont') #may need to load fonts
#require(latex2exp)
#require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

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

#require('ggplot2')
#require('extrafont') #may need to load fonts
#require(latex2exp)
#require(microbenchmark) #used for timing purposes

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
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')

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

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

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

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

m = 3
Delta = 0
omega = 4

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
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
# Figure 1: Asymptotic Normality Verification
################################################################################
################################################################################


################################################################################
#LEFT-TRUNCATION (S3) & PL Geometric
################################################################################

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

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

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

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
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

source('./code/rt_geom_sim_studies_RC_formulas.R')

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

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

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

# require('ggplot2')
# require('extrafont') #may need to load fonts
# require(latex2exp)
# require(microbenchmark) #used for timing purposes

rm(list=ls())
#problem set-up

m = 3
Delta = 0
omega = 4

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
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

# library('ggplot2')
# require('extrafont')

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
# Table 2: Robustness Simulation Study
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

m = 20
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
    write.csv(results, "./results/robust_none.csv")
    
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

m = 20
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

results = matrix(NA, )

cens_time = c(26, 32, 38)
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
    write.csv(results, "./results/robust_cens.csv")
    
  }
}



