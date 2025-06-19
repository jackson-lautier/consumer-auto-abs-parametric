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

#timing
microbenchmark(constrOptim(init, log_like_fn, NULL, ui=ui, ci=ci),
               c(optimize(P_constraint, c(0,1), tol = 1e-10)$minimum, 
                 mapply(g_tau_hat, c((Delta+1):(m+Delta)), p_hat)),
               thm_formulas(obs_data),
               times = 100)

