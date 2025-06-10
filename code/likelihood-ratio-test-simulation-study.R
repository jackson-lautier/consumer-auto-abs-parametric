require('ggplot2')
require('extrafont')

rm(list=ls())
#problem set-up

M = 5
Delta = 0
omega = 8

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 10
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.3
G = c(0.35, 0.25, 0.20, 0.15, 0.05)
THETA = c(p, G)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + M

source("./code/ime-formulas.R")
source('./code/rt_geom_sim_studies_RC_formulas.R')
source('./code/hsim-formulas.R')

#########################################################################
#LRT simulation study
n = 1000
replicates = 1000

results = c()
for(r in c(1:replicates)){
  
  set.seed(r)
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
  
  #ime: no restrictions
  haz_est = sapply(c( (Delta + 1) : xi), lnx)
  U_est = sapply(c( (Delta + 1) : xi), f_est)
  G_est = sapply(c((Delta + 1):(Delta + M)), g_est)
  l.1 = log_like_fn_0(c(U_est, G_est))
  
  #geometric is true: null hypothesis
  theta_hat = thm_formulas(obs_data)
  l.0 = -log_like_fn(theta_hat)
  
  results = append(results, -2 * (l.0 - l.1))
  
  if( (r/100) %in% c(1:(replicates/100))){
    print(r)
  }
  
}

deg_free = xi - (Delta + 1) - 1

df = data.frame("sim_result" = results)

sum(df$sim_result >= qchisq(0.95, deg_free)) / replicates


require('ggplot2')
ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dchisq, args = list(df = deg_free)) +
  xlab(expression( 'LRT Statistic:'~ Lambda[italic(n)] ) ) +
  ylab("Density Height") +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman", face = "italic"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.position = "bottom")

ggsave("LRT-simulation-study.pdf",height=4,width=6,device = cairo_pdf)





