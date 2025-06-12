rm(list=ls())
#problem set-up

M = 15
Delta = 5
omega = 24

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 35
tau = epsilon - (M + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.45

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

#calculate censoring rate
cur = c()
for(v in c((Delta + 1):(Delta + M))){
  for(u in c(v:omega)){
    prob = h_bar_star(u, v, THETA) * 1 * (v + tau == u)
    cur = append(cur, prob)
  }
}
cens.rate = sum(cur)
cens.rate

n = 1000
emp.cens = c()
for(r in c(1:1000)){
  
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
  
  emp.cens = append(emp.cens, (1 - sum(obs_data$Di)/n))
  
  if( (r/100) %in% c(1:10)){
    print(r)
  }
  
}

summary(emp.cens)

