rm(list=ls())
#problem set-up

m = 3
Delta = 0
omega = 4

#min = m + Delta + 1
#max = m + omega (set to this value for no censoring)
epsilon = 6
tau = epsilon - (m + Delta + 1)
xi = min(omega, epsilon - 1)

p = 0.3
G = c(0.5, 0.3, 0.2)
THETA = c(p, G)

minU = Delta + 1
maxU = min(omega, epsilon - 1)
minV = Delta + 1
maxV = Delta + m

source('./code/rt_geom_sim_studies_RC_formulas.R')

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



theta_hat = thm_formulas(obs_data)
l.0 = -log_like_fn(theta_hat)

#haz estimator functions from IME paper
Unx <- function(x) {
  ind_x <- ifelse( (obs_data$Yi <= x) & (x <= obs_data$Zi),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}

fnx <- function(x) {
  ind_x <- ifelse( (obs_data$Zi == x) & (obs_data$Di == 1),1,0 )
  return( (1/nrow(obs_data)) * (sum(ind_x)) )
}

lnx <- function(x) {
  
  return( fnx(x) / Unx(x) )
  
}


gnx <- function(x) {
  return(sum(obs_data$Yi == x) / nrow(obs_data))
}

f_est = function(u){
  
  idx = u - Delta
  
  if( u == (Delta + 1)){
    return( haz_est[idx] )
  }
  if( ((Delta + 1) < u) & (u <= xi ) ){
    return(haz_est[idx] * prod(1 - haz_est[1:(idx - 1)]))
  }
  
}

g_est = function(g){
  
  A = (sum(obs_data$Yi == g)/nrow(obs_data)) / 
    sum(sapply(c(g:xi), f_est))
  
  B = c()
  for(k in c((Delta + 1):(Delta + m))){
    B = append(B,
               (sum(obs_data$Yi == k)/nrow(obs_data)) / 
                 sum(sapply(c(k:xi), f_est)))
  }
  
  return( A / (sum(B)) )
  
}

f_X.1 = function(u, THETA){
  
  U = THETA[(minU - Delta):(maxU - Delta)]
  
  vec_idx = u - Delta
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( U[vec_idx] )
  }
  else{
    return(0)
  }
  
}


g_Y.1 = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = (v - Delta) + (maxU - minU) + 1
    return( THETA[vec_idx] )
  }
  
}

alpha.1 = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(xi))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y.1, THETA = THETA)
    res = append(res,
                 f_X.1(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

h_star.1 = function(u , v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
      (v <= u)){
    
    return(
      (f_X.1(u, THETA) * g_Y.1(v, THETA)) / alpha.1(THETA)
    )
  }
  else{
    return(0)
  }
  
}

S_X.1 = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (xi)) ){
    return( sum(sapply(c(u:xi), f_X.1, THETA)) )
  }
  else{
    return(0)
  }
  
}

h_bar_star.1 = function(u, v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
      (v <= u)){
    
    return(
      (S_X.1(u + 1, THETA) * g_Y.1(v, THETA)) / alpha.1(THETA)
    )
  }
  else{
    return(0)
  }
  
}


log_like_fn_0 = function(THETA){
  
  dat.1 = obs_data[obs_data$Di == 1,]
  dat.2 = obs_data[obs_data$Di == 0,]
  
  l1 = c()
  for(i in c(1:nrow(dat.1))){
    l1 = append(l1, h_star.1(u = dat.1$Zi[i],
                           v = dat.1$Yi[i],
                           THETA = THETA))
  }
  
  if(nrow(dat.2) > 0){
  l2 = c()
  for(i in c(1:nrow(dat.2))){
    l2 = append(l2, h_bar_star.1(u = dat.2$Zi[i],
                               v = dat.2$Yi[i],
                               THETA = THETA))
  }
  #else{ l2 = 1 }
  }
  
  return( sum(log(l1)) + sum(log(l2)) )
  
}


#null: h_* is true
haz_est = sapply(c( (Delta + 1) : xi), lnx)
U_est = sapply(c( (Delta + 1) : xi), f_est)
G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
l.1 = log_like_fn_0(c(U_est, G_est))

Q = -2 * (l.0 - l.1)


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
  G_est = sapply(c((Delta + 1):(Delta + m)), g_est)
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
#sum(df >= qchisq(0.95, deg_free)) / replicates

#write.csv(df, "sim-result-rc.csv")
require('ggplot2')
ggplot(df, aes(x=sim_result)) + 
  geom_density(color = "blue", linetype = "dashed") +
  stat_function(fun = dchisq, args = list(df = deg_free))



ggsave("LRT-example-RC.pdf",height=4,width=6,device = cairo_pdf)





