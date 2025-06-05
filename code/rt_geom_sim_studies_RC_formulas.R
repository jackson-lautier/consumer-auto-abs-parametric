#note THETA can be one dimensional 
f_X = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    return( p * (1 - p)^((u - (Delta + 1))) )
  }
  if( u == omega ){
    return( (1 - p)^((u - (Delta + 1))) )
  }
  
}

g_Y = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + M)) ){
    vec_idx = v + 1 - Delta
    return( THETA[vec_idx] )
  }
  
}

alpha = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + M + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + M))), g_Y, THETA = THETA)
    res = append(res,
                 f_X(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

h_star = function(u , v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + M)) &
      (v <= u)){
    
    return(
      (f_X(u, THETA) * g_Y(v, THETA)) / alpha(THETA)
    )
  }
  else{
    return(0)
  }
  
}

h_bar_star = function(u, v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + M)) &
      (v <= u)){
    
    return(
      (S_X(u + 1, THETA) * g_Y(v, THETA)) / alpha(THETA)
    )
  }
  else{
    return(0)
  }
  
}

#note THETA can be one dimensional 
df_dp = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    
    a = 1/p - (u - (Delta + 1))/(1 - p)
    b = f_X(u, THETA)
    
    return(b * a)
    
  }
  if( u == omega ){
    
    a = - (u - (Delta + 1))/(1 - p)
    b = f_X(u, THETA)
    
    return( b * a )
  }
  
}

d2f_dp2 = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    
    a = f_X(u, THETA)
    b = (u - (Delta + 1)) / (1 - p)
    c = -2/p + (u - Delta - 2) / (1 - p)
    
    return(a * b * c)
    
  }
  if( u == omega ){
    
    a = f_X(u, THETA)
    b = (u - (Delta + 1)) * (u - (Delta + 1) - 1)
    c = (1 - p)^(-2)
    
    return( a * b * c )
    
  }
  
}

S_X = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), f_X, THETA)) )
  }
  else{
    return(0)
  }
  
}

dS_dp = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    
    return(
      -(u - Delta - 1) * ((1 - p)^(u - Delta - 2))
    )
    
  }
  else{
    return(0)
  }
  
}

d2S_dp2 = function(u, THETA){

  p = THETA[1]

  if( ((Delta + 1) <= u) & (u <= (omega)) ){

    return(
      (u - Delta - 2) * (u - Delta - 1) * ((1 - p)^(u - Delta - 3))
    )

  }
  else{
    return(0)
  }

}


da_dp = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + M + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + M))), g_Y, THETA = THETA)
    res = append(res,
                 df_dp(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
  
}

d2a_dp2 = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + M + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + M))), g_Y, THETA = THETA)
    res = append(res,
                 d2f_dp2(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
  
}


da_dgv = function(v, THETA){
  
  return( sum( sapply(c(v:omega), f_X, THETA) ) )
  
}

d2a_dpdgv = function(v, THETA){
  
  return( sum( sapply(c(v:omega), df_dp, THETA) ) )
  
}

# u = c( (Delta + 1) : omega)
# reps = length(c( (Delta + 1) : (Delta + M)))
# 
# x_col = c()
# y_col = c()
# for(u in c( (Delta + 1) : omega)){
# 
#   for(v in c( (Delta + 1) : (Delta + M))){
#     if(v <= u){
#       x_col = append(x_col, u)
#       y_col = append(y_col, v)
#     }
#   }
# 
# }
# 
# h_prob = c()
# for(i in c(1:length(x_col))){
#   h_prob = append(h_prob, h_star(x_col[i], y_col[i], THETA))
# }
# 
# h_sum = c()
# for(i in c(1:length(x_col))){
#   h_sum = append(h_sum, sum(h_prob[1:i]))
# }
# 
# l_bound = c(0, h_sum[1:((length(h_prob)-1))])
# u_bound = h_sum
# 
# h_inv = data.frame("X" = x_col,
#                    "Y" = y_col,
#                    "l_bound" = l_bound,
#                    "u_bound" = u_bound)
# 
# h_sim = function(unif){
# 
#   X_i = h_inv$X[(h_inv$l_bound <= unif) & (h_inv$u_bound >= unif)]
#   Y_i = h_inv$Y[(h_inv$l_bound <= unif) & (h_inv$u_bound >= unif)]
# 
#   return(c(X_i, Y_i))
# 
# }

#likelihood function equation (4)
# log_like_fn = function(THETA){
#   
#   alp = alpha(THETA)
#   
#   n = nrow(obs_data)
#   
#   Li = c()
#   for(k in c((Delta + 1):(Delta + M))){
#     for(j in c(k:(omega))){
#       cnt = sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) ))
#       val = ( f_X(j, THETA) * g_Y(k, THETA) ) / alpha(THETA)
#       Li = append(Li, log( val^cnt ) )
#     }
#   }
#   
#   return( -sum(Li) )
#   
# }

log_like_fn = function(THETA){
  
  alp = alpha(THETA)
  
  n = nrow(obs_data)
  
  Li = c()
  for(i in c(1:n)){
    
    if(obs_data$Di[i] == 1){
      contrib = g_Y(obs_data$Yi[i], THETA) * f_X(obs_data$Zi[i], THETA)
      Li = append(Li, log(contrib))
    }
    
    if(obs_data$Di[i] == 0){
      contrib = g_Y(obs_data$Yi[i], THETA) * S_X(obs_data$Zi[i] + 1, THETA)
      Li = append(Li, log(contrib))
    }
  }
  
  val = -n * log(alp) + sum(Li)
  return( -val )
  
}

# h_dot_v = function(v){
#   
#   n = nrow(obs_data)
#   val = c()
#   for(u in c(v:omega)){
#     cur = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
#     val = append(val, cur/n)
#   }
#   
#   return(sum(val))
# }
# 
# h_u_dot = function(u){
#   
#   n = nrow(obs_data)
#   val = c()
#   for(v in c((Delta+1):(min(u,Delta+M)))){
#     cur = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
#     val = append(val, cur/n)
#   }
#   
#   return(sum(val))
# }
# 
# h_u_v = function(u,v){
#   
#   n = nrow(obs_data)
#   val = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
#   return(sum(val) / n)
#   
# }

P_constraint = function(p_input){
  
  n = nrow(obs_data)
  
  LHS = c()
  for(k in c((Delta + 1):(Delta + M))){
    A = gnv(k)
    B = sum(mapply(f_X, c(k:omega), p_input))
    C = sum(mapply(df_dp, c(k:omega), p_input))
    LHS = append(LHS, (A / B) * C )
  }
  
  RHS = c()
  for(i in c(1:n)){
    if(obs_data$Di[i] == 1){
      A = df_dp(obs_data$Zi[i], p_input)
      B = f_X(obs_data$Zi[i], p_input)
    }
    if(obs_data$Di[i] == 0){
      A = dS_dp(obs_data$Zi[i] + 1, p_input)
      B = S_X(obs_data$Zi[i] + 1, p_input)
    }
    
    RHS = append(RHS, A/B)
  
  }
  
  return( ((sum(RHS) / n) - sum(LHS))^2 )
  
}

g_tau_hat = function(v, p_input){
  
  v_min = Delta + 1
  v_max = M + Delta
  
  A = gnv(v) / S_X(v, p_input)
  B = mapply(gnv, c(v_min:v_max))
  C = mapply(S_X, c(v_min:v_max), p_input)
  return( A * (sum( B/C ))^(-1) )
  
}

gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

g_tau_MLE = function(v, a, b){
  
  v_min = Delta + 1
  v_max = M + Delta
  
  A = gnv(v) * (1 - (b/a))^(v - (Delta + 1))
  B = mapply(gnv, c(v_min:v_max))
  C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
  
  return( A * (sum( B * C ))^(-1) )
  
}

thm_formulas = function(obs_data){
  
  a1 = c()
  for(k in c((Delta+1):(Delta+M))){
    a1 = append(a1, (k - (Delta+1)) * gnv(k))
  }
  a2 = c()
  for(i in c(1:nrow(obs_data))){
    a2 = append(a2,
                (obs_data$Zi[i] - (Delta + 1)) * obs_data$Di[i])
  }
  a3 = c()
  for(i in c(1:nrow(obs_data))){
    a3 = append(a3,
                (obs_data$Zi[i] + 1 - (Delta + 1)) * (1 - obs_data$Di[i]))
  }
  a = sum(a1) - (1/nrow(obs_data)) * sum(a2) - (1/nrow(obs_data)) * sum(a3)
  
  b1 = c()
  for(i in c(1:nrow(obs_data))){
    b1 = append(b1, (obs_data$Zi[i] != omega) * obs_data$Di[i])
  }
  b = (1/nrow(obs_data)) * sum(b1)
  
  p_hat = b / (b - a)
  
  G_hat = c()
  for(k in c((Delta+1):(Delta+M))){
    G_hat = append(G_hat, g_tau_MLE(k, a, b) )
  }
  
  return(c(p_hat, G_hat))
  
}

#variance formulas
DF = function(u, THETA){
  
  num = d2f_dp2(u,THETA) * f_X(u,THETA) - (df_dp(u, THETA))^2
  den = f_X(u,THETA)^2
  
  return(num/den)
  
}

Yi_ind = function(Yi, v){
  
  return(1 * (Yi == v) )
  
}


psi = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + M))){
    a = Yi_ind(Yi, v)
    b = sum(sapply(c((v):(omega)), f_X, THETA))
    c = sum(sapply(c((v):(omega)), df_dp, THETA))
    
    lhs = append(lhs, (a/b)*c )
  }
  
  if(Di == 1){
    d = f_X(Zi, THETA)
    e = df_dp(Zi, THETA)
  }
  if(Di == 0){
    d = S_X(Zi + 1, THETA)
    e = dS_dp(Zi + 1, THETA)
  }
    
  rhs = e/d
  
  
  return(sum(lhs) - rhs)
  
}

dpsi_dp = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + M))){
    
    a = Yi_ind(Yi, v)
    
    b1 = sum(sapply(c((v):(omega)), d2f_dp2, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_dp, THETA))
    
    b = (b1 * b2 - b3^2) / (b2^2)
    
    lhs = append(lhs, a * b)
  }
  
  if(Di == 1){
    d = d2f_dp2(Zi, THETA) * f_X(Zi, THETA) - df_dp(Zi, THETA)^2
    e = f_X(Zi, THETA)^2
  }
  if(Di == 0){
    d = d2S_dp2(Zi + 1, THETA) * S_X(Zi + 1, THETA) - (dS_dp(Zi + 1, THETA))^2
    e = S_X(Zi + 1, THETA)^2
  }
  
  rhs = d/e
  
  return(sum(lhs) - rhs)
  
}





