#note THETA can be one dimensional 
f_X = function(u, THETA){
  
  p = THETA[1]
  b = THETA[2]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    return( p^((u - (Delta + 1))^b) - p^((u-Delta)^b) )
  }
  if( u == omega ){
    return( p^((u - (Delta + 1))^b) )
  }
  
}

g_Y = function(v, THETA){
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v + 2 - Delta
    return( THETA[vec_idx] )
  }
  
}

alpha = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y, THETA = THETA)
    res = append(res,
                 f_X(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

h_star = function(u , v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + m)) &
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
      (v <= (Delta + m)) &
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
  b = THETA[2]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    
    a1 = (u - (Delta + 1))^b
    a2 = p^(a1 - 1)
    a3 = (u - Delta)^b
    a4 = p^( a3 - 1)
    
    return(a1 * a2 - a3 * a4)
    
  }
  if( u == omega ){
    
    a1 = (u - (Delta + 1))^b
    a2 = p^(a1 - 1)
    
    return( a1 * a2 )
  }
  
}

d2f_dp2 = function(u, THETA){
  
  p = THETA[1]
  b = THETA[2]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    
    a1 = (u - (Delta + 1))^b
    a2 = p^(a1 - 2)
    a3 = (u - Delta)^b
    a4 = p^( a3 - 2)
    
    return(a1 * (a1 - 1) * a2 - a3 * (a3 - 1) * a4)
    
  }
  if( u == omega ){
    
    a1 = (u - (Delta + 1))^b
    a2 = p^(a1 - 2)
    
    return( a1 * (a1 - 1) * a2 )
  }
  
}

d2f_dpdb = function(u, THETA){
  
  p = THETA[1]
  b = THETA[2]
  
  if( u == Delta + 1){
    
    a1 = -p^((u - Delta)^b - 1)
    a2 = log(p) * log(u - Delta) * ((u - Delta)^(2 * b))
    a3 = (-p^((u - Delta)^b - 1)) * log(u - Delta) * ( (u - Delta)^b )
    
    return( a1 * a2 - a3  )
    
  }
  
  if( ((Delta + 2) <= u) & (u <= (omega - 1)) ){
    
    a1 = -p^((u - Delta)^b - 1)
    a2 = log(p) * log(u - Delta) * ((u - Delta)^(2 * b))
    a3 = (-p^((u - Delta)^b - 1)) * log(u - Delta) * ( (u - Delta)^b )
    
    a4 = u - (Delta + 1)
    
    a5 = (p^(a4^b - 1)) * log(p) * log(a4) * (a4^(2*b))
    a6 = (p^(a4^b - 1)) * log(a4) * (a4^b)
    
    
    return( a1 * a2 - a3 + a5 + a6 )
    
  }
  if( u == omega ){
    
    a4 = u - (Delta + 1)
    
    a5 = (p^(a4^b - 1)) * log(p) * log(a4) * (a4^(2*b))
    a6 = (p^(a4^b - 1)) * log(a4) * (a4^b)
    
    return( a5 + a6 )
  }
  
}

df_db = function(u, THETA){
  
  p = THETA[1]
  b = THETA[2]
  
  if( u == Delta + 1){
    
    return( 0 )
    
  }
  
  if( ((Delta + 1) < u) & (u <= (omega - 1)) ){
    
    a1 = p^((u - (Delta + 1))^b)
    a2 = log(p)
    a3 = (u - (Delta + 1))^b
    a4 = log( u - (Delta + 1) )
    
    a5 = p^((u - Delta)^b)
    a6 = log(p)
    a7 = (u - Delta)^b
    a8 = log( u - Delta )
    
    return(a1 * a2 * a3 * a4  - a5 * a6 * a7 * a8)
    
  }
  if( u == omega ){
    
    a1 = p^((u - (Delta + 1))^b)
    a2 = log(p)
    a3 = (u - (Delta + 1))^b
    a4 = log( u - (Delta + 1) )
    
    return( a1 * a2 * a3 * a4 )
  }
  
}

d2f_db2 = function(u, THETA){
  
  p = THETA[1]
  b = THETA[2]
  
  if( u == Delta + 1){
    
    return( 0 )
    
  }
  
  if( ((Delta + 1) < u) & (u <= (omega - 1)) ){
    
    a1 = p^((u - (Delta + 1))^b)
    a2 = log(p)
    a3 = (u - (Delta + 1))^b
    a4 = log( u - (Delta + 1) )
    
    a5 = p^((u - Delta)^b)
    a6 = log(p)
    a7 = (u - Delta)^b
    a8 = log( u - Delta )
    
    return(a1 * (a2^2) * (a3^2) * (a4^2) 
           + a1 * a2 * (a4^2) * a3
           - a5 * (a6^2) * (a7^2) * (a8^2)
           - a5 * a6 * (a8^2) * a7 )
    
  }
  if( u == omega ){
    
    a1 = p^((u - (Delta + 1))^b)
    a2 = log(p)
    a3 = (u - (Delta + 1))^b
    a4 = log( u - (Delta + 1) )
    
    return( a1 * (a2^2) * (a3^2) * (a4^2) 
            + a1 * a2 * (a4^2) * a3 )
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
    return( sum(sapply(c(u:omega), df_dp, THETA)) )
  }
  else{
    return(0)
  }
  
}

d2S_dp2 = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), d2f_dp2, THETA)) )
  }
  else{
    return(0)
  }
  
}

dS_db = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), df_db, THETA)) )
  }
  else{
    return(0)
  }
  
}

d2S_db2 = function(u, THETA){
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), d2f_db2, THETA)) )
  }
  else{
    return(0)
  }
  
}

d2S_dpdb = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega)) ){
    return( sum(sapply(c(u:omega), d2f_dpdb, THETA)) )
  }
  else{
    return(0)
  }
  
  
}


da_dp = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y, THETA = THETA)
    res = append(res,
                 df_dp(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
  
}

d2a_dp2 = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(omega))){
    v_end = min(u + 1 - Delta, Delta + m + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + m))), g_Y, THETA = THETA)
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
# reps = length(c( (Delta + 1) : (Delta + m)))
# 
# x_col = c()
# y_col = c()
# for(u in c( (Delta + 1) : omega)){
#   
#   for(v in c( (Delta + 1) : (Delta + m))){
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
#   for(v in c((Delta+1):(min(u,Delta+m)))){
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

P_constraint = function(THETA_input){
  
  n = nrow(obs_data)
  
  LHS1 = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = gnv(k)
    B = sum(sapply(c(k:omega), f_X, THETA = THETA_input))
    C = sum(sapply(c(k:omega), df_dp, THETA = THETA_input))
    LHS1 = append(LHS1, (A / B) * C )
  }
  
  RHS1 = c()
  for(i in c(1:n)){
    if(obs_data$Di[i] == 1){
      A = df_dp(obs_data$Zi[i], THETA_input)
      B = f_X(obs_data$Zi[i], THETA_input)
    }
    if(obs_data$Di[i] == 0){
      A = dS_dp(obs_data$Zi[i] + 1, THETA_input)
      B = S_X(obs_data$Zi[i] + 1, THETA_input)
    }
    
    RHS1 = append(RHS1, A/B)
    
  }
  
  LHS2 = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = gnv(k)
    B = sum(sapply(c(k:omega), f_X, THETA = THETA_input))
    C = sum(sapply(c(k:omega), df_db, THETA = THETA_input))
    LHS2 = append(LHS2, (A / B) * C )
  }
  
  RHS2 = c()
  for(i in c(1:n)){
    if(obs_data$Di[i] == 1){
      A = df_db(obs_data$Zi[i], THETA_input)
      B = f_X(obs_data$Zi[i], THETA_input)
    }
    if(obs_data$Di[i] == 0){
      A = dS_db(obs_data$Zi[i] + 1, THETA_input)
      B = S_X(obs_data$Zi[i] + 1, THETA_input)
    }
    
    RHS2 = append(RHS2, A/B)
    
  }
  
  return( ((sum(RHS1) / n) - sum(LHS1))^2 +
            ((sum(RHS2) / n) - sum(LHS2))^2)
  
}

g_tau_hat = function(v, p_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = gnv(v) / S_X(v, p_input)
  B = mapply(gnv, c(v_min:v_max))
  C = sapply(c(v_min:v_max), S_X, THETA = p_input)
  #C = mapply(S_X, c(v_min:v_max), p_input)
  return( A * (sum( B/C ))^(-1) )
  
}

gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

# g_tau_MLE = function(v){
#   
#   v_min = Delta + 1
#   v_max = m + Delta
#   
#   A = h_dot_v(v) * (1 - (b/a))^(v - (Delta + 1))
#   B = mapply(h_dot_v, c(v_min:v_max))
#   C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
#   
#   return( A * (sum( B * C ))^(-1) )
#   
# }

# thm_formulas = function(obs_data){
#   
#   a1 = c()
#   for(k in c((Delta+1):(Delta+m))){
#     a1 = append(a1, (k - (Delta+1)) * h_dot_v(k))
#   }
#   a2 = c()
#   for(j in c((Delta + 1):omega)){
#     a2 = append(a2, (j - (Delta + 1)) * h_u_dot(j))
#   }
#   b1 = c()
#   for(j in c((Delta + 1):(omega - 1))){
#     b1 = append(b1, h_u_dot(j))
#   }
#   a = sum(a1) - sum(a2)
#   b = sum(b1)
#   
#   p_hat = b / (b - a)
#   
#   G_hat = c()
#   for(k in c((Delta+1):(Delta+m))){
#     A = h_dot_v(k) * (1 - (b/a))^(k - (Delta + 1))
#     B = mapply(h_dot_v, c((Delta+1):(Delta+m)))
#     C = (1 - (b/a))^(c((Delta+1):(Delta+m)) - (Delta + 1))
#     G_hat = append(G_hat, A * (sum( B * C ))^(-1) )
#   }
#   
#   return(c(p_hat, G_hat))
#   
# }

#variance formulas
Yi_ind = function(Yi, v){
  
  return(1 * (Yi == v) )
  
}


psi1 = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + m))){
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

dpsi1_dp = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + m))){
    
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

psi2 = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + m))){
    a = Yi_ind(Yi, v)
    b = sum(sapply(c((v):(omega)), f_X, THETA))
    c = sum(sapply(c((v):(omega)), df_db, THETA))
    
    lhs = append(lhs, (a/b)*c )
  }
  
  if(Di == 1){
    d = f_X(Zi, THETA)
    e = df_db(Zi, THETA)
  }
  if(Di == 0){
    d = S_X(Zi + 1, THETA)
    e = dS_db(Zi + 1, THETA)
  }
  
  rhs = e/d
  
  
  return(sum(lhs) - rhs)
  
}

dpsi2_db = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + m))){
    
    a = Yi_ind(Yi, v)
    
    b1 = sum(sapply(c((v):(omega)), d2f_db2, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_db, THETA))
    
    b = (b1 * b2 - b3^2) / (b2^2)
    
    lhs = append(lhs, a * b)
  }
  
  if(Di == 1){
    d = d2f_db2(Zi, THETA) * f_X(Zi, THETA) - df_db(Zi, THETA)^2
    e = f_X(Zi, THETA)^2
  }
  if(Di == 0){
    d = d2S_db2(Zi + 1, THETA) * S_X(Zi + 1, THETA) - (dS_db(Zi + 1, THETA))^2
    e = S_X(Zi + 1, THETA)^2
  }
  
  rhs = d/e
  
  return(sum(lhs) - rhs)
  
}

dpsi1_db = function(Yi, Zi, Di, THETA){
  
  lhs = c()
  rhs = c()
  
  for(v in c((Delta + 1):(Delta + m))){
    
    a = Yi_ind(Yi, v)
    
    b1 = sum(sapply(c((v):(omega)), d2f_dpdb, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_dp, THETA))
    b4 = sum(sapply(c((v):(omega)), df_db, THETA))
    
    b = (b1 * b2 - b3 * b4) / (b2^2)
    
    lhs = append(lhs, a * b)
  }
  
  if(Di == 1){
    d = d2f_dpdb(Zi, THETA) * f_X(Zi, THETA) -
      df_dp(Zi, THETA) * df_db(Zi, THETA)
    e = f_X(Zi, THETA)^2
  }
  if(Di == 0){
    d = d2S_dpdb(Zi + 1, THETA) * S_X(Zi + 1, THETA) -
      dS_dp(Zi + 1, THETA) * dS_db(Zi + 1, THETA)
    e = S_X(Zi + 1, THETA)^2
  }
  
  rhs = d/e
  
  return(sum(lhs) - rhs)
  
}



