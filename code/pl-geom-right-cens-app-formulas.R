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
  
  if( ((Delta + 1) <= v) & (v <= (Delta + m)) ){
    vec_idx = v + 1 - Delta
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

gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

g_tau_MLE = function(v, a, b){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = gnv(v) * (1 - (b/a))^(v - (Delta + 1))
  B = mapply(gnv, c(v_min:v_max))
  C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
  
  return( A * (sum( B * C ))^(-1) )
  
}

thm_formulas = function(obs_data){
  
  a1 = c()
  for(k in c((Delta+1):(Delta+m))){
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
  for(k in c((Delta+1):(Delta+m))){
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

dpsi_dp = function(Yi, Zi, Di, THETA){
  
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