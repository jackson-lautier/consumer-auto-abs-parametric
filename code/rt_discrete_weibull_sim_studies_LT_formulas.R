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


#likelihood function equation (4)
log_like_fn = function(THETA){
  
  alp = alpha(THETA)
  
  n = nrow(obs_data)
  
  Li = c()
  for(k in c((Delta + 1):(Delta + m))){
    for(j in c(k:(omega))){
      cnt = sum(( (obs_data$Yi == k ) & (obs_data$Xi == j) ))
      val = ( f_X(j, THETA) * g_Y(k, THETA) ) / alpha(THETA)
      Li = append(Li, cnt * log( val ) )
    }
  }
  
  return( -sum(Li) )
  
}

h_dot_v = function(v){
  
  n = nrow(obs_data)
  val = c()
  for(u in c(v:omega)){
    cur = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
    val = append(val, cur/n)
  }
  
  return(sum(val))
}

h_u_dot = function(u){
  
  n = nrow(obs_data)
  val = c()
  for(v in c((Delta+1):(min(u,Delta+m)))){
    cur = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
    val = append(val, cur/n)
  }
  
  return(sum(val))
}

h_u_v = function(u,v){
  
  n = nrow(obs_data)
  val = sum( (obs_data$Yi == v) & (obs_data$Xi == u) )
  return(sum(val) / n)
  
}

P_constraint = function(P_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  G_dummy = 1 / rep(v_max - v_min + 1, v_max - v_min + 1)
  THETA = c(P_input, G_dummy)
  
  LHS1 = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = h_dot_v(k)
    B = sum(sapply(c(k:omega), f_X, THETA))
    C = sum(sapply(c(k:omega), df_dp, THETA))
    LHS1 = append(LHS1, (A / B) * C )
  }
  
  RHS1 = c()
  for(k in c((Delta + 1):(Delta + m))){
    for(j in c(k:omega)){
      A = h_u_v(j,k)
      B = f_X(j, THETA)
      C = df_dp(j, THETA)
      RHS1 = append(RHS1, (A / B) * C)
    }
  }
  
  LHS2 = c()
  for(k in c((Delta + 1):(Delta + m))){
    A = h_dot_v(k)
    B = sum(sapply(c(k:omega), f_X, THETA))
    C = sum(sapply(c(k:omega), df_db, THETA))
    LHS2 = append(LHS2, (A / B) * C )
  }
  
  RHS2 = c()
  for(k in c((Delta + 1):(Delta + m))){
    for(j in c(k:omega)){
      A = h_u_v(j,k)
      B = f_X(j, THETA)
      C = df_db(j, THETA)
      RHS2 = append(RHS2, (A / B) * C)
    }
  }
  
  a1 = (sum(RHS1) - sum(LHS1))^2
  a2 = (sum(RHS2) - sum(LHS2))^2

  return( a1 + a2 )
  
}

g_tau_hat = function(v, p_input){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = h_dot_v(v) / S_X(v, p_input)
  B = mapply(h_dot_v, c(v_min:v_max))
  C = sapply(c(v_min:v_max), S_X, THETA = p_input)
  #C = mapply(S_X, c(v_min:v_max), p_input)
  return( A * (sum( B/C ))^(-1) )
  
}

gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

g_tau_MLE = function(v){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = h_dot_v(v) * (1 - (b/a))^(v - (Delta + 1))
  B = mapply(h_dot_v, c(v_min:v_max))
  C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
  
  return( A * (sum( B * C ))^(-1) )
  
}


#variance formulas
Wi = function(Xi, Yi, u, v){

  return(1 * ((Xi == u) & (Yi == v)) )

}


psi1 = function(Xi, Yi, THETA){

  lhs = c()
  rhs = c()
  for(v in c((Delta + 1):(Delta + m))){
    a = sum(sapply(c((v):(omega)), Wi, Xi = Xi, Yi = Yi, v = v))
    b = sum(sapply(c((v):(omega)), f_X, THETA))
    c = sum(sapply(c((v):(omega)), df_dp, THETA))

    cur = c()
    for(u in c(v:omega)){
      cur = append(cur,
                   (Wi(Xi, Yi, u, v) / f_X(u, THETA)) * df_dp(u, THETA))
    }

    lhs = append(lhs, (a/b)*c )
    rhs = append(rhs, sum(cur))
  }

  return(sum(lhs) - sum(rhs))

}

psi2 = function(Xi, Yi, THETA){
  
  lhs = c()
  rhs = c()
  for(v in c((Delta + 1):(Delta + m))){
    a = sum(sapply(c((v):(omega)), Wi, Xi = Xi, Yi = Yi, v = v))
    b = sum(sapply(c((v):(omega)), f_X, THETA))
    c = sum(sapply(c((v):(omega)), df_db, THETA))
    
    cur = c()
    for(u in c(v:omega)){
      cur = append(cur,
                   (Wi(Xi, Yi, u, v) / f_X(u, THETA)) * df_db(u, THETA))
    }
    
    lhs = append(lhs, (a/b)*c )
    rhs = append(rhs, sum(cur))
  }
  
  return(sum(lhs) - sum(rhs))
  
}

DF11 = function(u, THETA){
  
  num = d2f_dp2(u,THETA) * f_X(u,THETA) - (df_dp(u, THETA))^2
  den = f_X(u,THETA)^2
  
  return(num/den)
  
}

dpsi1_dp = function(Xi, Yi, THETA){

  lhs = c()
  rhs = c()
  for(v in c((Delta + 1):(Delta + m))){

    a = sum(sapply(c((v):(omega)), Wi, Xi = Xi, Yi = Yi, v = v))

    b1 = sum(sapply(c((v):(omega)), d2f_dp2, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_dp, THETA))

    b = (b1 * b2 - b3^2) / (b2^2)

    cur = c()
    for(u in c(v:omega)){
      cur = append(cur,
                   (Wi(Xi, Yi, u, v) * DF11(u, THETA)))
    }

    lhs = append(lhs, a*b )
    rhs = append(rhs, sum(cur))
  }

  return(sum(lhs) - sum(rhs))

}

DF12 = function(u, THETA){
  
  num = d2f_dpdb(u,THETA) * f_X(u,THETA) - df_dp(u, THETA) * df_db(u,THETA)
  den = f_X(u,THETA)^2
  
  return(num/den)
  
}

dpsi1_db = function(Xi, Yi, THETA){
  
  lhs = c()
  rhs = c()
  for(v in c((Delta + 1):(Delta + m))){
    
    a = sum(sapply(c((v):(omega)), Wi, Xi = Xi, Yi = Yi, v = v))
    
    b1 = sum(sapply(c((v):(omega)), d2f_dpdb, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_dp, THETA))
    b4 = sum(sapply(c((v):(omega)), df_db, THETA))
    
    b = (b1 * b2 - b3 * b4) / (b2^2)
    
    cur = c()
    for(u in c(v:omega)){
      cur = append(cur,
                   (Wi(Xi, Yi, u, v) * DF12(u, THETA)))
    }
    
    lhs = append(lhs, a*b )
    rhs = append(rhs, sum(cur))
  }
  
  return(sum(lhs) - sum(rhs))
  
}

DF21 = function(u, THETA){
  
  return( DF12(u, THETA) )
  
}

dpsi2_dp = function(Xi, Yi, THETA){
  
  return( dpsi1_db(Xi, Yi, THETA) )
  
}

DF22 = function(u, THETA){
  
  num = d2f_db2(u,THETA) * f_X(u,THETA) - (df_db(u, THETA))^2
  den = f_X(u,THETA)^2
  
  return(num/den)
  
}

dpsi2_db = function(Xi, Yi, THETA){
  
  lhs = c()
  rhs = c()
  for(v in c((Delta + 1):(Delta + m))){
    
    a = sum(sapply(c((v):(omega)), Wi, Xi = Xi, Yi = Yi, v = v))
    
    b1 = sum(sapply(c((v):(omega)), d2f_db2, THETA))
    b2 = sum(sapply(c((v):(omega)), f_X, THETA))
    b3 = sum(sapply(c((v):(omega)), df_db, THETA))
    
    b = (b1 * b2 - b3^2) / (b2^2)
    
    cur = c()
    for(u in c(v:omega)){
      cur = append(cur,
                   (Wi(Xi, Yi, u, v) * DF22(u, THETA)))
    }
    
    lhs = append(lhs, a*b )
    rhs = append(rhs, sum(cur))
  }
  
  return(sum(lhs) - sum(rhs))
  
}


