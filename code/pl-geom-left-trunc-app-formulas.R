#from theorem 3.4
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

#theorem formulas LT only
f_X = function(u, THETA){
  
  p = THETA[1]
  
  if( ((Delta + 1) <= u) & (u <= (omega - 1)) ){
    return( p * (1 - p)^((u - (Delta + 1))) )
  }
  if( u == omega ){
    return( (1 - p)^((u - (Delta + 1))) )
  }
  
}

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


DF = function(u, THETA){
  
  num = d2f_dp2(u,THETA) * f_X(u,THETA) - (df_dp(u, THETA))^2
  den = f_X(u,THETA)^2
  
  return(num/den)
  
}

Wi = function(Xi, Yi, u, v){
  
  return(1 * ((Xi == u) & (Yi == v)) )
  
}


psi = function(Xi, Yi, THETA){
  
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