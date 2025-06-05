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
  for(k in c((Delta + 1):(Delta + M))){
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
  
  if( ((Delta + 1) <= v) & (v <= (Delta + M)) ){
    vec_idx = (v - Delta) + (maxU - minU) + 1
    return( THETA[vec_idx] )
  }
  
}

alpha.1 = function(THETA){
  
  res = c()
  for(u in c((Delta + 1):(xi))){
    v_end = min(u + 1 - Delta, Delta + M + 1)
    g_sum = sapply(c((Delta + 1):(min(u, Delta + M))), g_Y.1, THETA = THETA)
    res = append(res,
                 f_X.1(u, THETA) * sum(g_sum))
    
  }
  
  return(sum(res))
  
}

h_star.1 = function(u , v, THETA){
  
  if( ((Delta + 1) <= u) &
      (u <= omega) &
      ((Delta + 1) <= v) &
      (v <= (Delta + M)) &
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
      (v <= (Delta + M)) &
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
  }else{ l2 = 1 }
  
  a1 = ifelse(l1 == 0, 0, log(l1))
  a2 = ifelse(l2 == 0, 0, log(l2))
  
  return( sum(a1) + sum(a2) )
  
}
