#p.val = seq(0.965, 0.97, by = 5e-05)
#b.val = seq(0.910, 0.925, by = 0.00015)

p.val = seq(0.900, 0.999, by = 0.00099)
b.val = seq(0.750, 1.500, by = 0.0075)

P_constraint_2 = function(p.param, b.param){
  
  P_input = c(p.param, b.param)
  
  v_min = Delta + 1
  v_max = m + Delta
  
  G_dummy = 1 / rep(v_max - v_min + 1, v_max - v_min + 1)
  THETA_input = c(P_input, G_dummy)
  
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
  
  return( min(0.3,
              ((sum(RHS1) / n) - sum(LHS1))^2 +
            ((sum(RHS2) / n) - sum(LHS2))^2)  )
  
}

P_constraint_3 <- Vectorize(P_constraint_2)

z <- outer(p.val,b.val,P_constraint_3);

persp(p.val, b.val, z, theta = -25, phi = 30, col = "lightblue")
contour(p.val, b.val, z)
filled.contour(p.val, b.val, z,
               xlab = "p", ylab = "b")

which(z == min(z), arr.ind=TRUE)
p.val[which(z == min(z), arr.ind=TRUE)[1]]
b.val[which(z == min(z), arr.ind=TRUE)[2]]
P_constraint_2(0.99108, 1.245)

which(z <= 0.0299, arr.ind=TRUE)
