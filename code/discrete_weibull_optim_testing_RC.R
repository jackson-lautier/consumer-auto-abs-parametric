p.start = 0.50
p.end = 0.999
seq.length = 100
p.val = seq(p.start, p.end, by = (p.end - p.start)/seq.length)
b.start = 1.00
b.end = 2.00
b.val = seq(b.start, b.end, by = (b.end - b.start)/seq.length)


P_constraint_2 = function(p.param, b.param){
  
  P_input = c(p.param, b.param)
  
  v_min = Delta + 1
  v_max = M + Delta
  
  G_dummy = 1 / rep(v_max - v_min + 1, v_max - v_min + 1)
  THETA_input = c(P_input, G_dummy)
  
  n = nrow(obs_data)
  
  LHS1 = c()
  for(k in c((Delta + 1):(Delta + M))){
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
  for(k in c((Delta + 1):(Delta + M))){
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

persp(p.val, b.val, z, theta = -35, phi = 40, col = "lightblue")
contour(p.val, b.val, z)
filled.contour(p.val, b.val, z,
               xlab = "p", ylab = "b")

which(z == min(z), arr.ind=TRUE)
p.val[which(z == min(z), arr.ind=TRUE)[1]]
b.val[which(z == min(z), arr.ind=TRUE)[2]]
P_constraint_2(p.val[which(z == min(z), arr.ind=TRUE)[1]],
               b.val[which(z == min(z), arr.ind=TRUE)[2]])

p.con.min = 0.001
which(z <= p.con.min, arr.ind=TRUE)

p.val[which(z <= p.con.min, arr.ind=TRUE)[,1]]
b.val[which(z <= p.con.min, arr.ind=TRUE)[,2]]

#p.val = 0.99487
#b.val = 1.4025
#P-constraint = 0.007265869

#25mo loans (aart-2017-25mo)
#p.val = 0.98388
#b.val = 1.198
#P-constraint = 0.001201568

#50mo loans (aart-2017-50mo-2)
#p.val = 0.98514
#b.val = 1.1374
#P-constraint = 0.007265869

#50mo loans (aart-2017-50mo)
#p.val = 0.9874
#b.val = 1.168
#P-constraint = 0.007257382

#37mo loans (aart-2017-37mo)
#p.val = 0.9895
#b.val = 1.271
#P-constraint = 0.001575376

#37mo loans (aart-2017-37mo) - super-prime only
#p.val = 0.98428
#b.val = 1.15
#P-constraint = 0.001521468

#38mo loans (aart-2017-38mo) - super-prime only
#p.val = 0.992
#b.val = 1.336
#P-constraint = 2.259926e-05

#25mo loans (aart-2019-25mo)
#p.val = 0.98258
#b.val = 1.2792
#P-constraint = 0.001995938
