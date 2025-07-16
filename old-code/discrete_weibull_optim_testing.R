p.val = seq(0.949, 0.959, by = 0.0001)
b.val = seq(0.88, 0.96, by = 8e-04)

pairs <- expand.grid(p.val, b.val)

res = c()
for(i in c(1:nrow(pairs))){
  
  param = as.numeric(pairs[i,])
  
  res = append(res, P_constraint(param))
  
  if( (i/1000) %in% c(1:10)){
    print(i)
  }
  
}

plot(c(1:nrow(pairs)), res, type = "l")

x <- seq(-10, 10, by = 1)
y <- seq(-10, 10, by = 1)

f <- function(x,y){return(-x^2 - y^2)}

P_constraint_2 = function(p.param, b.param){
  
  P_input = c(p.param, b.param)
  
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
  
  return( min(0.3, a1 + a2) )
  
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
P_constraint_2(0.959, 0.9208)

which(z <= 0.0043, arr.ind=TRUE)
