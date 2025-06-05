u = c( (Delta + 1) : omega)
reps = length(c( (Delta + 1) : (Delta + m)))

x_col = c()
y_col = c()
for(u in c( (Delta + 1) : omega)){
  
  for(v in c( (Delta + 1) : (Delta + m))){
    if(v <= u){
      x_col = append(x_col, u)
      y_col = append(y_col, v)
    }
  }
  
}

h_prob = c()
for(i in c(1:length(x_col))){
  h_prob = append(h_prob, h_star(x_col[i], y_col[i], THETA))
}

h_sum = c()
for(i in c(1:length(x_col))){
  h_sum = append(h_sum, sum(h_prob[1:i]))
}

l_bound = c(0, h_sum[1:((length(h_prob)-1))])
u_bound = h_sum

h_inv = data.frame("X" = x_col,
                   "Y" = y_col,
                   "l_bound" = l_bound,
                   "u_bound" = u_bound)

h_sim = function(unif){
  
  X_i = h_inv$X[(h_inv$l_bound <= unif) & (h_inv$u_bound >= unif)]
  Y_i = h_inv$Y[(h_inv$l_bound <= unif) & (h_inv$u_bound >= unif)]
  
  return(c(X_i, Y_i))
  
}