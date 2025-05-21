f_star <- function(x) {
  res = sum( (obs_data$C == 0) * (obs_data$Z == x))
  return(res/n)
}

est_haz <- function(x) {
  num = sum( (obs_data$C == 0) * (obs_data$Z == x))
  den = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(num/den)
}

est_C <- function(x) {
  ans = sum( (x >= obs_data$Y) * (x <= obs_data$Z))
  return(ans/n)
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}