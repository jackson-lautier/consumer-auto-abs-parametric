require('lubridate')

#aart 2017 - 25 month loans

dir.create("./data-clean/") #to store results

rm(list=ls())

source("./code/default_time.R")

loan_term_c = 25 #left-truncation only; no right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

summary(aart$reportingPeriodBeginningLoanBalanceAmount)

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,] 


#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #2 loans

aart = aart[!(aart$assetNumber %in% bad_data),]

table(aart$Xc)
#remove loans with terms greater than 26 months (possible extensions)
aart$Xc = ifelse(aart$Xc >= loan_term_c + 1, loan_term_c + 1, aart$Xc)
table(aart$Xc)

obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

write.csv(obs_data, './data-clean/aart-2017-25mo.csv')

aart.2017.25mo.parameters = data.frame("delta" = delta,
                                       "m" = M,
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2017.25mo.parameters,
          './data-clean/aart-2017-25mo-trapezoid-dim.csv')

#aart 2017 - 50 month loans

rm(list=ls())

source("./code/default_time.R")

loan_term_c = 50 #left-truncation and right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,]

#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart = aart[aart$initialLoanAge <= 18,] #clean up likely loan extensions
aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #1 loan to remove

aart = aart[!(aart$assetNumber %in% bad_data),]

aart = aart[aart$risk_cat_ir == "super_prime",]
table(aart$Xc)

aart$Xc = ifelse(aart$Xc >= loan_term_c - 2, loan_term_c - 2, aart$Xc) #make 48 month loans
table(aart$Xc)


obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

write.csv(obs_data, './data-clean/aart-2017-50mo.csv')

aart.2017.50mo.parameters = data.frame("delta" = delta,
                                       "m" = M,
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2017.50mo.parameters,
          './data-clean/aart-2017-50mo-trapezoid-dim.csv')


#aart 2017 - 73 month loans
rm(list=ls())

source("./code/default_time.R")

loan_term_c = 73 #left-truncation and right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,]

#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

#aart = aart[aart$initialLoanAge <= 18,] #clean up likely loan extensions
aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #212 loan to remove

aart = aart[!(aart$assetNumber %in% bad_data),]

aart = aart[aart$risk_cat_ir == "super_prime",]
table(aart$Xc)

aart$Xc = ifelse(aart$Xc >= loan_term_c + 1, loan_term_c + 1, aart$Xc) #make 75 month loans
table(aart$Xc)


obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

write.csv(obs_data, './data-clean/aart-2017-73mo.csv')

aart.2017.73mo.parameters = data.frame("delta" = delta,
                                       "m" = min(M, max(obs_data$Z) - delta),
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2017.73mo.parameters,
          './data-clean/aart-2017-73mo-trapezoid-dim.csv')











################################################################################

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

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

lam_hat = vector()
for (i in c((delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))


obs_data$D = 1 - obs_data$C
len_obs_window = 43
e = len_obs_window + (M + delta)
tau = e - (M + delta + 1)
Delta = delta; m = M
omega = max(aart$Xc)

names(obs_data) = c("Xi", "Yi", "Ci", "Di")

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

a1 = c()
for(k in c((Delta+1):(Delta+m))){
  a1 = append(a1, (k - (Delta+1)) * h_dot_v(k))
}
a2 = c()
for(j in c((Delta + 1):omega)){
  a2 = append(a2, (j - (Delta + 1)) * h_u_dot(j))
}
b1 = c()
for(j in c((Delta + 1):(omega - 1))){
  b1 = append(b1, h_u_dot(j))
}
a = sum(a1) - sum(a2)
b = sum(b1)

p_hat = b / (b - a)

G_hat = c()
for(k in c((Delta+1):(Delta+m))){
  A = h_dot_v(k) * (1 - (b/a))^(k - (Delta + 1))
  B = mapply(h_dot_v, c((Delta+1):(Delta+m)))
  C = (1 - (b/a))^(c((Delta+1):(Delta+m)) - (Delta + 1))
  G_hat = append(G_hat, A * (sum( B * C ))^(-1) )
}

p + geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                     y = p_hat, yend = p_hat), color = "red",
                 linetype = "dashed")

#confidence intervals

#estimate parameters
THETA_est = c(p_hat, G_hat)

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

#Vn estimate
Vi = c()
for(i in c(1:n)){
  
  Xi = obs_data$Xi[i]
  Yi = obs_data$Yi[i]
  
  Vi = append(Vi, (psi(Xi, Yi, THETA_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)
#(0.02262015, 0.03991475)

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")

df$year_lt = "AART-2017-25M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1725.csv")

################################################################################
# left-truncation: 2019 DATA
################################################################################


rm(list=ls())

source("./code/default_time.R")

loan_term_c = c(25)
len_obs_window = 46 #num. mnths in obs. window from CRC paper

path = "./data/"
aart <- read.csv(paste(path,'aart193_compiledr.csv',sep=""))

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,] 

#calculate remaining payments
aart_trust_start_date = "08-01-2019"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #24 'bad' loans of 2,195

aart = aart[!(aart$assetNumber %in% bad_data),]

table(aart$Xc)
aart$Xc = ifelse(aart$Xc >= 26, 26, aart$Xc)
table(aart$Xc)

obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

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

lam_hat = vector()
for (i in c((delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))


obs_data$D = 1 - obs_data$C
len_obs_window = 46
e = len_obs_window + (M + delta)
tau = e - (M + delta + 1)
Delta = delta; m = M
omega = max(aart$Xc)

names(obs_data) = c("Xi", "Yi", "Ci", "Di")

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

a1 = c()
for(k in c((Delta+1):(Delta+m))){
  a1 = append(a1, (k - (Delta+1)) * h_dot_v(k))
}
a2 = c()
for(j in c((Delta + 1):omega)){
  a2 = append(a2, (j - (Delta + 1)) * h_u_dot(j))
}
b1 = c()
for(j in c((Delta + 1):(omega - 1))){
  b1 = append(b1, h_u_dot(j))
}
a = sum(a1) - sum(a2)
b = sum(b1)

p_hat = b / (b - a)

G_hat = c()
for(k in c((Delta+1):(Delta+m))){
  A = h_dot_v(k) * (1 - (b/a))^(k - (Delta + 1))
  B = mapply(h_dot_v, c((Delta+1):(Delta+m)))
  C = (1 - (b/a))^(c((Delta+1):(Delta+m)) - (Delta + 1))
  G_hat = append(G_hat, A * (sum( B * C ))^(-1) )
}

p + geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                     y = p_hat, yend = p_hat), color = "red",
                 linetype = "dashed")

#confidence intervals

#estimate parameters
THETA_est = c(p_hat, G_hat)

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

#Vn estimate
Vi = c()
for(i in c(1:n)){
  
  Xi = obs_data$Xi[i]
  Yi = obs_data$Yi[i]
  
  Vi = append(Vi, (psi(Xi, Yi, THETA_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)
#(0.02262015, 0.03991475) #2017

df$year_lt = "AART-2019-25M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1925.csv")


################################################################################
#RIGHT-CENSORING
################################################################################

# require('lubridate')
# require('ggplot2')
# require('extrafont') #may need to load fonts

rm(list=ls())

source("./code/default_time.R")

loan_term_c = 50 #left-truncation and right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,]

#calculate remaining payments
aart_trust_start_date = "06-01-2017"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart = aart[aart$initialLoanAge <= 18,] #clean up likely loan extensions
aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #24 'bad' loans of 2,195

aart = aart[!(aart$assetNumber %in% bad_data),]

aart = aart[aart$risk_cat_ir == "super_prime",]
table(aart$Xc)

aart$Xc = ifelse(aart$Xc >= loan_term_c - 2, loan_term_c - 2, aart$Xc) #make 48 month loans
table(aart$Xc)


obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

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

lam_hat = vector()
for (i in c((delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))

obs_data$D = 1 - obs_data$C
len_obs_window = 43
e = len_obs_window + (M + delta)
tau = e - (M + delta + 1)
Delta = delta; m = M
omega = max(aart$Xc)

names(obs_data) = c("Zi", "Yi", "Ci", "Di")

#censoring check for omega
obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0), 1, obs_data$Di)

#from corollary 4.2.2
gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

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


g_tau_MLE = function(v, a, b){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = gnv(v) * (1 - (b/a))^(v - (Delta + 1))
  B = mapply(gnv, c(v_min:v_max))
  C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
  
  return( A * (sum( B * C ))^(-1) )
  
}


G_hat = c()
for(k in c((Delta+1):(Delta+m))){
  G_hat = append(G_hat, g_tau_MLE(k, a, b) )
}

#confidence intervals
THETA_est = c(p_hat, G_hat)

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

Vi = c()
for(i in c(1:n)){
  
  Zi = obs_data$Zi[i]
  Yi = obs_data$Yi[i]
  Di = obs_data$Di[i]
  
  Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)

df$year_lt = "AART-2017-50M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1750.csv")

df = df[-1,]

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")


################################################################################
# left-truncation & right-censoring: 2019 DATA
################################################################################

rm(list=ls())

source("./code/default_time.R")

loan_term_c = 50 #left-truncation and right-censoring
len_obs_window = 46 #num. mnths in obs. window

path = "./data/"
aart <- read.csv(paste(path,'aart193_compiledr.csv',sep=""))

date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
min(date); max(date)


min(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
mean(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
median(as.numeric(aart$obligorCreditScore), na.rm = TRUE)
max(as.numeric(aart$obligorCreditScore), na.rm = TRUE)

min(aart$originalLoanTerm); max(aart$originalLoanTerm)

aart <- aart[aart$originalLoanTerm == loan_term_c,]

#calculate remaining payments
aart_trust_start_date = "08-01-2019"
date <- paste(aart$originationDate,"-01",sep="")
date <- as.Date(date, "%m/%Y-%d")
age = interval(date,as.Date(aart_trust_start_date,"%m-%d-%Y")) %/% months(1)
aart$initialLoanAge = age

aart = aart[aart$initialLoanAge <= 18,] #clean up likely loan extensions
aart$remainingTermtoMaturityNumber = aart$originalLoanTerm - aart$initialLoanAge

#create credit risk categories
aart$risk_cat_ir <- as.factor(
  ifelse(aart$originalInterestRatePercentage<0.05,"super_prime",
         ifelse(aart$originalInterestRatePercentage<0.10,"prime",
                ifelse(aart$originalInterestRatePercentage<0.15,"near_prime",
                       ifelse(aart$originalInterestRatePercentage<0.20,"subprime","deep_subprime")))))

delta = loan_term_c - max(aart$remainingTermtoMaturityNumber) 
M = loan_term_c - min(aart$remainingTermtoMaturityNumber) - delta
T_start = M + delta + aart$remainingTermtoMaturityNumber - loan_term_c 
Y = M + delta - T_start + 1

######################################################################
######################################################################
######################################################################
# algorithm to find loan outcomes (def, repay, cens)
######################################################################
######################################################################
######################################################################
X = vector()
C = vector()
D = vector()
R = vector()

for (j in c(1:nrow(aart))) {
  c_bond = default_time(aart[j,])
  X = append(X, c_bond[1])
  C = append(C, c_bond[2])
  R = append(R, c_bond[3])
  D = append(D, c_bond[4])
}

#shift back to the original timeline
Xc = M + delta + X - T_start + 1

######################################################################
######################################################################
aart = cbind(aart,Y,X,Xc,C,D,R)

a_cens = aart[aart$C == 1,]
n = nrow(a_cens)
check = c()

for (i in c(1:n)) {
  b_dat = a_cens[i,]
  
  final_bal = as.numeric(b_dat[1,paste("BAL",len_obs_window,sep="")])
  check = append(check,
                 ifelse(is.na(final_bal),"check",0))
}

bad_data = a_cens$assetNumber[check == "check"]
length(bad_data) #24 'bad' loans of 2,195

aart = aart[!(aart$assetNumber %in% bad_data),]

aart = aart[aart$risk_cat_ir == "super_prime",]
table(aart$Xc)
aart$Xc = ifelse(aart$Xc >= loan_term_c - 2, loan_term_c - 2, aart$Xc) #make 48 month loans
table(aart$Xc)


obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

#get unique observations of Z
#note: Z = min(X_i, C_i)
z = sort(unique(obs_data$Z))

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

lam_hat = vector()
for (i in c((delta+1):(max(z)))) {
  lam_hat = append(lam_hat,est_haz(i))
}

Var_est <- function(x){
  num = f_star(x) * (est_C(x) - f_star(x))
  den = (est_C(x))^3
  return(num/den)
}

Var_hat = vector()
for (i in c((delta+1):(max(z)))) {
  Var_hat = append(Var_hat, Var_est(i))
}

est_dist_true = data.frame(
  "Age" = c((delta+1):(max(z))),
  "lam_hat" = lam_hat,
  "Var_hat" = Var_hat
)

est_dist_true$lam_hat[nrow(est_dist_true)] = 1
est_dist_true$Var_hat[nrow(est_dist_true)] = 0

#remove zero estimates of lambda - OK?
est_dist_true = est_dist_true[est_dist_true$lam_hat > 0,]

age = est_dist_true$Age
lam_hat = est_dist_true$lam_hat
Var_hat = est_dist_true$Var_hat

CI_lower_log = log(lam_hat) - qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)
CI_upper_log = log(lam_hat) + qnorm(0.975) * sqrt( (Var_hat/(lam_hat)^2) / n)

#plotting via hazard rates
est_dist = data.frame(
  "Age" = age,
  "lam_hat" = lam_hat,
  "Est_Var" = Var_hat,
  "CI_lower" = exp(CI_lower_log),
  "CI_upper" = exp(CI_upper_log)
)

est_dist$lam_hat[nrow(est_dist)] = 1


df = data.frame("age" = est_dist$Age, "lam_hat" = est_dist$lam_hat,
                "ci_low" = est_dist$CI_lower, "ci_high" = est_dist$CI_upper)
#remove lam_hat = 1
df = df[c(1:(nrow(df)-1)),]

p <-
  ggplot() +
  geom_line(data=df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  #facet_wrap(vars(window)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))

obs_data$D = 1 - obs_data$C
len_obs_window = 46
e = len_obs_window + (M + delta)
tau = e - (M + delta + 1)
Delta = delta; m = M
omega = max(aart$Xc)

names(obs_data) = c("Zi", "Yi", "Ci", "Di")

#censoring check for omega
obs_data[(obs_data$Zi == omega) & (obs_data$Di == 0),]
obs_data$Di = ifelse((obs_data$Zi == omega) & (obs_data$Di == 0), 1, obs_data$Di)

#from RC theorem (add #)
gnv = function(v){
  return((1/nrow(obs_data)) * sum(obs_data$Yi == v))
}

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


g_tau_MLE = function(v, a, b){
  
  v_min = Delta + 1
  v_max = m + Delta
  
  A = gnv(v) * (1 - (b/a))^(v - (Delta + 1))
  B = mapply(gnv, c(v_min:v_max))
  C = (1 - (b/a))^(c(v_min:v_max) - (Delta + 1))
  
  return( A * (sum( B * C ))^(-1) )
  
}


G_hat = c()
for(k in c((Delta+1):(Delta+m))){
  G_hat = append(G_hat, g_tau_MLE(k, a, b) )
}

#confidence intervals
THETA_est = c(p_hat, G_hat)

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

Vi = c()
for(i in c(1:n)){
  
  Zi = obs_data$Zi[i]
  Yi = obs_data$Yi[i]
  Di = obs_data$Di[i]
  
  Vi = append(Vi, (psi(Yi, Zi, Di, THETA_est))^2 )
  
}

Vn = (1/n) * sum(Vi)
sigma = sqrt(1 / Vn)

#CI estimates
CI_low = THETA_est[1] - qnorm(0.975) * sigma / sqrt(n)
CI_upp = THETA_est[1] + qnorm(0.975) * sigma / sqrt(n)
c(CI_low, CI_upp)

df$year_lt = "AART-2019-50M"
df$phat = p_hat
df$p_upp = CI_upp
df$p_low = CI_low

write.csv(df, "./results/df1950.csv")

p + 
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = p_hat, yend = p_hat), color = "red",
               linetype = "dashed") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_low, yend = CI_low), color = "red",
               linetype = "solid") +
  geom_segment(aes(x = df$age[1], xend = df$age[length(df$age)],
                   y = CI_upp, yend = CI_upp), color = "red",
               linetype = "solid")


################################################################################
# construct the figure
################################################################################
# require('cowplot')
# require('ggplot2')
# require('extrafont')

rm(list=ls())

df1 = read.csv('./results/df1725.csv')
df1 = df1[,-1]
df2 = read.csv('./results/df1925.csv')
df2 = df2[,-1]
df3 = read.csv('./results/df1750.csv')
df3 = df3[,-1]
df4 = read.csv('./results/df1950.csv')
df4 = df4[,-1]

plot_df = rbind(df1, df2)
plot_df = plot_df[plot_df$age >= 8,]

p1 <-
  ggplot() +
  geom_line(data=plot_df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  geom_line(data=plot_df, aes(x=age, y=phat), color="red", linetype = "dashed") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=p_low, ymax=p_upp),
              fill="red", alpha=0.35) +
  facet_grid(cols=vars(year_lt)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))# +

plot_df = rbind(df3, df4)
plot_df = plot_df[-1,] #remove NA row
plot_df = plot_df[plot_df$age >= 7,]

p2 <-
  ggplot() +
  geom_line(data=plot_df, aes(x=age, y=lam_hat), color="blue") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=ci_low, ymax=ci_high),
              fill="lightblue", alpha=0.5) +
  geom_line(data=plot_df, aes(x=age, y=phat), color="red", linetype = "dashed") +
  geom_ribbon(data=plot_df, aes(x=age, ymin=p_low, ymax=p_upp),
              fill="red", alpha=0.35) +
  facet_grid(cols=vars(year_lt)) +
  xlab("Loan Age") + ylab("Estimated Hazard Rate") +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9,family="Times New Roman"),
        strip.text=element_text(size=9,family="Times New Roman"),
        axis.text=element_text(size=9,family="Times New Roman"))# +

set_null_device(cairo_pdf)
plot_grid(p1, p2, nrow=2)
ggsave("./results/aart_comp.pdf",height=4,width=6,device=cairo_pdf)
file.remove('./Rplot001.pdf')











