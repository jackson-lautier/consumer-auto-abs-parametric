######################################################################################
######################################################################################
######################################################################################
# Data processing scripts to produce the data used in the manuscript:
# 
# "ESTIMATING THE TIME-TO-EVENT DISTRIBUTION FOR LOAN-LEVEL
# DATA WITHIN A CONSUMER AUTO LOAN ASSET-BACKED SECURITY"
#
# LAUTIER, POZDNYAKOV, YAN
# 2025
#
#R version 4.3.2 (2023-10-31 ucrt) -- "Eye Holes"
#RStudio 2023.03.0+386 "Cherry Blossom" Release
#(3c53477afb13ab959aeb5b34df1f10c237b256c3, 2023-03-09) for Windows
######################################################################################
######################################################################################
######################################################################################
######################################################################################
# INSTRUCTIONS
#
# supporting files:
# ".\raw-data\aart173_compiledr.csv"
# ".\raw-data\aart193_compiledr.csv"
#
# "./code/default_time.R"
#
#The code must be run sequentially downwards.
#As the new, cleaned files are prepared, they will be saved in a new
#folder 'processed-data' in the wd.
#For data analysis, proceed directly to 'data_analysis.R'.
#
#
######################################################################################
######################################################################################
######################################################################################
######################################################################################
require('lubridate')


######################################################################################
######################################################################################
######################################################################################
######################################################################################

#where processed data will be stored
dir.create('./processed-data/')

######################################################################################
#aart 2017 - 25 month loans
######################################################################################
rm(list=ls())

source("./code/default_time.R")

loan_term_c = 25 #left-truncation only; no right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

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

write.csv(obs_data, './processed-data/aart-2017-25mo.csv')

aart.2017.25mo.parameters = data.frame("delta" = delta,
                                       "m" = M,
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2017.25mo.parameters,
          './processed-data/aart-2017-25mo-trapezoid-dim.csv')

######################################################################################
#aart 2017 - 50 month loans
######################################################################################

rm(list=ls())

source("./code/default_time.R")

loan_term_c = 50 #left-truncation and right-censoring
len_obs_window = 43 #num. mnths in obs. window

path = "./raw-data/"
aart <- read.csv(paste(path,'aart173_compiledr.csv',sep=""))

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

#aart$Xc = ifelse(aart$Xc >= loan_term_c - 2, loan_term_c - 2, aart$Xc) #make 48 month loans
aart$Xc = ifelse(aart$Xc >= loan_term_c, loan_term_c, aart$Xc) #make 50 month loans
table(aart$Xc)


obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

write.csv(obs_data, './processed-data/aart-2017-50mo.csv')

aart.2017.50mo.parameters = data.frame("delta" = delta,
                                       "m" = M,
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2017.50mo.parameters,
          './processed-data/aart-2017-50mo-trapezoid-dim.csv')

######################################################################################
#aart 2019 - 25 month loans
######################################################################################

#2019 25mo loans
rm(list=ls())

source("./code/default_time.R")

loan_term_c = c(25)
len_obs_window = 46 #num. mnths in obs. window from CRC paper

path = "./raw-data/"
aart <- read.csv(paste(path,'aart193_compiledr.csv',sep=""))

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
length(bad_data) #2

aart = aart[!(aart$assetNumber %in% bad_data),]

table(aart$Xc)
aart$Xc = ifelse(aart$Xc >= 26, 26, aart$Xc)
table(aart$Xc)

obs_data <- data.frame(aart$Xc,aart$Y,aart$C)
names(obs_data)[names(obs_data) == 'aart.Xc'] <- 'Z'
names(obs_data)[names(obs_data) == 'aart.Y'] <- 'Y'
names(obs_data)[names(obs_data) == 'aart.C'] <- 'C'
n = nrow(obs_data)

write.csv(obs_data, './processed-data/aart-2019-25mo.csv')

aart.2019.25mo.parameters = data.frame("delta" = delta,
                                       "m" = M,
                                       "omega" = max(obs_data$Z),
                                       "e" = len_obs_window + (M + delta))

write.csv(aart.2019.25mo.parameters,
          './processed-data/aart-2019-25mo-trapezoid-dim.csv')