library(reshape2) #data melting
library(cowplot) #preparing the plot

################################################################################
#formulas
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

################################################################################
#search round 1
p.start = 0.50
p.end = 0.999
seq.length = 100
p.val = seq(p.start, p.end, by = (p.end - p.start)/seq.length)
b.start = 1.00
b.end = 2.00
b.val = seq(b.start, b.end, by = (b.end - b.start)/seq.length)

z <- outer(p.val,b.val,P_constraint_3);
rownames(z) = p.val
colnames(z) = b.val

# Generate data
dat.plot <- melt(z)
names(dat.plot) <- c("p1", "p2", "z")

# Basic plot
p1 <-
  ggplot(dat.plot, aes(p1, p2, z = z)) +  
  stat_contour(geom="polygon", aes(fill=..level..)) +
  xlab( expression( italic(p)[1] ) ) +
  ylab( expression( italic(p)[2] ) ) +
  ylim( c(min(b.val), max(b.val)) ) +
  xlim( c(min(p.val), max(p.val)) ) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.title=element_blank(), 
        legend.position = "bottom",
        #legend.key.width = unit(dev.size()[1] / 40, "inches"),
        legend.key.height = unit(0.2, "cm"))

#get minimum value
which(z == min(z, na.rm = TRUE), arr.ind=TRUE)
p.val[which(z == min(z, na.rm = TRUE), arr.ind=TRUE)[1]]
b.val[which(z == min(z, na.rm = TRUE), arr.ind=TRUE)[2]]
P_constraint_2(p.val[which(z == min(z, na.rm = TRUE), arr.ind=TRUE)[1]],
               b.val[which(z == min(z, na.rm = TRUE), arr.ind=TRUE)[2]])

#search round 2
p.start = 0.95
p.end = 0.999
seq.length = 100
p.val = seq(p.start, p.end, by = (p.end - p.start)/seq.length)
b.start = 1.00
b.end = 1.25
b.val = seq(b.start, b.end, by = (b.end - b.start)/seq.length)

z <- outer(p.val,b.val,P_constraint_3);
rownames(z) = p.val
colnames(z) = b.val

# Generate data
dat.plot <- melt(z)
names(dat.plot) <- c("p1", "p2", "z")

# Basic plot
p2 <-
  ggplot(dat.plot, aes(p1, p2, z = z)) +  
  stat_contour(geom="polygon", aes(fill=..level..)) +
  xlab( expression( italic(p)[1] ) ) +
  ylab( expression( italic(p)[2] ) ) +
  ylim( c(min(b.val), max(b.val)) ) +
  xlim( c(min(p.val), max(p.val)) ) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.title=element_blank(), 
        legend.position = "bottom",
        #legend.key.width = unit(dev.size()[1] / 40, "inches"),
        legend.key.height = unit(0.2, "cm"))

#get minimum value
which(z == min(z), arr.ind=TRUE)
p.val[which(z == min(z), arr.ind=TRUE)[1]]
b.val[which(z == min(z), arr.ind=TRUE)[2]]
P_constraint_2(p.val[which(z == min(z), arr.ind=TRUE)[1]],
               b.val[which(z == min(z), arr.ind=TRUE)[2]])

p.con.min = 0.02
which(z <= p.con.min, arr.ind=TRUE)

p.val[which(z <= p.con.min, arr.ind=TRUE)[,1]]
b.val[which(z <= p.con.min, arr.ind=TRUE)[,2]]

#search round 3
p.start = 0.975
p.end = 0.990
seq.length = 100
p.val = seq(p.start, p.end, by = (p.end - p.start)/seq.length)
b.start = 1.04
b.end = 1.25
b.val = seq(b.start, b.end, by = (b.end - b.start)/seq.length)

z <- outer(p.val,b.val,P_constraint_3);
rownames(z) = p.val
colnames(z) = b.val

# Generate data
dat.plot <- melt(z)
names(dat.plot) <- c("p1", "p2", "z")

# Basic plot
p3 <-
  ggplot(dat.plot, aes(p1, p2, z = z)) +  
  stat_contour(geom="polygon", aes(fill=..level..)) +
  xlab( expression( italic(p)[1] ) ) +
  ylab( expression( italic(p)[2] ) ) +
  ylim( c(min(b.val), max(b.val)) ) +
  xlim( c(min(p.val), max(p.val)) ) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.title=element_blank(), 
        legend.position = "bottom",
        #legend.key.width = unit(dev.size()[1] / 40, "inches"),
        legend.key.height = unit(0.2, "cm"))

#get minimum value
which(z == min(z), arr.ind=TRUE)
p.val[which(z == min(z), arr.ind=TRUE)[1]]
b.val[which(z == min(z), arr.ind=TRUE)[2]]
P_constraint_2(p.val[which(z == min(z), arr.ind=TRUE)[1]],
               b.val[which(z == min(z), arr.ind=TRUE)[2]])

p.con.min = 0.015
which(z <= p.con.min, arr.ind=TRUE)

p.val[which(z <= p.con.min, arr.ind=TRUE)[,1]]
b.val[which(z <= p.con.min, arr.ind=TRUE)[,2]]

#search round 4
p.start = 0.984
p.end = 0.990
seq.length = 100
p.val = seq(p.start, p.end, by = (p.end - p.start)/seq.length)
b.start = 1.12
b.end = 1.24
b.val = seq(b.start, b.end, by = (b.end - b.start)/seq.length)

z <- outer(p.val,b.val,P_constraint_3);
rownames(z) = p.val
colnames(z) = b.val

# Generate data
dat.plot <- melt(z)
names(dat.plot) <- c("p1", "p2", "z")

# Basic plot
p4 <-
  ggplot(dat.plot, aes(p1, p2, z = z)) +  
  stat_contour(geom="polygon", aes(fill=..level..)) +
  xlab( expression( italic(p)[1] ) ) +
  ylab( expression( italic(p)[2] ) ) +
  ylim( c(min(b.val), max(b.val)) ) +
  xlim( c(min(p.val), max(p.val)) ) +
  theme_bw() +
  theme(axis.title.x=element_text(size=9, family="Times New Roman"),
        axis.title.y=element_text(size=9, family="Times New Roman"),
        axis.text.x=element_text(size=9, family="Times New Roman"),
        axis.text.y=element_text(size=9, family="Times New Roman"),
        legend.text=element_text(size=9, family="Times New Roman"),
        legend.title=element_blank(), 
        legend.position = "bottom",
        #legend.key.width = unit(dev.size()[1] / 40, "inches"),
        legend.key.height = unit(0.2, "cm"))

#get minimum value
which(z == min(z), arr.ind=TRUE)
p.val[which(z == min(z), arr.ind=TRUE)[1]]
b.val[which(z == min(z), arr.ind=TRUE)[2]]
P_constraint_2(p.val[which(z == min(z), arr.ind=TRUE)[1]],
               b.val[which(z == min(z), arr.ind=TRUE)[2]])

plot_grid(p1, p2, p3, p4, nrow = 2)

ggsave("dw-optimization-supplement.pdf",height=4,width=6,device = cairo_pdf)

