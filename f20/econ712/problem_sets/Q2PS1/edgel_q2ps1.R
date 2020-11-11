###
### This file is used to generate charts for Econ 712 problem set #1 of Q2
###
### Econ 712: Macroeconomics I
### Fall 2020
###
### Date created:  10 Nov 2020
### Last modified: 10 Nov 2020
### Author: Danny Edgel
###

### clear workspace
rm(list=ls())

### declare package libraries
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/Q2PS1/')


###
### Question 3c
###____________________________________________

### choose number of elements used to graph
n <- 1000

### define constants
R1    <- 1.05   
R2    <- 1.1
delta <- .2
kstar <- 100
Pi0   <- 1000

# ### define marginal profit and investment
# I_ss <- seq(10,20,length.out = n) # vector of possible steady-state I values
# k_ss <- (R+delta-1)*I_ss + kstar     # steady-state capital at each steady-state I

### define a vector of possible investment amounts at each level of capital
k <- seq(100,110,length.out = n+1)
V <- rep(0,times = n)
I <- rep(0,times = n)
ss1 <- ss2 <- rep(0,times = n)
for ( i in 1:n ){
  V[i] <- Pi0-.5*(k[i]-kstar)^2 - .5*(k[i+1]-(1-delta)*k[i])^2
  I[i] <- k[i+1]-(1-delta)*k[i]
  
  ## identify steady-state values
  ss1[i] <- (k[i]-kstar)/(R1+delta-1)
  ss2[i] <- (k[i]-kstar)/(R2+delta-1)
}
k <- k[1:n]

dat <- data.frame(
  k   = c(k,k,k),
  v   = c(I,ss1,ss2),
  grp = c(rep('Resource Constraint',times = n),
          rep('Saddle @ R=1.05',times = n),
          rep('Saddle @ R=1.10',times = n))
)

png(file='figure3c.png', width = 450, height = 250)
ggplot(dat) + 
  geom_line(aes(k,v, group = grp, color = grp)) + 
  geom_vline(xintercept = k[abs(ss1-I) == min(abs(ss1-I))],
             linetype = 'dashed') +
  geom_segment(data = data.frame(x1 = k[abs(ss1-I) == min(abs(ss1-I))],
                                 y1 = ss2[abs(ss1-I) == min(abs(ss1-I))],
                                 x2 = k[abs(ss2-I) == min(abs(ss2-I))],
                                 y2 = ss2[abs(ss2-I) == min(abs(ss2-I))]),
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = 'red', arrow = arrow(length = unit(0.2,"cm")), size = 1)+
  theme_classic() + ylab('I') + theme(legend.title = element_blank()) 
dev.off()

