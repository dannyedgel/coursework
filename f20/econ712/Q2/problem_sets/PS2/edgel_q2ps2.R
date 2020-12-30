###
### This file is used to generate charts for Econ 712 problem set #2 of Q2
###
### Econ 712: Macroeconomics I
### Fall 2020
###
### Date created:  17 Nov 2020
### Last modified: 17 Nov 2020
### Author: Danny Edgel
###

### clear workspace
rm(list=ls())

### declare package libraries
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/Q2PS2/')


###
### Question 3
###____________________________________________

### vector for number of scenarios (automating a, b, and c)
runs <- c(1) #,2,3)

### choose number of elements used to graph
n <- 10000

### define constants
beta  <- 0.95
delta <- 0.1
z     <- 1
gamma <- 2

### set labor supply equal to one
N <- 1

### guess maximum capital level; define grid
maxkap <- 100
minkap <- 0
kgrid  <- seq(minkap,maxkap, length.out = n) 


### loop through scenarios
for ( run in runs ){
  
  ## make current run's adjustment (defaults to 3(a))
  if ( run == 2 ){
    gamma <- 1.01
  }else if ( run == 3 ){
    gamma <- 2
    z     <- 1.2
  }
  
  ## initialize variables
  v  <- tv <- decis <- tdecis <- rep(0, times = n)
  cons <- util <- vint <- array(rep(0,times = n*n), c(n,n))
  
  ## iteratre on Bellman equation
  while test == 0{
    cons <- w*N + array(rep((1 + r - delta)*kgrid,times = n),c(n,n)) - 
      array(rep(kgrid,times = n),c(n,n))
  }
  
} ## end current run


# 
# dat <- data.frame(
#   k   = c(k,k,k),
#   v   = c(I,ss1,ss2),
#   grp = c(rep('Resource Constraint',times = n),
#           rep('Saddle @ R=1.05',times = n),
#           rep('Saddle @ R=1.10',times = n))
# )
# 
# png(file='figure3c.png', width = 450, height = 250)
# ggplot(dat) + 
#   geom_line(aes(k,v, group = grp, color = grp)) + 
#   geom_vline(xintercept = k[abs(ss1-I) == min(abs(ss1-I))],
#              linetype = 'dashed') +
#   geom_segment(data = data.frame(x1 = k[abs(ss1-I) == min(abs(ss1-I))],
#                                  y1 = ss2[abs(ss1-I) == min(abs(ss1-I))],
#                                  x2 = k[abs(ss2-I) == min(abs(ss2-I))],
#                                  y2 = ss2[abs(ss2-I) == min(abs(ss2-I))]),
#                aes(x = x1, y = y1, xend = x2, yend = y2),
#                color = 'red', arrow = arrow(length = unit(0.2,"cm")), size = 1)+
#   theme_classic() + ylab('I') + theme(legend.title = element_blank()) 
# dev.off()
# 
