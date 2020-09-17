
rm(list=ls())

filepath <- 'C:/Users/edgel/Google Drive/UW-Madison/Fall 2020/Econ 712 (Macro)/Problem Sets/PS1/'

library(ggplot2)

r <- 0.01 # interest rate
d <- 1 # constant dividend
dim <- 299 # terminal period t <- 99
pstar <- 100


tvector <- 0:dim
at      <- rep(1+r,dim+1)
for ( x in 0:dim ){at[x+1]<-at[x+1]^tvector[x+1]}
p0_0    <- pstar
p0_1    <- pstar - 1
p0_2    <- pstar + 1

p_t_0 <- rep(pstar,dim+1)
p_t_1 <- (p0_1 - pstar)*at + pstar
p_t_2 <- (p0_2 - pstar)*at + pstar

pt0 <- p_t_0[1:dim]
pt0_p1 <- p_t_0[2:(dim+1)]
pt1 <- p_t_1[1:dim]
pt1_p1 <- p_t_1[2:(dim+1)]
pt2 <- p_t_2[1:dim]
pt2_p1 <- p_t_2[2:(dim+1)]

dat  <- data.frame(pt0,pt0_p1,pt1,pt1_p1,pt2,pt2_p1)
tdat <- data.frame(tvector,p_t_0,p_t_1,p_t_2)

## Problem 2:
png(file=paste0(filepath,'problem2_phaseplot.png'), width = 300, height = 244)
ggplot(dat) +
  geom_point(aes(pt0,pt0_p1)) +
  geom_line(aes(pt1,pt1_p1), color='red') +
  geom_line(aes(pt2,pt2_p1), color='blue') +
  xlab(expression(p[t])) + ylab(expression(p[t+1])) +
  theme_classic()
dev.off()

png(file=paste0(filepath,'problem2_timeplot.png'), width = 350, height = 170)
ggplot(tdat) +
  geom_line(aes(tvector,p_t_0)) +
  geom_line(aes(tvector,p_t_1), color='red') +
  geom_line(aes(tvector,p_t_2), color='blue') +
  xlab('t') + ylab(expression(p[t])) +
  theme_classic()
dev.off()

### Problem 3:
dim <- 99
t <- 0:dim

p0 <- 100 # initial price
pvector <- rep(0,dim+1) # creating a vector of price from t<-0 to t<-99
pvector[1] <- p0 # giving value to the first element of the price vector

p0_1 <- 90
p0_2 <- 110

p1 <- p2 <- pvector
p1[1] <- p0_1
p2[1] <- p0_2


for (n in 2:(dim+1) ){ # starting from t <- 1 to t <- 100 
    pvector[n] <- (1+r)*pvector[n-1]-d # updating the price in the next period 
                                       # with the first-order difference equation
    p1[n] <- (1+r)*p1[n-1]-d
    p2[n] <- (1+r)*p2[n-1]-d
}


tdat <- data.frame(t,pvector,p1,p2)
tdat$pvector <- pvector

png(file=paste0(filepath,'problem3_timeplot.png'), width = 350, height = 170)
ggplot(tdat) +
  geom_line(aes(t,pvector)) +
  geom_line(aes(t,p1),color='blue') +
  geom_line(aes(t,p2),color='red') +
  ggtitle('Price Dynamics') +
  xlab('Time t') + ylab(expression(Price~p[t]))+
  scale_x_continuous(breaks = seq(0,100,10)) +
  theme_classic()
dev.off()


### Problem 4:
r1 <- 0.01 # initial discount rate
r2 <- 0.02 # second discount rate

pstar1 <- d / r1 # initial pstar
pstar2 <- d / r2 # second pstar

pvector <- rep(0,dim+1) # creating a vector of price from t<-0 to t<-99
pvector[dim+1] <- pstar2 # giving value to the last element of the price vector

for (n in (dim+1):20 ){ 
  ### Beginning in period 49, investors price 
  #if ( n < 50 ){ pvector[n] <- (1+r1)*pvector[n-1]-d
  #}else{         pvector[n] <- (1+r2)*pvector[n-1]-d }
  pvector[n-1]<-(d+pvector[n])/(1+r2)
  if ( n >= 50 ){ pvector[n-1]<-(d+pvector[n])/(1+r2)
  }else{ pvector[n-1]<-(d+pvector[n])/(1+r1) }
}

### replace pre-t=20 entries with r=0.01 steady state value
pvector[1:19] <- pstar1


### plot

png(file=paste0(filepath,'problem4_timeplot.png'), width = 350, height = 170)
tdat <- data.frame(t,pvector)
ggplot(tdat,aes(t,pvector)) +
  geom_line(color='blue') +
  ggtitle('Price Dynamics') +
  xlab('Time t') + ylab(expression(Price~p[t])) +
  scale_x_continuous(breaks = seq(0,100,10)) +
  scale_y_continuous(breaks = seq(20,120,20)) + ylim(40,110) +
  theme_classic()
dev.off()



