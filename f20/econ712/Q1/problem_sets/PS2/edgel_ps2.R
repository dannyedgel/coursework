###
### This file is used to complete problems for macro problem set #2
###
### Econ 712: Macroeconomics I
### Fall 2020
###
### Date created:  14 Sep 2020
### Last modified: 16 Sep 2020
### Author: Danny Edgel
###

### clear workspace
rm(list=ls())

### declare package libraries
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/Fall 2020/Econ 712 (Macro)/Problem Sets/PS2/')


###                                                           ###
### Problem 1                                                 ###
###___________________________________________________________###

### Define parameters
z     <- 1
alpha <- 0.3
delta <- 0.1
beta  <- 0.97

### Define solved steady-state capital and consumption variables
k_ss <- ((1/beta + delta - 1)/(alpha*z))^(1/(alpha-1))
c_ss <- z*(k_ss^alpha)-delta*k_ss

### report steady state variables:
cat(paste0('Steady state k: ',k_ss,'\nSteady state c: ',c_ss))


### define the steady-state Jacobian matrix for this system
J <- matrix(
  c(alpha*z*(k_ss^(alpha-1))+1-delta,
    beta*(alpha-1)*alpha*z*c_ss*((z*(k_ss^alpha)+(1-delta)*k_ss-c_ss)^(alpha-2))
      *(alpha*z*(k_ss^(alpha-1))+1-delta),
    -1,
    beta*(1-delta+alpha*z*((z*(k_ss^alpha)+(1-delta)*k_ss-c_ss)^(alpha-1)))-
      beta*(alpha-1)*alpha*z*c_ss*((z*(k_ss^alpha)+(1-delta)*k_ss-c_ss)^(alpha-2))
    ), nrow = 2
)

### use J's eigenvalues and eigenvectors to diagonalize J
diag <- matrix(c(eigen(J)$values[1],0,0,eigen(J)$values[2]),nrow=2)
P    <- as.matrix(data.frame(eigen(J)$vectors[,1],eigen(J)$vectors[,2]))
Pinv <- solve(P)

### Report PJP^-2
P
diag
Pinv
P %*% diag %*% Pinv
J

### Part 5: Determine evolution of k and c after a positive productivity shock
dim  <- 20    # number of time periods observed
tvec <- 1:dim # time vector
k    <- k_ss  # capital vector (entry 1)
c    <- c_ss  # consumption vector (entry 1)

z0   <- z + 0.1 # new productivity factor > z

k0_ss <- ((1/beta + delta - 1)/(alpha*z0))^(1/(alpha-1)) # new steady-state capital
c0_ss <- z0*(k0_ss^alpha)-delta*k0_ss # new steady-state consumption


# Repeat the Jacobian for the new steady state to extract new saddle path
J0 <- matrix(
  c(alpha*z0*(k0_ss^(alpha-1))+1-delta,
    beta*(alpha-1)*alpha*z0*c0_ss*((z0*(k0_ss^alpha)+(1-delta)*k0_ss-c0_ss)^(alpha-2))
    *(alpha*z0*(k0_ss^(alpha-1))+1-delta),
    -1,
    beta*(1-delta+alpha*z0*((z0*(k0_ss^alpha)+(1-delta)*k0_ss-c0_ss)^(alpha-1)))-
      beta*(alpha-1)*alpha*z0*c0_ss*((z0*(k0_ss^alpha)+(1-delta)*k0_ss-c0_ss)^(alpha-2))
  ), nrow = 2
)


diag0 <- matrix(c(eigen(J0)$values[1],0,0,eigen(J0)$values[2]),nrow=2)
P0    <- as.matrix(data.frame(eigen(J0)$vectors[,1],eigen(J0)$vectors[,2]))
Pinv0 <- solve(P0)

P0
diag0
Pinv0
P0 %*% diag0 %*% Pinv0
J0

# in each period, calculate capital and consumption using the last period's
# capital and consumption and the new productivity factor
for ( t in 2:dim ){
  if ( t < 5 ){
    k <- c(k, k_ss)
    c <- c(c, c_ss)
  }else{
    k <- c(
      k,
      P0[1,2]*(diag0[2,2]^t)*((k_ss-k0_ss)/(P0[1,2]*(diag0[2,2]^5)))+k0_ss
    )
    c <- c(
      c,
      P0[2,2]*(diag0[2,2]^t)*((k_ss-k0_ss)/(P0[1,2]*(diag0[2,2]^5)))+c0_ss
    )
  }
}



# plot capital against consumption
#dat <- data.frame(t,c,k,c0,k0,c0_ss,k0_ss)
dat <- data.frame(tvec,c,k)
x <- ggplot(dat)

png(file='problem1_timeplot_k.png', width = 525, height = 350)
x +
  geom_line(aes(tvec,k),color='blue') +
  geom_line(aes(tvec,rep(k_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(k0_ss,times=dim)),linetype='dashed') +
  ggtitle('Capital after a positive productivity shock') +
  ylab(expression(Capital~(k[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(k_ss-.1,k0_ss+.1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()

png(file='problem1_timeplot_c.png', width = 525, height = 350)
x +
  geom_line(aes(tvec,c),color='red') +
  geom_line(aes(tvec,rep(c_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(c0_ss,times=dim)),linetype='dashed') +
  ggtitle('Consumption after a positive productivity shock') +
  ylab(expression(Consumption~(c[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(c_ss-.03,c0_ss+.03), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()


### Part 6: computing consumption and capital using the nonlinear definitions

k    <- k_ss  # capital vector (entry 1)
c    <- c_ss  # consumption vector (entry 1)

for ( t in 2:dim ){
  if ( t < 5 ){
    k <- c(k, z*(k[t-1]^alpha)+(1-delta)*k[t-1]-c[t-1])
    c <- c(c, beta*c[t-1]*(1-delta+alpha*z*(z*(k[t-1]^alpha) + (1-delta)*k[t-1]-c[t-1])^(alpha-1)))
  }else if ( t == 5 ){
    k <- c(
      k,
      z * (k[t-1])^alpha + (1-delta) * k[t-1] - c[t-1]
    )
    c <- c(
      c,
      P0[2,2]*(diag0[2,2]^5)*((k_ss-k0_ss)/(P0[1,2]*(diag0[2,2]^5)))+c0_ss
    )
  }else{
    k <- c(
      k,
      z0 * (k[t-1])^alpha + (1-delta) * k[t-1] - c[t-1]
    )
    c <- c(
      c,
      beta * c[t-1] * ( 1 - delta + alpha*z0 * (z0 * (k[t-1]^alpha) + (1-delta) * k[t-1] - c[t-1] )^(alpha-1) )
    )
  }
}


dat <- data.frame(tvec,c,k)
x <- ggplot(dat)

png(file='problem1.6_timeplot_k.png', width = 525, height = 350)
x +
  geom_line(aes(tvec,k),color='blue') +
  geom_line(aes(tvec,rep(k_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(k0_ss,times=dim)),linetype='dashed') +
  ggtitle('Capital after a positive productivity shock (direct estimation)') +
  ylab(expression(Capital~(k[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(k_ss-.1,k0_ss+.1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()

png(file='problem1.6_timeplot_c.png', width = 525, height = 350)  
x +
  geom_line(aes(tvec,c),color='red') +
  geom_line(aes(tvec,rep(c_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(c0_ss,times=dim)),linetype='dashed') +
  ggtitle('Consumption after a positive productivity shock (direct estimation)') +
  ylab(expression(Consumption~(c[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(c_ss-.03,c0_ss+.03), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()
  

### Part 7: computing consumption and capital using the shooting method
N <- 20000
c.cand <- seq(from=c_ss,to=c0_ss,length.out = N)              # vector of candidate values
D      <- data.frame(cand = c.cand, dist = rep(0, times = N)) # norm distance vector w/
                                                              # associated c_t0 candidate

M <- 30 # length of the trajectory vector

for ( i in 1:N ){
  c_t0      <- D$cand[i]
  
  # fill in trajectory vector
  trj <- matrix(rep(0,times=2*M), nrow = 2)
  trj[1,1] <- k_ss
  trj[2,1] <- c_t0
  
  j <- 1
  while ( j < M ){
    j <- j + 1
    k <- trj[1,j-1]
    c <- trj[2,j-1]
    
    if ( sum(is.na(trj[,j])) == 0 & j <= M ){
      
      trj[1,j] <- z0*(k^alpha)+(1-delta)*k-c
      trj[2,j] <- beta*c*(1-delta+alpha*z0*(z0*(k^alpha) + (1-delta)*k-c)^(alpha-1))
      
    }else{ j<-M+1 }
    
  } # end trajectory loop
  
  D$dist[i] <- norm(trj[,j-1] - matrix(c(k0_ss,c0_ss),nrow=2))
} # end candidate loop

# Extract optimal c_t0 from distance matrix
cstar <- D$cand[D$dist == min(D$dist,na.rm = T)]
cstar <- cstar[is.na(cstar) == F] 

### implement shooting method to plot evolution of capital and consumption
### directly

k    <- k_ss  # capital vector (entry 1)
c    <- c_ss  # consumption vector (entry 1)


for ( t in 2:dim ){
  if ( t < 5 ){
    k <- c(k, z*(k[t-1]^alpha)+(1-delta)*k[t-1]-c[t-1])
    c <- c(c, beta*c[t-1]*(1-delta+alpha*z*(z*(k[t-1]^alpha) + (1-delta)*k[t-1]-c[t-1])^(alpha-1)))
  }else if ( t == 5 ){
    k <- c(
      k,
      z * (k[t-1])^alpha + (1-delta) * k[t-1] - c[t-1]
    )
    c <- c(
      c,
      cstar
    )
  }else{
    k <- c(
      k,
      z0 * (k[t-1])^alpha + (1-delta) * k[t-1] - c[t-1]
    )
    c <- c(
      c,
      beta * c[t-1] * ( 1 - delta + alpha*z0 * (z0 * (k[t-1]^alpha) + (1-delta) * k[t-1] - c[t-1] )^(alpha-1) )
    )
  }
}


dat <- data.frame(tvec,c,k)
x <- ggplot(dat)

png(file='problem1.6b_timeplot_k.png', width = 525, height = 350)
x +
  geom_line(aes(tvec,k),color='blue') +
  geom_line(aes(tvec,rep(k_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(k0_ss,times=dim)),linetype='dashed') +
  ggtitle('Capital after a positive productivity shock (shooting method)') +
  ylab(expression(Capital~(k[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(k_ss-.1,k0_ss+.1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()

png(file='problem1.6b_timeplot_c.png', width = 525, height = 350)  
x +
  geom_line(aes(tvec,c),color='red') +
  geom_line(aes(tvec,rep(c_ss,times=dim)),linetype='dashed') +
  geom_line(aes(tvec,rep(c0_ss,times=dim)),linetype='dashed') +
  ggtitle('Consumption after a positive productivity shock (shooting method)') +
  ylab(expression(Consumption~(c[t]))) +
  xlab('Time (t)') +
  scale_y_continuous(limits = c(c_ss-.03,c0_ss+.03), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1-.1,20+.1), expand = c(0, 0)) +
  theme_classic()
dev.off()
