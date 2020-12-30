
library( ggplot2 )
library( ggpubr )
library( rootSolve )

setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/PS6/')

rm(list = ls())

### LAST PROBLEM
n <- 1000
alpha <- seq(0,.5,length.out = n)
beta <- seq(0,1,length.out = n)

alpha <- alpha[2:(n-1)]
beta  <- beta[2:(n-1)]

n <- length(alpha)

### social planner's allocation
t_sp <- 1/2
ell_sp <- (2/3)*(1-alpha)
l_sp <- 1-ell_sp
c_sp <- (1-t_sp)*ell_sp
g_sp <- t_sp*ell_sp

u_sp <- log(l_sp) + log(alpha + c_sp) + log(alpha + g_sp)

### calculate tau for Ramsey Equilibrium at each alpha
foc <- function(t, a){ 
  
  u1 <- -2/(1-t+a)
  u2 <- (1-3*a-2*t)/(2*a+(1-3*a)*t-t^2)
  u3 <- 2/(1-t)
  
  return(u1 + u2 + u3) 
}

t_r <- c()
for (al in alpha) {
  t_r <- c(
    t_r,
    min(uniroot.all(foc, interval = c(beta[2],beta[n-1]), a = al))
  )
}



### RE allocation
t <- t_r

ell_r <- (1-t-alpha)/(2*(1-t))
l_r <- 1-ell_r
c_r <- (1/2)*(1-t-alpha)
g_r <- t*ell_r

u_r <- log(l_r) + log(alpha + c_r) + log(alpha + g_r)

### NE allocation
t_n <- t <- (1/2)
l_n <- (1/2)*(1+2*alpha)
ell_n <- 1- l_n
c_n <- (1-t_n)*ell_n
g_n <- t_n*ell_n

u_n <- log(l_n) + log(alpha + c_n) + log(alpha + g_n)

### RE deviation
u_dev <- log(1-ell_r) + log(alpha + (1-(1/2))*ell_r) + log(alpha + (1/2)*ell_r)


### calculate the value difference between staying on the equilibrium and
### deviating from it at each possible value of alpha and beta

v <- matrix( nrow = n, ncol = n )

c <- 0
for ( b in beta ){
  c <- c + 1
    
  ### calculate benefit to government to staying on equilibrium path
  v[,c] <- (b/(1-b))*(u_r - u_n) - (u_dev-u_r)
}

bmin <- rep(0,times=n)
for ( i in 1:n ){
  bmin[i] <- min(beta[v[i,]>=0])
}

### Plot RE allocation at each alpha

figure1 <-ggplot(data = data.frame(alpha = c(alpha,alpha,alpha,alpha),
                                   value = c(t_r,c_r,l_r,g_r),
                                   param = c(
                                     rep('tax rate',times = n),
                                     rep('consumption', times = n),
                                     rep('leisure',times = n),
                                     rep('public good', times = n))
                                   )
                 ) +
  geom_line(aes(alpha,value,group=param,color=param)) + 
  labs( title = 'Ramsey equilibrium allocation at each \u03b1',
        y = element_blank(), x = expression(alpha),
        color = 'Parameter' ) + theme_classic() + 
  scale_y_continuous( limits = c(0,1) ) 

png(file='figure1.png', width = 350, height = 250)
figure1
dev.off()


### plot comparisons between each of the allocations
dat <- rbind(
  data.frame(
    alpha = c(alpha,alpha,alpha),
    value = c(l_n,l_r,l_sp),
    param = rep('Leisure',times = n),
    equilibrium = c(rep('NE', times = n),
                    rep('Ramsey', times = n),
                    rep('Social Planner', times = n)
                    )
    ),
    data.frame(
    alpha = c(alpha,alpha,alpha),
    value = c(c_n,c_r,c_sp),
    param = rep('Consumption',times = n),
    equilibrium = c(rep('NE', times = n),
                    rep('Ramsey', times = n),
                    rep('Social Planner', times = n)
    )
  ),
  data.frame(
    alpha = c(alpha,alpha,alpha),
    value = c(g_n,g_r,g_sp),
    param = rep('Public Good',times = n),
    equilibrium = c(rep('NE', times = n),
                    rep('Ramsey', times = n),
                    rep('Social Planner', times = n)
    )
  )
)

figure2 <- ggplot(dat,aes(alpha,value,group=equilibrium,color=equilibrium)) + 
  labs( title = 'Allocations by Type',
        x = expression(alpha),
        y = element_blank() ) +
  scale_x_continuous(limits = c(0,0.5), breaks = c(0,0.25,0.5)) +
  geom_line() + facet_wrap(~ param) + theme_classic()





png(file='figure2.png', width = 400, height = 150)
  figure2
dev.off()



### plot value-beta relationship at the bounds and midpoint of alpha
### (note: upper bound graphed separately because of scale issues)

figure3 <-ggplot(data = data.frame(beta = bmin,alpha)[2:n,]) +
  geom_line(aes(alpha,beta), color = 'red') + 
  labs( title = 'Minimum \u03b2 at each \u03b1',
        subtitle = 'Required to support Ramsey equilibrium',
        y = '\u03b2', x = '\u03b1' ) + theme_classic()

png(file='figure3.png', width = 350, height = 300)
  figure3
dev.off()
