###
### This file is used to generate charts for Econ 712 problem set #3
###
### Econ 712: Macroeconomics I
### Fall 2020
###
### Date created:  22 Sep 2020
### Last modified: 22 Sep 2020
### Author: Danny Edgel
###

### clear workspace
rm(list=ls())

### declare package libraries
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/PS3/')


###
### Question 2a
###____________________________________________

### define constants
w1 <- 0  # good 1 endowment
w2 <- 2  # good 2 endowment

max_bc_len <- 100
max_p_len  <- 1000

## generate data frame of all possible price ratios
all_p <- seq(from=1/100,to=20,by=0.01)

### define utility function
Uf <- function(c1,c2) u =10*c1-4*c1^2+4*c2-c2^2

## initialize and fill a list of budget constraints
bc.list <- list()

for ( i in 1:length(all_p) ){
  
  print(paste0(
    'Generating budget constraint for price ratio ',i,
    ' of ',length(all_p)
  ))
  
  p <- all_p[i]
  
  ## fill in consumption of each good 
  u <- c1 <- c2 <- c()
  k <- 1 # j loop counter
  
  ## from (w1,w2) to (0,c2_max)
  for ( j in seq(w1,0,length.out = max_bc_len) ){
    c1 <- c(c1,j)
    c2 <- c(
      c2,
      w2 + p*w1 - p*c1[k]
    )
    u <- c(
      u,
      10*c1[k] - 4*(c1[k]^2) + 4*c2[k] - c2[k]^2
    )
    
    k <- k + 1
  } ## end first j loop
  
  ## from (w1,w2) to (c1_max,0)
  for ( j in seq(w2,0,length.out = max_bc_len) ){
    c2 <- c(c2,j)
    c1 <- c(
      c1,
      w1 + (1/p)*w2 - (1/p)*c2[k]
    )
    u <- c(
      u,
      10*c1[k] - 4*(c1[k]^2) + 4*c2[k] - c2[k]^2
    )
    
    k <- k + 1
  } ## end second j loop
  
  ### sort data on c1 (the x variable)
  p.dat <- unique(data.frame(c1,c2,u))
  p.dat <- p.dat[order(p.dat$c1),]
  
  ### add current price ratio, consumption, and utility to list
  bc.list[[as.character(p)]] <- p.dat
  
} ## end i loop


### generate a few indifference curves
util <- c(0.5,1,10)
u.df <- data.frame(c1=seq(0,6,length=100),c2=seq(0,6, length=100))
u  <- as.data.frame(outer(u.df$c1,u.df$c2,function(c1,c2)Uf(c1,c2))) %>% #convert the matrix to data frame
  rownames_to_column() %>% #get row coordinates
  gather(key, value, -rowname) %>% #convert to long format
  mutate(key = as.numeric(gsub("V", "", key))/10, #convert the column names to numbers
         rowname = as.numeric(rowname)/10)
# 
# u  <- outer(u.df$c1,u.df$c2,function(c1,c2)Uf(c1,c2))
# contour(u.df$c1,u.df$c2,u,levels=util) + lines(bc.list[[10]]$c1,bc.list[[10]]$c2)
# 
#u.df <- u.df[u.df$u %in% c(0.5, 1, 4),]


### export plot
png(file='2a_IC-and-BC.png', width = 350, height = 350)

ggplot() + 
  geom_line(data = bc.list[['0.5']],aes(c1,c2), color= 'firebrick1', size = .75) +
  geom_line(data = bc.list[['1']],aes(c1,c2), color= 'firebrick4', size = .75) +
  geom_contour(data = u,aes(x=rowname,y=key,z=value), breaks = c(8,9,7)) + 
  coord_cartesian(ylim=c(0,3), xlim = c(0,3)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab(expression(c[2])) + xlab(expression(c[1])) +
  ggtitle('Indifference curves (blue) vs. budget constraints (reds)') 
  
dev.off()


### solve for optimal allocation with each price ratio
c1 <- c()
c2 <- c()
for ( x in bc.list ){
  c1 <- c(
    c1,
    x$c1[x$u == max(x$u)]
  )
  c2 <- c(
    c2,
    x$c2[x$u == max(x$u)]
  )
} ## end i loop

## derive x1 and x2 from c1 and c2
x1 <- c1 - w1
x2 <- c2 - w2

dat <- data.frame(x1,x2)
png(file='2a_offer-curve.png', width = 350, height = 350)

ggplot(dat[dat$x2<=0,]) + 
  stat_smooth(aes(x1,x2), color= 'firebrick1', size = .75) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim=c(-2,0.5), xlim = c(-0.5,2)) +
  ylab(expression(x[2])) + xlab(expression(x[1])) +
  ggtitle('Offer Curve') +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0)

dev.off()


###
### Question 2b
###____________________________________________
### define constants
w1 <- 1  # good 1 endowment
w2 <- 0  # good 2 endowment

max_bc_len <- 100
max_p_len  <- 1000

## generate data frame of all possible price ratios
all_p <- seq(from=1/100,to=5,by=0.01)

### define utility function
Uf <- function(c1,c2) apply(data.frame(x1=2*c1 + c2,x2=c1 + 2*c2),1,FUN=min)

## initialize and fill a list of budget constraints
bc.list <- list()

for ( i in 1:length(all_p) ){
  
  print(paste0(
    'Generating budget constraint for price ratio ',i,
    ' of ',length(all_p)
  ))
  
  p <- all_p[i]
  
  ## fill in consumption of each good 
  u <- c1 <- c2 <- c()
  k <- 1 # j loop counter
  
  ## from (w1,w2) to (0,c2_max)
  for ( j in seq(w1,0,length.out = max_bc_len) ){
    c1 <- c(c1,j)
    c2 <- c(
      c2,
      w2 + p*w1 - p*c1[k]
    )
    u <- c(
      u,
      Uf(c1[k],c2[k])
    )
    
    k <- k + 1
  } ## end first j loop
  
  ## from (w1,w2) to (c1_max,0)
  for ( j in seq(w2,0,length.out = max_bc_len) ){
    c2 <- c(c2,j)
    c1 <- c(
      c1,
      w1 + (1/p)*w2 - (1/p)*c2[k]
    )
    u <- c(
      u,
      Uf(c1[k],c2[k])
    )
    
    k <- k + 1
  } ## end second j loop
  
  ### sort data on c1 (the x variable)
  p.dat <- unique(data.frame(c1,c2,u))
  p.dat <- p.dat[order(p.dat$c1),]
  
  ### add current price ratio, consumption, and utility to list
  bc.list[[as.character(p)]] <- p.dat
  
} ## end i loop


### generate a few indifference curves
util <- c(0.5,1,10)
u.df <- data.frame(c1=seq(0,6,length=100),c2=seq(0,6, length=100))
u  <- as.data.frame(outer(u.df$c1,u.df$c2,function(c1,c2)Uf(c1,c2))) %>% #convert the matrix to data frame
  rownames_to_column() %>% #get row coordinates
  gather(key, value, -rowname) %>% #convert to long format
  mutate(key = as.numeric(gsub("V", "", key))/10, #convert the column names to numbers
         rowname = as.numeric(rowname)/10)
# 
# u  <- outer(u.df$c1,u.df$c2,function(c1,c2)Uf(c1,c2))
# contour(u.df$c1,u.df$c2,u,levels=util) + lines(bc.list[[10]]$c1,bc.list[[10]]$c2)
# 
#u.df <- u.df[u.df$u %in% c(0.5, 1, 4),]
contour(outer(u.df$c1,u.df$c2,function(c1,c2)Uf(c1,c2)))

### export plot
png(file='2b_IC-and-BC.png', width = 350, height = 350)

ggplot() + 
  geom_line(data = bc.list[['0.5']],aes(c1,c2), color= 'firebrick1', size = .75) +
  geom_line(data = bc.list[['1']],aes(c1,c2), color= 'firebrick4', size = .75) +
  geom_contour(data = u,aes(x=rowname,y=key,z=value), breaks = c(0.75,1,2)) + 
  coord_cartesian(ylim=c(0,3), xlim = c(0,3)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab(expression(c[2])) + xlab(expression(c[1])) +
  ggtitle('Indifference curves (blue) vs. budget constraints (reds)') 

dev.off()


### solve for optimal allocation with each price ratio
c1 <- c()
c2 <- c()
for ( x in bc.list ){
  c1 <- c(
    c1,
    x$c1[x$u == max(x$u)]
  )
  c2 <- c(
    c2,
    x$c2[x$u == max(x$u)]
  )
} ## end i loop

## derive x1 and x2 from c1 and c2
x1 <- c1 - w1
x2 <- c2 - w2

dat <- unique(data.frame(x1,x2))

dat <- dat %>%
  mutate(x1 = round(x1,digits = 2),
         x2 = round(x2,digits = 2)
         )

dat <- unique(dat)
dat <- dat[order(-dat$x2),]
#dat <- dat %>% group_by(x1) %>% summarise(x2 = max(x2))


png(file='2b_offer-curve.png', width = 350, height = 350)

ggplot(dat) + 
  geom_path(aes(x1,x2), color= 'firebrick1', size = .75) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0),lim=c(0,3)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim=c(-0.5,2.5), xlim = c(-2.5,0.5)) +
  ylab(expression(x[2])) + xlab(expression(x[1])) +
  ggtitle('Offer Curve') +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0)

dev.off()



###
### Question 2c
###____________________________________________

### define constants
w1 <- 1  # good 1 endowment
w2 <- 10  # good 2 endowment

max_bc_len <- 100
max_p_len  <- 1000

## generate data frame of all possible price ratios
all_p <- seq(from=1/10,to=40,by=0.01)

### define utility function
Uf <- function(c1,c2) apply(data.frame(x1=2*c1 + c2,x2=c1 + 2*c2),1,FUN=min)

## initialize and fill a list of budget constraints
bc.list <- list()

for ( i in 1:length(all_p) ){
  
  print(paste0(
    'Generating budget constraint for price ratio ',i,
    ' of ',length(all_p)
  ))
  
  p <- all_p[i]
  
  ## fill in consumption of each good 
  u <- c1 <- c2 <- c()
  k <- 1 # j loop counter
  
  ## from (w1,w2) to (0,c2_max)
  for ( j in seq(w1,0,length.out = max_bc_len) ){
    c1 <- c(c1,j)
    c2 <- c(
      c2,
      w2 + p*w1 - p*c1[k]
    )
    u <- c(
      u,
      Uf(c1[k],c2[k])
    )
    
    k <- k + 1
  } ## end first j loop
  
  ## from (w1,w2) to (c1_max,0)
  for ( j in seq(w2,0,length.out = max_bc_len) ){
    c2 <- c(c2,j)
    c1 <- c(
      c1,
      w1 + (1/p)*w2 - (1/p)*c2[k]
    )
    u <- c(
      u,
      Uf(c1[k],c2[k])
    )
    
    k <- k + 1
  } ## end second j loop
  
  ### sort data on c1 (the x variable)
  p.dat <- unique(data.frame(c1,c2,u))
  p.dat <- p.dat[order(p.dat$c1),]
  
  ### add current price ratio, consumption, and utility to list
  bc.list[[as.character(p)]] <- p.dat
  
} ## end i loop



### solve for optimal allocation with each price ratio
c1 <- c()
c2 <- c()
for ( x in bc.list ){
  c1 <- c(
    c1,
    x$c1[x$u == max(x$u)]
  )
  c2 <- c(
    c2,
    x$c2[x$u == max(x$u)]
  )
} ## end i loop

## derive x1 and x2 from c1 and c2
x1 <- c1 - w1
x2 <- c2 - w2

dat <- unique(data.frame(x1,x2))

dat <- dat %>%
  mutate(x1 = round(x1,digits = 4),
         x2 = round(x2,digits = 4)
  )

dat <- unique(dat)
dat <- dat[order(-dat$x1),]

png(file='2c_offer-curve.png', width = 350, height = 350)

ggplot(dat) + 
  geom_line(aes(x1,x2), color= 'firebrick1', size = .75) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim=c(-12,22), xlim = c(-12,22)) +
  ylab(expression(x[2])) + xlab(expression(x[1])) +
  ggtitle('Offer Curve') +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0)

dev.off()
