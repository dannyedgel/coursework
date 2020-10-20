###
### This file is used to graph house prices in one period against those in the
### previous period for question two of Econ 712 PS7
###
### Date created:  20 October 2020
### Last modified: 20 October 2020
### Author: Danny Edgel
###

### clear workspace
rm( list = ls() )

### set current directory
setwd( 'C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/PS7')

### declare graphing packages
library(ggplot2)

### set parameters
alpha <- .75
beta  <- .75
y     <- 1


### define bounds of prices for graphing purposes
pmin <- 0.001
pmax <- y - .2  ## note: last element of pvec will not be used so that p<y
n    <- 1000    ## number of elements in pvec


### check that assumptions hold
1+alpha > beta*y


### define price vector
pvec <- seq(pmin, pmax, length.out = n+1)[1:n]

### define vector for next period's price
pnext <- rep(0, times = n)

for ( i in 1:n ){
  pt <- pvec[i]
  
  ### use law of motion to determine next period's price
  pnext[i] <- pt/(beta*(y-pt)) - alpha/beta
}

### graph the price in each period against the price in the next period
### (also add data for 45 degree line)
dat <- data.frame( pt    = c(pvec,pvec), 
                   pnext = c(pnext,pvec),
                   grp   = c(rep('motion',times = n),
                             rep('45 degree',times = n)
                             )
)

gr <- ggplot(dat) +
  geom_line(aes(pt,pnext, color = grp, group = grp,)) + 
  scale_color_manual(values = c('black','red'),
                     labels = expression(
                       45~degree,p[t+1]
                     )
                     ) +
  xlab(expression(p[t])) + ylab(expression(p[t+1])) +
  theme_classic() + theme(legend.title = element_blank())




png(file='priceplot.png', width = 320, height = 250)
  gr
dev.off()
