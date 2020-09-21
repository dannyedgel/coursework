###
### This file is used to generate charts for Econ 711 problem set #2
###
### Econ 711: Microeconomics I
### Fall 2020
###
### Date created:  20 Sep 2020
### Last modified: 20 Sep 2020
### Author: Danny Edgel
###

### clear workspace
rm(list=ls())

### declare package libraries
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ711/problem_sets/PS2/')


###
### Question 2
###____________________________________________

### define constants
B     <- 1
dim   <- 100

### define variables
y1     <- -dim:0
y2     <- B*(-y1)^(2/3)
poly.y <- rep(-2,times = dim+1)

### generate plot
dat <- data.frame(y1,y2)
p <- ggplot(dat,aes(y1,y2)) +
  geom_line(color='blue',size=1.5) +
  geom_polygon(data = data.frame(x1 = c(y1,y1[dim+1],0:(-dim),y1[1]),
                                 x2 = c(y2,-2,poly.y,y2[1]) ), 
               aes(x1,x2), fill = 'steelblue2' ) +
  geom_vline(xintercept=0, size = 1) +
  geom_hline(yintercept=0, size = 1) +
  ylab(expression(y[2])) +
  xlab(expression(y[1])) +
  scale_x_continuous(lim = c(-dim,1), expand = c(0, 0)) +
  scale_y_continuous(lim = c(-2.0001,max(y2)), position='right',
                     expand = c(0, 0)) +
  theme(axis.ticks = element_blank(),
        axis.text  = element_blank(),
        panel.background = element_blank()
        )
p


png(file='problem2a_prodset.png', width = 350, height = 170)
p
dev.off()
