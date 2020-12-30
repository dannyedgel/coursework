
library(ggplot2)

pinv <- as.matrix(data.frame(c(1,-4),c(-2,7)))
p    <- as.matrix(data.frame(c(-7,-4),c(-2,-1)))
w    <- as.matrix(data.frame(c(1,2),c(-1,3)))
p
pinv
w
t <- p %*% w %*% pinv
x <- as.matrix(data.frame(c(1,2)))
t %*% x
t %*% p %*% x
pinv %*% t %*% p %*% x
w %*% x


dataset <- c(1,2)

for ( d in dataset ){
  p1 <- c(7,4)
  p2 <- c(5,5)
  p3 <- c(4,8)
  
  if ( d == 1 ){
    y1 <- c(-20,40)
    y2 <- c(-50,60)
    y3 <- c(-70,90)
  }else{
    y1 <- c(-20,40)
    y2 <- c(-40,70)
    y3 <- c(-70,90)
  }
  
  
  p <- list(p1,p2,p3)
  y <- list(y1,y2,y3)
  
  pi <- matrix(nrow=3,ncol=3)
  
  for ( r in 1:3 ){
    for ( c in 1:3 ){
      pi[r,c] <- p[[r]] %*% y[[c]]
    }
  }
  
  print(paste0('Dataset ',d,':'))
  print(pi)
  print('')
}

y[[4]] <- c(0,0)

dat <- data.frame(y1=rep(0,times=length(y)),y2=rep(0,times=length(y)))

for ( x in 1:length(y) ){
  dat$y1[x] <- y[[x]][1]
  dat$y2[x] <- y[[x]][2]
}


dat <- rbind(
  data.frame(y1=c(-80,-80),y2=c(-25,max(dat$y2))),
  data.frame(y1=c(-70,-40,-20,0),y2=c(90,70,40,0)),
  data.frame(y1=c(0),y2=c(-25))
)

filepath <- 'C:/Users/edgel/Google Drive/UW-Madison/Fall 2020/Econ 711 (Micro)/Problem Sets/PS1/'
png(file=paste0(filepath,'problem2b_prodset.png'), width = 350, height = 170)
ggplot(dat,aes(y1,y2)) +
  geom_line(dat[2:(nrow(dat)-1),],mapping=aes(y1,y2)) + 
  geom_point(dat[3:(nrow(dat)-1),],mapping=aes(y1,y2)) +
  geom_polygon(fill='blue',alpha=0.5) +
  theme_classic() + 
  scale_y_continuous(limits = c(-25,100), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-80,20), expand = c(0, 0)) +
  ylab('Output') + xlab('Input') +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0)
dev.off()


dat <- rbind(
  data.frame(y1=c(-80,-80),y2=c(-25,290/3)),
  data.frame(y1=c(-70,-40,-20),y2=c(90,70,40)),
  data.frame(y1=c(70/3),y2=c(-25))
)

png(file=paste0(filepath,'problem2c_prodset.png'), width = 350, height = 170)
ggplot(dat,aes(y1,y2)) +
  geom_line(dat[2:nrow(dat),],mapping=aes(y1,y2)) + 
  geom_point(dat[3:(nrow(dat)-1),],mapping=aes(y1,y2)) +
  geom_polygon(fill='blue',alpha=0.5) +
  theme_classic() + 
  scale_y_continuous(limits = c(-25,110), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-80,30), expand = c(0, 0)) +
  ylab('Output') + xlab('Input') +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0)
dev.off()

