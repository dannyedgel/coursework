# Pavel Brendler, Anton Babkin, 10/9/15.
# Modified by Jason Choi, 10/18/18.

# Translated to R by Danny Edgel, 10/06/2020

# This program evaluates macroeconomic consequences of eliminating social security in the U.S.
# The model is a simplified version of the model by Conesa and Krueger (1999).

# You are asked to understand the logic of the code, complete the missing parts
# and use it to run a policy experiment.
# Look for these comments: # <------------------- YOUR CODE GOES HERE

# This code is not written very efficiently, so it is not fast.
# Instead, it is very explicit, so that you can follow and understand it easier.


# Clear the memory
rm( list = ls() )

# install packages for dot product and graphing
library(pracma)
library(ggplot2)

### set working directory
setwd('C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/PS5')

### run model with and without social security
for ( ss in 1:2 ){
    
##########################################################################
# Parameters
##########################################################################

# Demographics
J  <- 66                      # life-span
JR <- 46                      # age of retirement
tR <- J-JR+1                  # length of retirement
tW <- JR-1                    # length of working life
n  <- 0.011                   # Population growth

# Preferences
beta  <- 0.97                  # discount factor
sigma <- 2                    # coefficient of relative risk aversion
gamma <- 0.42                 # weight on consumption

# Production
alpha <- 0.36                 # production elasticity of capital
delta <- 0.06                 # rate of depreciation

# Social Security tax rate
if ( ss == 1 ){ tau <- 0.11 } # Use this for model with SS
if ( ss == 2 ){ tau <- 0    } # Use this for model without SS


# Measure of each generation
mass <- rep(1, times = J)
for ( ik0 in 2:J ){
	mass[ik0]<-mass[ik0-1]/(1+n)
}

# Normalized measure of each generation (sums up to 1)
mass <- mass / sum(mass)

# Age-efficiency profile
e <- rep(0, times = tW)
e[1:11]  <- seq(from = 0.6,  to = 1,    length.out = 11)
e[11:21] <- seq(from = 1,    to = 1.08, length.out = 11)
e[21:31] <- seq(from = 1.08, to = 1.12, length.out = 11)
e[31:41] <- seq(from = 1.12, to = 1.06, length.out = 11)
e[41:45] <- seq(from = 1.06, to = 1.02, length.out = 5)

#  Capital grid
#  Be careful when setting maxkap across experiments! This value should not be binding!
#  To see if it's binding, check the optimal capital decision kp(k) of
#  every cohort.
#  If you have problems with convergence, increasing nk might help.

maxkap <- 14                               # maximum value of capital grid
minkap <- 0.01                             # minimum value of capital grid
nk     <- 180                              # number of grid points
inckap <- ((maxkap-minkap)/(nk-1)^2)       # distance between points
aux    <- 1:nk
kap    <- (minkap+inckap*(aux-1)^2)        # capital grid
# This formula makes a non-uniform grid with more density at lower range.

neg <- -1e10                              # very small number


#########################################################################
# Initial guesses for interest rate, wages and pension benefits
# Comment out the appropriate lines!

# Social Security
if ( ss == 1 ){
    K0<-3.1392
    L0<-0.3496
}

# No Social Security
if ( ss == 2){
    K0<-3.9288
    L0<- 0.3663
}



################################################################
## Loop over capital and labor
################################################################

tolk   <- 1e-3           # Numerical tolerance for capital
tollab <- 1e-3           # Numerical tolerance for labor
nq     <- 100            # Max number of iterations

q  <- 0                  # Counter for iterations
K1 <- K0 + 10
L1 <- L0 + 10

cat('\nComputing equilibrium prices... \n')
while ( q < nq & ( abs(K1-K0) > tolk | abs(L1-L0) > tollab ) ){
    q <- q + 1

    cat('\nIteration ',q,' out of ',nq,'\n')

    # Prices
    r0 <- alpha*((L0/K0)^(1-alpha))-delta  # <--- MBK = MCk
    w0 <- (1-alpha)*((K0/L0)^alpha)        # <--- MBL = MCL

    # Pension benefit
    b <- (tau*w0*L0) / sum(mass[JR:J])     # <--- Government BC is satisfied

    ############################################################
    # BACKWARD INDUCTION
    ############################################################

    # Initialization
    v      <- matrix(0, nrow = nk, ncol = J) # value function of agents
    kapopt <- matrix(1, nrow = nk, ncol = J) # optimal savings of agents
    
    # (store INDEX of k' in capital grid, not k' itself!)
    labopt <-  matrix(1, nrow = nk, ncol = tW) # optimal labor supply


    # Retired households

    # Last period utility
    cons <- ((1 + r0) * kap) + b                    # last period consumption (vector!)
    util <- (cons^(1-sigma))/(1-sigma)       # last period utility (vector!)
    v[,J] <- util                            # last period indirect utility (vector!)

    for ( j in seq(J-1,tW+1,by = -1) ){ # age
        for ( ik0 in 1:nk ){        # assets today

            # Initialize right-hand side of Bellman equation
            vmin <- neg
            ik1  <- 0

            # Loop over all k's in the capital grid to find the value,
            # which gives max of the right-hand side of Bellman equation

            while ( ik1 < nk ){ 	# assets tomorrow
                ik1  <- ik1 + 1
                kap0 <- kap[ik0] # current asset holdings
                kap1 <- kap[ik1] # future asset holdings

                # Instantaneous utility
                cons <- b + (1+r0)*kap0 - kap1 # <--- budget constraint satisfied

                if ( cons <= 0 ){
                    util <- neg
                }else{
                    util <- (cons^(1-sigma))/(1-sigma) # <--- retired instantaneous utility
                }

                # Right-hand side of Bellman equation
                v0 <- util + beta*v[ik1,j+1] # <--- RHS of Bellman

                # Store indirect utility and optimal saving
                if ( v0 > vmin ){
                    v[ik0,j]      <- v0
                    kapopt[ik0,j] <- ik1
                    vmin          <- v0
                }
            }
        }
    }

    # Working households
    for ( j in seq( tW, 1, by = -1 ) ){           # age
        for ( ik0 in 1:nk ){          # assets today

            # Initialize right-hand side of Bellman equation
            vmin <- neg
            ik1  <- 0

            # Loop over all k's in the capital grid to find the value,
            # which gives max of the right-hand side of Bellman equation

            while ( ik1 < nk ){ 	# assets tomorrow
                ik1 <- ik1+1

                kap0 <- kap[ik0] # current asset holdings
                kap1 <- kap[ik1] # future asset holdings

                # Optimal labor supply
                lab <- gamma - ((1-gamma)*(((1+r0)*kap0)-kap1)/((1-tau)*w0*e[j]))
                
                #(gamma*(1-tau)*e[j]*w0 - (1-gamma)*((1+r0)*kap0-kap1))/
                #    ((1-tau)*e[j]*w0) # <--- individual labor supply

                # Check feasibility of labor supply
                if ( lab > 1 ){
                    lab <- 1
                }else if ( lab < 0 ){
                    lab <- 0
                }

                # Instantaneous utility
                cons <- (1-tau)*w0*e[j]*lab+(1+r0)*kap0 - kap1 # <--- working-age BC

                if ( cons <= 0 ){
                    util <- neg
                }else{
                    util <- (((cons^gamma)*((1-lab)^(1-gamma)))^(1-sigma))/
                        (1-sigma) # <--- instantaneous utility in working age
                }

                # Right-hand side of Bellman equation
                v0 <- util + beta*v[ik1,j+1] # <--- utility plus discounted future utility

                # Store indirect utility, optimal saving and labor
                if ( v0 > vmin ){
                    v[ik0,j] <- v0
                    kapopt[ik0,j] <- ik1
                    labopt[ik0,j] <- lab
                    vmin <- v0
                }
            }
        }
    }


    

    ############################################################################
    # Aggregate capital stock and employment                                  #
    ############################################################################

    # Initializations
    kgen   <- rep(0, times = J)  # capital supply k(j+1) for each generation j
    ikgen  <- rep(0, times = J)  # grid index of k(j+1)
    labgen <- rep(0, times = tW) # labor supply l(j) for each generation j

    # Use decision rules to iteratively fill in kgen and labgen
    ik0 <- 1                # starting capital of j <- 1, kap(ik0) <- 0
    for ( j in 1:J ){               # iterations over cohort
        
        # capital decision kp(k)
        ik1      <- kapopt[ik0, j]
        ikgen[j] <- ik1
        kgen[j]  <- kap[ik1]

        # labor decision l(k)
        if ( j <= tW ){
            labgen[j] <- labopt[ik0, j]
        }

        # update k <- kp
        ik0 <- ik1
    }

    K1 <- dot(kgen,mass)                     # dot product of vectors
    L1 <- dot((labgen * e),mass[1:tW])      # dot product of vectors

    # Update the guess on capital and labor
    K0 <- 0.9*K0 + 0.1*K1
    L0 <- 0.9*L0 + 0.1*L1

    # Display results
    cat('\nCapital: ',  K0,
        '\nLabor: ',    L0,
        '\nPension: ',  b  ,
        '\n\nCapital deviation: ', abs(K1-K0),
        '\nLabor deviation: ',     abs(L1-L0),
        '\n'
        )
} ## end while loop

# Display equilibrium results
cat(
    'K0: ', K0, '\n',
    'L0: ', L0, '\n',
    'w:  ', w0, '\n',
    'r:  ', r0, '\n',
    'b:  ', b
)

# Check if any cohort wants to save on the upper bound of capital grid
if ( sum(kgen == kap[length(kap)]) > 0 ){
    print('Capital decision on upper bound; increase!')
}


## Plots
# Value function for a retired agent
age <- 50
dat <- data.frame(kap,value = v[,age])

if ( tau !=  0 ){
    filename <- 'value_ss.png'
}else if ( tau == 0 ){
    filename <- 'value_noss.png'
}

png(file=filename, width = 350, height = 350)
print({
    ggplot(dat) +
        geom_line(aes(kap,value), color = 'blue') +
        xlab('asset holdings, k') + 
        ylab(expression(paste('value function, ',V[50],'(k)'))) +
        scale_x_continuous(breaks = seq(0,14,by=2)) +
        ggtitle(paste0('value function of a retired agent at age ',age)) +
        theme_classic()
})
dev.off()
# figure(1)
# age = 50
# plot(kap,v(:,age))
# xlabel('asset holdings, k','FontSize',14)
# ylabel('value function, V_{50}(k)','FontSize',14)
# title(['value function of a retired agent at age ', num2str(age)],'FontSize',14)


# Savings of a working agent
age <- 20
dat <- data.frame(kap,savings = kap[kapopt[,age]], type = 'Saving')
dat <- rbind(
    dat,
    data.frame(kap, savings = kap, type = '45 deg. line')
)

if ( tau !=  0 ){
    filename <- 'savingk_ss.png'
}else if ( tau == 0 ){
    filename <- 'savingk_noss.png'
}

png(file=filename, width = 425, height = 350)
print({
    ggplot(dat) +
    geom_line(aes(kap,savings, group = type, color = type, linetype = type)) +
    ggtitle(paste0("saving k' of a working agent at age ",age))+
    xlab('asset holdings, k') +  ylab("saving, k'") +
    scale_linetype_manual(values=c("dashed","solid")) +
    scale_color_manual(values = c('red','blue')) +
    scale_x_continuous(breaks = seq(0,14,by=2)) +
    scale_y_continuous(breaks = seq(0,14,by=2)) +
    guides(size = 'legend', nrow = 2, label = T) + theme_classic() +
    theme(legend.title = element_blank())
})

dev.off()

# figure(2)
# age = 20
# plot(kap,kap(kapopt(:,age)),'k-',kap,kap,'r--')
# xlabel('asset holdings, k','FontSize',14)
# ylabel('saving, k''','FontSize',14)
# legend('saving','45 degree line','Location','East')
# title(['saving k'' of a working agent at age ', num2str(age)],'FontSize',14)


## Saving Model Results
if ( tau !=  0 ){
    # Solve the model with social security and
    save.image('ss.RData')
}else if ( tau == 0 ){
    # Solve the model without social security and
    save.image('no_ss.RData')
}


} ## end inital SS/no-SS loop

#############################################################################################
## Compare Model Results ## RUN THIS PART OF THE CODE AFTER YOU'VE SOLVED BOTH SS AND WITHOUT SS
#############################################################################################
rm( list = ls() )

load('ss.RData')
kgen_ss <- kgen

load('no_ss.RData')
kgen_noss <- kgen

dat <- data.frame(x = (20 + 1):(20 + J), kgen = kgen_ss, 
                  type = 'With SS' )
dat <- rbind(
    dat,
    data.frame(x = (20 + 1):(20 + J), kgen = kgen_noss, 
               type = 'Without SS')
)

png(file='wealth.png', width = 450, height = 350)
ggplot(dat) +
    geom_line(aes(x,kgen, group = type, color = type)) + 
    xlab('(real-life) age') +
    ylab('wealth') +
    scale_linetype_manual(values=c("solid","solid")) +
    scale_color_manual(values = c('red','blue')) +
    scale_x_continuous(breaks = seq(21,86,by=4)) +
    scale_y_continuous(breaks = seq(0,12.5,by=1)) +
    guides(size = 'legend', nrow = 2, label = T) + theme_classic() +
    theme(legend.title = element_blank())
dev.off()

# figure(3)
# plot1 <- plot(20+1:20+J,kgen_ss,'b-',20+1:20+J,kgen_noss,'ro')
# set(plot1,'LineWidth',1.5)
# xlabel('(real-life) age','FontSize',14)
# ylabel('wealth','FontSize',14)
# AX <- legend('with Social Security','without Social Security','Location','NorthWest')
# LEG <- findobj(AX,'type','text')
# set(LEG,'FontSize',14)
# 
# # save figure to pdf
# set(3, 'PaperSize', [5 5])
# set(3, 'PaperPositionMode', 'manual')
# set(3, 'PaperPosition', [0 0 5 5])
# print(3, 'fig_wealth', '-dpdf')

## Welfare comparison
rm( list = ls() )

load('ss.Rdata')

# Newborn generation
V1_ss <- v[1,1]

# Aggregate welfare
V_ss <- rep(0,times = J)
V_ss[1] <- v[1]
for ( j in 2:J ){
    ik0 <- ikgen[j-1]
    V_ss[j] <- v[ik0,j]
}
    
W_ss <- dot(V_ss,mass)


load('no_ss.RData')

# Newborn generation
V1_noss <- v[1,1]

# Aggregate welfare
V_noss    <- rep(0, times = J)
V_noss[1] <- v[1,1]
for ( j in 2:J ){
    ik0       <- ikgen[j-1]
    V_noss[j] <- v[ik0,j]
}
    
W_noss <- dot(V_noss,mass)

### display each of the values calculated above
cat(
    '\nNewborn generation: ',
    '\n\nV1 w/ SS: ', V1_ss,
    '\nW w/ SS:  ', W_ss,
    '\nV1 w/o SS:', V1_noss,
    '\nW w/o SS: ', W_noss
)
