# investigate priors

# simple prior investigation using dens+HPDI
# looking at outcome scale using generative model

############################
## R0
############################
rnorm(n=10000,1.5,1) %>% dens(show.HPDI = c(0.95))

############################
## sigma
############################
rexp(n=10000, 10) %>% dens()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Dirichlet priors ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# dirichlet prior on levelled predictor importance
library(gtools)
set.seed(1805)
delta <- rdirichlet( 10 , alpha=rep(2,3) )
str(delta)
h <- 3
plot( NULL , xlim=c(1,3) , ylim=c(0,0.9) , xlab="index" , ylab="probability" )
for ( i in 1:nrow(delta) ) lines( 1:3 , delta[i,] , type="b" ,
                                  pch=ifelse(i==h,16,1) , lwd=ifelse(i==h,4,1.5) ,
                                  col=ifelse(i==h,"black",col.alpha("black",0.7)) )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Understand prior on normal through sigma using exp  ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rexp(10,rate=10) -> mse
rnorm(n = 2000,0.80,mse[1]) %>% dens(show.HPDI = 0.9)
rnorm(n = 2000,0.80,mse[2]) %>% dens(show.HPDI = 0.9,add=T)
rnorm(n = 2000,0.80,mse[3]) %>% dens(show.HPDI = 0.9,add=T)
rnorm(n = 2000,0.80,mse[4]) %>% dens(show.HPDI = 0.9,add=T)
rnorm(n = 2000,0.80,mse[5]) %>% dens(show.HPDI = 0.9,add=T)
