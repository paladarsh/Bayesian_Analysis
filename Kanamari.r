rm(list=ls())
library(rjags)
library(ggplot2)
C=191
W = c(40, 34, 29, 19, 14, 9, 8)
q5_a = 3.99
q5_b = 114.26
af = c(-0.6619, -99, -2.5608, -3.2435, -3.1099, -2.9789, -3.0185, -3.0201, -3.1487, -3.2690, -3.5202)
bf = c(0.7684,  -99,  1.7937,  1.6653,  1.5797, 1.5053,   1.3729, 1.2879, 1.1071, 0.9339, 0.6642)
cf = c(-0.0277, -99,  0.1082,  0.1088,  0.1147, 0.1011,  0.0815, 0.0778, 0.0637, 0.0533, 0.0289)
vf = c(      0, -99,  0.2788,  0.3423,  0.4007, 0.4133, 0.3884,  0.3391, 0.2829, 0.2246, 0.1774)
# fertility model constants
m = c(0, 1.3921, 1.5924, 1.2291, 0.4523, -0.8868, -3.4372)
X = matrix(c(   0, 0.2744, 0.5408, 0.7319, 0.8774,  1.0405,  1.5172, 
                0, 0.3178, 0.5108, 0.5107, 0.3483,  0.0529, -0.7236),
            nrow=7,ncol=2)
data=list(C=C,W=W, q5_a=q5_a, q5_b=q5_b, af=af, bf=bf, cf=cf, vf=vf,m=m,X=X)
# X(7*2),beta()
model_string <- textConnection("model{
# Likelihood
    h = log(q5);
    gamma = m+X%*%beta;
    for (i in 1:7){
      tot[i]=exp(gamma[i]);
    }
    for (i in 1:7){
      phi[i] = exp(gamma[i]) / sum(tot);
    }
    Fx[1] = 0;                                                # F10
    for (i in 2:8){
      Fx[i]   = TFR * phi[i-1] / 5;              # F15...F45 
    }
    mx[1]  = exp(af[1] + bf[1]*h + cf[1]*(h*h) + vf[1]*k);
    mx[2]  = -0.25 * (mx[1] + log(1-q5));               #  recalculate 1_mu_4 = -1/4 log(l[5]/l[1])

    for (i in 3:11) {
      mx[i]  = exp(af[i] + bf[i]*h + cf[i]*(h*h) + vf[i]*k);
    }

    lx[1]  = 1;                                          # x=0
    lx[2]  = lx[1] * exp(-mx[1]);                        # x=1
    for (i in 3:12){
      lx[i]  = lx[i-1] * exp(-5*mx[i-1]);  # x=5,10,...50
    }
  
    Lx[1]= 1* (lx[1]+lx[2])/2 + 4*(lx[2]+lx[3])/2 ;   # 5L0
    for (i in 2:10){
      Lx[i]  = 5* (lx[i+1]+lx[1+2])/2 ;                   # 5Lx
    }
  
    for (i in 1:7){
      Kx[i]  = (Lx[i+2]/Lx[i+3] * Fx[i] + Fx[i+1]) * Lx[1]/2 ;
    }
    Kstar  = inprod(W, Kx);
    C ~ dpois(Kstar);

  # PRIORS
  beta[1]  ~ dnorm(0,1); 
  beta[2]  ~ dnorm(0,1); 
  
  q5 ~ dbeta(q5_a, q5_b);   # 90% prior prob. that q5 is between 1/2 and 2x estimated q5
  TFR ~ dunif(1,11);
    
  k ~ dnorm(0,1);
  }"
  );

burn <- 100
n.iter <- 10100
thin <- 4
n.chains <- 2
inits <- list(q5 = rbeta(n=1, shape1=3.99, shape2=114.26),
              k = rnorm(n=1, mean=0,sd=1),
              beta = rnorm(2),
              TFR = runif(n=1,min=1, max=5)
);
model <- jags.model(model_string, data = data, inits=inits, n.chains = n.chains, quiet = TRUE);
update(model, burn, progress.bar = "none")
params=c("TFR")
samples <- coda.samples(model, variable.names = params, n.iter = n.iter, progress.bar = "none")
summary(samples)
plot(samples)
acf(samples[[1]])
