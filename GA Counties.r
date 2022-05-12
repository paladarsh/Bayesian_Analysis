rm(list=ls())
library(rjags)
library(dplyr)
library(LearnBayes)
library(ggplot2)


#-------------- read data prepared in external programs
TFR.df = read.csv('GA Counties 2010 Census data.csv', stringsAsFactors = FALSE)

NCHS   = read.csv('NCHS fertility data.csv', stringsAsFactors = FALSE)
Oasis  = read.csv('GA Oasis data.csv', stringsAsFactors = FALSE)
load(file='svd.constants.RData')
m = svd.constants$m
X = svd.constants$X
#------------- MORTALITY model ------------------
q5_hat = Oasis$q5  # estimates for 2006-2014 from GA public data website

# calculate a and b coeffs for each q5 (2 x 159)
ab = sapply(q5_hat, function(this.q) {
  LearnBayes::beta.select( list(x= this.q/2, p=.05), list(x=this.q*2, p=.95))
})

q5_a = ab[1,]
q5_b = ab[2,]

##--- Wilmoth et al. coefficients from Pop Studies
wilmoth = 
  read.csv(text = '
           age,am,bm,cm,vm,af,bf,cf,vf
           0,  -0.5101, 0.8164,-0.0245,     0,-0.6619, 0.7684,-0.0277,     0
           1,      -99,    -99,    -99,   -99,    -99,    -99,    -99,   -99
           5,  -3.0435, 1.5270, 0.0817,0.1720,-2.5608, 1.7937, 0.1082,0.2788
           10, -3.9554, 1.2390, 0.0638,0.1683,-3.2435, 1.6653, 0.1088,0.3423
           15, -3.9374, 1.0425, 0.0750,0.2161,-3.1099, 1.5797, 0.1147,0.4007
           20, -3.4165, 1.1651, 0.0945,0.3022,-2.9789, 1.5053, 0.1011,0.4133
           25, -3.4237, 1.1444, 0.0905,0.3624,-3.0185, 1.3729, 0.0815,0.3884
           30, -3.4438, 1.0682, 0.0814,0.3848,-3.0201, 1.2879, 0.0778,0.3391
           35, -3.4198, 0.9620, 0.0714,0.3779,-3.1487, 1.1071, 0.0637,0.2829
           40, -3.3829, 0.8337, 0.0609,0.3530,-3.2690, 0.9339, 0.0533,0.2246
           45, -3.4456, 0.6039, 0.0362,0.3060,-3.5202, 0.6642, 0.0289,0.1774
           50, -3.4217, 0.4001, 0.0138,0.2564,-3.4076, 0.5556, 0.0208,0.1429
           55, -3.4144, 0.1760,-0.0128,0.2017,-3.2587, 0.4461, 0.0101,0.1190
           60, -3.1402, 0.0921,-0.0216,0.1616,-2.8907, 0.3988, 0.0042,0.0807
           65, -2.8565, 0.0217,-0.0283,0.1216,-2.6608, 0.2591,-0.0135,0.0571
           70, -2.4114, 0.0388,-0.0235,0.0864,-2.2949, 0.1759,-0.0229,0.0295
           75, -2.0411, 0.0093,-0.0252,0.0537,-2.0414, 0.0481,-0.0354,0.0114
           80, -1.6456, 0.0085,-0.0221,0.0316,-1.7308,-0.0064,-0.0347,0.0033
           85, -1.3203,-0.0183,-0.0219,0.0061,-1.4473,-0.0531,-0.0327,0.0040
           90, -1.0368,-0.0314,-0.0184,     0,-1.1582,-0.0617,-0.0259,     0
           95, -0.7310,-0.0170,-0.0133,     0,-0.8655,-0.0598,-0.0198,     0
           100,-0.5024,-0.0081,-0.0086,     0,-0.6294,-0.0513,-0.0134,     0
           105,-0.3275,-0.0001,-0.0048,     0,-0.4282,-0.0341,-0.0075,     0
           110,-0.2212,-0.0028,-0.0027,     0,-0.2966,-0.0229,-0.0041,     0
           ')

af = wilmoth$af[1:11]  # keep age 0,1,...45 
bf = wilmoth$bf[1:11]  # keep age 0,1,...45 
cf = wilmoth$cf[1:11]  # keep age 0,1,...45 
vf = wilmoth$vf[1:11]  # keep age 0,1,...45 
ncounty = nrow(TFR.df)  
## construct the matrix of constants for cumulative hazard calculations
n = c(1,5, rep(5,9))    # widths of life table age intervals for x=0,1,5,10...45
cs_constants = matrix(0, 11, 12)
for (j in 1:11) cs_constants[1:j,j+1] = head(n,j)
cs_constants

## construct the constants for the trapezoidal approx of L0...L45 from a row of l0,l1,l5,l10,...,l50
trapez_constants = matrix(0, 12, 10, dimnames=list(paste0('l', c(0,1,seq(5,50,5))), paste0('L', seq(0,45,5))))
trapez_constants[c('l0','l1','l5'), 'L0'] = c( 1/2, 5/2, 4/2)
for (j in 2:10) trapez_constants[j+1:2, j] = 5/2

model_string <- textConnection("model{
# Likelihood
  gamma =  m + beta%*%X ;

  for(i in 1:n){
    for(j in 1:7){
      expg[i,j]=exp(gamma[i,j]);
    }
    colsums_phi[i] = sum(expg[i,1:7]);
  }
  for(i in 1:n){
    for(j in 1:7){
      phi[i,j]=(1/colsums_phi[i])*expg[i,j];
    }
  }


  h = log(q5);
  for (i in 1:n){
    tp[i,1:11]=af+h[i]*bf+h[i]*h[i]*cf+k[i]*vf;
    for (j in 1:11){
      tmp[i,j]=exp(tp[i,j]);
    }
  }
  mx[1:n,1]=tmp[1:n,1];
  mx[1:n,2]=-0.25*tmp[1:n,1] - 0.25 * log(1-q5);
  for(i in 3:11){
    mx[1:n,i]=tmp[1:n,i];
  }
  tmplx=-mx%*%cs_constants;
  for (i in 1:n){
    for (j in 1:12){
      lx[i,j]=exp(tmplx[i,j]);
    }
  }
  Lx=lx%*%trapez_constants;
  for (j in 1:7){
    Sx[1:n,j]=Lx[1:n,j+2]/Lx[1:n,j+3];
  }
  for (i in 1:n){
    Fx[i,1]=0.0;
    for (j in 1:7){
      Fx[i,j+1]=0.20*TFR[i]*phi[i,j];
    }
  }
  for (j in 1:7) {
    Kx[1:n,j] = (Sx[1:n,j] * Fx[1:n,j] + Fx[1:n,j+1]) * Lx[1:n,1] / 2 ; 
  }
  for (i in 1:n){
    expected_C[i] = sum(W[i,1:7]*Kx[i,1:7]);
  }

  # PRIORS
  for (i in 1:n){
    C[i] ~ dpois(expected_C[i]);
    beta[i, 1]  ~ dnorm(0,1); 
    beta[i, 2]  ~ dnorm(0,1); 
    
    TFR[i] ~ dnorm(mu, 1/sig_sq);  
    
    q5[i] ~ dbeta(q5_a[i], q5_b[i]);
    
    k[i] ~ dnorm(0,1);
  }

  }"
);

burn <- 100
thin <- 4
n.chains <- 2


C=TFR.df$C
W = as.matrix(select(TFR.df, contains('W')))
mu=runif(1,0,4)
sig_sq=runif(1,0,2)
m = matrix(m, nrow=ncounty, ncol=7, byrow=TRUE)
data=list(n=ncounty,C=C,W=W, q5_a=q5_a, q5_b=q5_b, af=af,bf=bf,cf=cf,vf=vf,cs_constants=cs_constants,trapez_constants=trapez_constants,m=m,X=X,mu=mu,sig_sq=sig_sq);
inits <- list(q5 = rbeta(n=ncounty, shape1=q5_a, shape2=q5_b),
              k = rnorm(n=ncounty, mean=0,sd=1),
              beta = matrix(runif(ncounty*2, min=-.10, max=.10) , ncounty,2),
              TFR = pmax( .10, TFR.df$iTFR + rnorm(ncounty, 0, sd=.10))
);
params=c("TFR")


model <- jags.model(model_string, data = data, inits=inits, n.chains = n.chains, quiet = FALSE);
update(model, burn, progress.bar = "text")
samples <- coda.samples(model, variable.names = params, n.iter = 600, progress.bar = "text")


summ=summary(samples)
Q=as.data.frame(summ$quantiles)

names(Q)=c('Q2.5','Q25','Q50','Q75','Q97.5')
Q$index  = rank(Q$Q50)
Q$iTFR   = TFR.df$iTFR
ggplot(data=Q, aes(x=Q50, y=index)) + 
  xlim(0.5,3.5) + xlab('County TFR') + ylab('Index of County') +
  geom_segment(aes(x=Q2.5, xend=Q97.5, y=index, yend=index), col='grey') +
  geom_point(size=1) +
  geom_point(aes(x=iTFR, y=index), shape='x',size=3, col='red') +
  theme_bw() 

save(model,samples,summ,Q, file=paste0('xyz.RData'))
