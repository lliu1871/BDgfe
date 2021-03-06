rm(list=ls())
N=30 #number of duplication times
m=1 # order of simulated tree among 100 trees 

###for sub  
#mu=0.8

mu=-0.1
alpha=0.6
beta<-0.9
gamma<-5
lambda<-0.3
time1=1e-10
current_time <- 10
start_ntaxa=2
end_ntaxa=end_ntaxa=start_ntaxa+N   
M=30 # the number of samples when simulating data

##### a sharp sub curve ##########
mu=-0.01
alpha=0.6
beta<-1.1
gamma<-1
lambda<-0.2
time1=1e-10
current_time <- 10
start_ntaxa=2
end_ntaxa=end_ntaxa=start_ntaxa+N   #define number of duplicates
#M=30 # the duplication

library(cubature)
#########################################################
#   subfunctionalization - age dependent
#########################################################
### first calculate the expectation of loss rate since the loss rate is depend on the age, which is also a random variable,
### then substitute the loss rate in the formulae with the expectation of loss rate

### 1. expectation of loss rate phi_t
### den_age: the density function of age
### phi_t: expectation of loss rate 

mu_t=function(mu,alpha, beta,gamma,time)
{
  z = mu+beta*exp(gamma-alpha*time)/(1+exp(gamma-alpha*time))
  #print(paste(z,time,sep="-"))
  if(z<=1e-5)return(1e-5)
  else return(z)
}

#lik=x=seq(0.01,9.99,length=100)
#for(i in 1:length(x))lik[i]=mu_t(mu,alpha,beta,gamma,x[i])
#plot(x,lik)
#integral_1 = function(alpha, gamma,age, time2)
#{
#  integral=function(x) #function under intergral in P(t,T)
#  {
#    -log(exp(gamma-alpha*x)+1)/alpha
#  }
  
#  return(integral(time2)-integral(time2-age))
#}

#den_age = function(mu, alpha, beta, gamma, lambda,time1, age,time2)
#{
#  integrand=function(x) 
#  {
#    exp(-(mu + lambda)*x)*exp(-beta*integral_1(alpha,gamma,  x, time2))
#  }
#  integral <-integrate(integrand,time1,time2)$value
#  #integral <-adaptIntegrate(integrand,time1,time2)$integral
#  nume = integrand(age) #define the numeritor of density
#  return(nume/integral)
#}


#phi_t = function(mu, alpha, beta, gamma, lambda,time1, time2) # mean loss rate from 0 to time2
#{
#  integrand=function(x) #function under intergral in P(t,T)
#  {
#    mu_t(mu,alpha, beta,gamma,x)*den_age(mu, alpha, beta, gamma, lambda,time1, x,time2)
#  }
#  integral <-integrate(integrand,time1,time2)$value
#  if(integral <= 0) 
#  {
#    error("the integral in P_tz is <= 0")
#    integral=1e-5
#  }
#  return(integral)
#}

phi_t = function(mu, alpha, beta, gamma, lambda,time1, time2) # mean loss rate from 0 to time2
{
  integrand_1=function(x) #function under intergral in P(t,T)
  {
    beta*exp(gamma-alpha*x)/(1+exp(gamma-alpha*x))*exp(-(lambda+mu)*x - beta/alpha*log(exp(gamma-alpha*(time2 - x)) + 1))
  }
  integral_1 <-integrate(integrand_1,time1,time2)$value
  
  integrand_2=function(x) #function under intergral in P(t,T)
  {
    exp(-(lambda+mu)*x - beta/alpha*log(exp(gamma-alpha*(time2 - x)) + 1))
  }
  integral_2 <-integrate(integrand_2,time1,time2)$value
  
  
  if(integral_1 <= 0 | integral_2<=0) 
  {
    error("the integral in P_tz is <= 0")
    integral=1e-5
  }
  return(mu+integral_1 / integral_2)
}


rou <- function(mu, alpha, beta, gamma, lambda,time1, time2) #time2 > time1
{
  integrand=function(x) #function under intergral in P(t,T)
  {
    phi_t(mu, alpha, beta, gamma, lambda,0, x)
  }
  integral <-integrate(Vectorize(integrand),time1,time2)$value
  
  return(integral-lambda*(time2-time1))
}

P_tz <- function(mu, alpha, beta, gamma, lambda,time1, time2)
{
  integrand=function(time) #function under intergral in P(t,T)
  {
    phi_t(mu, alpha, beta, gamma, lambda, 0, time)*exp(rou(mu, alpha, beta, gamma, lambda,time1, time2=time))
  }
  integral <-integrate(Vectorize(integrand),time1,time2)$value
  if(integral <= 0) 
  {
    error("the integral in P_tz is <= 0")
    integral=1e-5
  }
  result <- 1/(1+integral)
  return (result)
}


u_t <- function(mu, alpha, beta, gamma, lambda,time1, time2)
{
  if(time2 < time1) print("time2 must be greater than time1")
  result <- 1 - P_tz(mu, alpha, beta, gamma, lambda,time1=time1, time2=time2) * exp(rou(mu, alpha, beta, gamma, lambda, time1=time1, time2=time2))
  if(result > 1 | result < 0) print("u_t is wrong")
  return(result)
}


yita <- function(mu, alpha, beta, gamma, lambda, time1, time2,T)
{
  result <- 1- (1-u_t(mu, alpha, beta, gamma, lambda,time1=time1, time2=time2))/(1-(1- P_tz (mu, alpha, beta, gamma, lambda,time1=time2, time2=T))*u_t(mu, alpha, beta, gamma, lambda,time1=time1, time2=time2)) 
  if(result > 1 | result < 0) print("yita is wrong")
  return(result)
}


logcon_pdf_time_interval_single <- function(mu, alpha, beta, gamma, lambda, time1, x, current_time, start_ntaxa, end_ntaxa)  #for simulation
{
  
  if(end_ntaxa < start_ntaxa) print("end_ntaxa must be greater than start_ntaxa")
  
  loglike <- log(lambda)+log(P_tz(mu, alpha, beta, gamma, lambda, x,current_time))+log(end_ntaxa-start_ntaxa)+log(1-u_t(mu, alpha, beta, gamma, lambda, x,current_time))+(end_ntaxa-start_ntaxa-1)*log(u_t(mu, alpha, beta, gamma, lambda, x,current_time))-(end_ntaxa-start_ntaxa)*log(u_t(mu, alpha, beta, gamma, lambda, time1, current_time))
  return(loglike)
}

logcon_pdf_time_interval <- function(mu, alpha, beta, gamma, lambda,time1, x, current_time) 
{
  
  #check if x is between time1 and current_time
  x <- c(time1, x)  #adding time1 to the duplication time list
  ndup <- length(x)-1
  end_ntaxa <- ndup + 2
  x <- sort(x)
  if(x[1] < time1 | x[ndup] > current_time) print("the duplication times are not between time1 and current_time")
  
  loglike <- log(factorial(end_ntaxa-2))-(end_ntaxa-2)*log(u_t(mu, alpha, beta, gamma, lambda,time1, current_time))
  for(i in 1:ndup)
  {  
    loglike <- loglike + log(lambda)+log(1-u_t(mu, alpha, beta, gamma, lambda, x[i+1], current_time))+log(P_tz(mu, alpha, beta, gamma, lambda,x[i+1], current_time))
  }
  return(loglike)
}

#############################
## plot of likelihood

lik.sub.age.sharp=x=seq(0.01,9.99,length=100)
for(i in 1:length(x))
{
  lik.sub.age.sharp[i]=logcon_pdf_time_interval_single(mu,alpha, beta,gamma,lambda,time1, x[i], current_time, start_ntaxa, end_ntaxa)  #for simulation
}
plot(x,exp(lik.sub.age),type="l",lty=3,lwd=2)
lines(x,exp(lik.sub.age),lty=3,lwd=2)
legend("topright",c("Non","Neo","Sub"),lty=1:3)
sub.0.4=lik.sub.age
sub.0.3=lik.sub.age
sub.0.2=lik.sub.age
sub.0.1=lik.sub.age

#####################################
# simulation for sub
#####################################

sim_duplication_time <- function(mu,alpha, beta,gamma,lambda, time1, current_time, start_ntaxa, end_ntaxa)
{
  ndup <- end_ntaxa - start_ntaxa 
  duptime <- rep(1e-10, ndup+1)
  for(i in 1:ndup)
  {
    starttime <- duptime[i]
    ntaxa <- start_ntaxa + (i-1)
    fn=function(x)
    {
      return(-exp(logcon_pdf_time_interval_single(mu,alpha, beta,gamma, lambda, time1=starttime,x, current_time, start_ntaxa=ntaxa, end_ntaxa)))
    }
    optim=optim(1, fn, lower=starttime+1e-5, upper = current_time-1e-10, method="L-BFGS-B")
    upper_bound =(-optim$value)
    
    stop = 0
    while(stop==0)
    {
      t <- runif(1, starttime, current_time)
      if(runif(1)*upper_bound < exp(logcon_pdf_time_interval_single(mu,alpha, beta,gamma,lambda, time1=starttime, t, current_time, start_ntaxa=ntaxa, end_ntaxa=end_ntaxa)))
      {
        time<- t
        stop <- 1
      }
    }
    duptime[i+1] <- time
  }
  duptime[-1]
}

### simulating data
set.seed(1)
dup=c(10,20,30,40,60,80,100)
#dup=c(10)
###for sub  
#mu=0.8
for(i in 1:length(dup))
{
  N=dup[i]
  m=1
  mu=-0.1
  alpha=0.6
  beta<-0.9
  gamma<-5
  lambda<-0.2
  time1=1e-10
  current_time <- 10
  start_ntaxa=2
  end_ntaxa=start_ntaxa+N   #define number of duplicates
  M=30 # the duplication
  
  #################### generate tree
  #set.seed(1)
  tree=array(dim = c(M,end_ntaxa-start_ntaxa))
  for(i in 1:M)
  {
    
    tree[i,] <- sim_duplication_time (mu,alpha, beta,gamma,lambda, time1, current_time, start_ntaxa, end_ntaxa)
  }
  write.table(tree,paste("tree-sub-",N,".txt",sep=""))
  
}

#ave.sub=apply(tree,2,mean)


################################################################################
# 
# @brief General Markov chain Monte Carlo algorithm using Metropolis-Hastings.
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-11-06, version 1.1
#
# @param    likelihoodFunction     function      the log-likelihood function
# @param    priors                 list          a list of functions of the log-prior density per parameter
# @param    parameters             vector        initial set of paramete values
# @param    logTransform           vector        should be a log-transform be used for parameter proposals
# @param    iterations             scalar        the number of iterations
# @param    burnin                 scalar        number of burnin iterations
# @param    thining                scalar        number of iterations between samples
# @param    verbose                boolean       should we print information during the MCMC?
# @return                          list          samples from the posterior
#
################################################################################
require(coda)

tess.mcmc <- function(likelihoodFunction,priors,parameters,current_time,logTransforms,delta,iterations,burnin=round(iterations/3),thining=1,verbose=FALSE) {
  
  OPTIMIZATIONS <- min(ceiling(burnin/20),10)
  PRINT_FREQ_B <- min(burnin,20)
  PRINT_FREQ_C <- min(iterations,20)
  
  # create a list for the samples
  chain = array(dim = c(floor(iterations/thining)+1,length(parameters))) #reserve memory for the chain, for large chains we might consider writing to a file instead of storing in memory
  
  
  # pre-compute current posterior probability
  pp <- likelihoodFunction(parameters)
  for ( j in 1:length(parameters) ) {
    pp <- pp + priors[[j]](parameters[j])
  }
  
  if ( verbose == TRUE ) {
    cat("Burning-in the chain ...\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*")
  }
  
  # this is the tuning parameter:  delta <- rep(1,length(parameters))
  accepted <- rep(0,length(parameters))
  tried <- rep(0,length(parameters))
  
  for (i in 1:(burnin+iterations)) {
    
    if ( verbose == TRUE ) {
      if ( i <= burnin ) {
        if ( i %% (burnin/PRINT_FREQ_B) == 0 ) {
          cat("**")
        }
      } else if (i == (burnin+1) ) {
        cat("\nFinished burnin period!\n\n")
        cat("Running the chain ...\n")
        cat("0--------25--------50--------75--------100\n")
        cat("*")
      } else {
        if ( (i-burnin) %% (iterations/PRINT_FREQ_C) == 0 ) {
          cat("**")
        }
      }
    }
    
    
    # propose new values
    #for ( j in 1:length(parameters) ) {
    for ( j in c(2,5) ) {  # only update lambda and beta
      # increase our tried counter
      tried[j] <- tried[j] + 1
      
      if ( logTransforms[j] == TRUE ) {
        if (parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }
        
        eta           <- log(parameters[j]) ### propose a new value for parameter[j]
        
        if(j==2)
        {
          new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
          new_val       <- exp(new_eta) 
          # calculate the Hastings ratio
          hr            <- log(new_val / parameters[j])
          parameters[j] <- new_val
          while(mu_t(parameters[j-1],parameters[j],parameters[j+1],parameters[j+2],current_time)==parameters[j+3])
          {
            new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
            new_val       <- exp(new_eta)
            hr            <- log(new_val / parameters[j])
            parameters[j] <- new_val
          }
        }
        else if(j==3)
        {
          new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
          new_val       <- exp(new_eta) 
          # calculate the Hastings ratio
          hr            <- log(new_val / parameters[j])
          parameters[j] <- new_val
          while(mu_t(parameters[j-2],parameters[j-1],parameters[j],parameters[j+1],current_time)==parameters[j+2])
          {
            new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
            new_val       <- exp(new_eta)
            hr            <- log(new_val / parameters[j])
            parameters[j] <- new_val
          }
        }
        else if(j==4)
        {
          new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
          new_val       <- exp(new_eta) 
          # calculate the Hastings ratio
          hr            <- log(new_val / parameters[j])
          parameters[j] <- new_val
          while(mu_t(parameters[j-3],parameters[j-2],parameters[j-1],parameters[j],current_time)==parameters[j+1])
          {
            new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
            new_val       <- exp(new_eta)
            hr            <- log(new_val / parameters[j])
            parameters[j] <- new_val
          }
        }
        else if(j==5)
        {
          new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
          new_val       <- exp(new_eta) 
          # calculate the Hastings ratio
          hr            <- log(new_val / parameters[j])
          parameters[j] <- new_val
          while(mu_t(parameters[j-4],parameters[j-3],parameters[j-2],parameters[j-1],current_time)==parameters[j])
          {
            new_eta = runif(1,max(log(1e-5),eta-delta[j]),eta+delta[j])
            new_val       <- exp(new_eta)
            hr            <- log(new_val / parameters[j])
            parameters[j] <- new_val
          }
        }
     
        new_val       <- exp(new_eta)
        hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- 0.0
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        if ( is.finite(new_pp) ) {
          new_pp        <- new_pp + likelihoodFunction(parameters)
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
          accepted[j] <- accepted[j] + 1
        } else {
          parameters[j] <- exp(eta)
        }
      } 
      else {
        
        eta           <- parameters[j] ### propose a new value for parameter[j]
        new_eta=runif(1,eta-delta[j],eta+delta[j])
                
        new_val       <- new_eta
        hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- 0.0
        for ( k in 1:length(parameters) ) {
          new_pp <- new_pp + priors[[k]](parameters[k])
        }
        if ( is.finite(new_pp) ) {
          new_pp        <- new_pp + likelihoodFunction(parameters)
        }
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
          accepted[j] <- accepted[j] + 1
        } else {
          parameters[j] <- eta
        }
        
      }
      
    }
    
    # sample the parameter
    if (i >= burnin) {
      if ( (i-burnin) %% thining == 0 ) {
        chain[(i-burnin)/thining+1,] <- parameters
      }
    }
    
  }
  
  if ( verbose == TRUE ) {
    cat("\nFinished MCMC!\n\n")
    cat("Parameter\t| delta\t\t| Acceptance Probability\n")
    cat("===============================================================\n")
    for ( j in 1:length(parameters) ) {
      cat(sprintf("%i\t\t| %.3f\t\t| %.3f\n",j,delta[j],accepted[j]/tried[j]))
    }
  }
  
  
  return(as.mcmc(chain)) #return a mcmc object, used by coda to plot
}

###############construct parameters in tess.mcmc
tree=read.table(paste("tree-sub-",N,".txt",sep=""),header=T)
loglike <- function(para)
{
  loglike=logcon_pdf_time_interval(para[1], para[2], para[3], para[4],para[5], time1=0, as.matrix(tree[m,]), current_time)
  return(loglike)    
}

muf=function(mu)dunif(mu, min=-0.5, max=0.5, log = T)
alphaf=function(alpha)dunif(alpha, min=0, max=1, log = T)
betaf=function(beta)dunif(beta, min=0.5, max=1.5, log = T)
gammaf=function(gamma)dunif(gamma, min=0, max=10, log = T) 
lambdaf=function(lambda)dunif(lambda, min=0, max=0.5, log = T)

#priors=list(muf,alphaf,betaf,gammaf,lambdaf)
priors=list(muf,alphaf,betaf,gammaf,lambdaf)
parameters=c(-0.1,0.6,0.9,5,0.2)
logTransforms=c(FALSE,rep(TRUE,4)) #whether do transform when propose new value
delta=c(0.02,0.1,0.1,0.1,0.2)
iterations=10
ptm <- proc.time()
para.chain=tess.mcmc(loglike,priors,parameters,current_time,logTransforms,delta,iterations)
proc.time() - ptm
summary(para.chain)
write(t(as.matrix(summary(para.chain)[1]$"statistics"[,1])),paste("est-sub",N,m,sep="-"))

#write.table(para.chain,file=paste("sub_chain_",N,".txt",sep=""),sep="\t")

#pdf(paste("sub_chain_",N,".pdf",sep=""))
#plot(para.chain)
#dev.off()

##################################  AIC
##############################################

dup=c(10,20,30,40,60,80,100)
sub.tnon.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-non-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],para[3],para[4],para[5],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=5
    para=c(-0.1,0.6,9,0.5,0.2)
    sub.tnon.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write(sub.tnon.aic,"sub.tnon.aic")


dup=c(10,20,30,40,60,80,100)
sub.tneo.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-neo-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],para[3],para[4],para[5],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=5
    para=c(-0.1,0.6,9,0.5,0.2)
    sub.tneo.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write(sub.tneo.aic,"sub.tneo.aic")

dup=c(10,20,30,40,60,80,100)
sub.tsub.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-sub-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],para[3],para[4],para[5],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=5
    para=c(-0.1,0.6,9,0.5,0.2)
    sub.tsub.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write(sub.tsub.aic,"sub.tsub.aic")


