rm(list=ls())
N=30 #number of duplication times
m=1 # order of simulated tree among 100 trees 
mu=0.8
lambda<-0.2
time1=1e-10
current_time <- 10
start_ntaxa=2
end_ntaxa=start_ntaxa+N  
M=30 # number of samples

#########################################################
#   nonfunctionalization
#########################################################

rou <- function(mu, lambda,time1, time2) #time2 > time1
{
  if(time2 < time1) print("time2 must be greater than time1")
  y<- (mu-lambda)*(time2-time1)
  return (y)
}

P_tz <- function(mu, lambda,time1, time2) #lambda!=mu
{
  if(mu==lambda)integral=mu*(time2-time1)
  else {integral =mu/(mu-lambda)*(exp((mu-lambda)*(time2-time1))-1)}
  if(integral <= 0) 
  {
    
    print(paste(mu,lambda,integral,sep="-"))
  }
  result <- 1/(1+integral)
  return (result)
}

u_t <- function(mu, lambda, time1, time2)
{
  if(time2 < time1) print("time2 must be greater than time1")
  result <- 1 - P_tz(mu,lambda,time1=time1, time2=time2) * exp(rou(mu, lambda, time1=time1, time2=time2))
  if(result > 1 | result < 0) {print("u_t is wrong");print(lambda)}
  return(result)
  
}


yita = function(mu,lambda,time1,time2,T)
{
  result = 1 - (1-u_t(mu,lambda,time1,T))/(1-u_t(mu,lambda,time2,T))
  if(result > 1 | result < 0) print("yita")
  return(result)
}


logcon_pdf_time_interval_single <- function(mu, lambda, time1, x, current_time, start_ntaxa, end_ntaxa)  #for simulation
{
  
  if(end_ntaxa < start_ntaxa) print("end_ntaxa must be greater than start_ntaxa")
  
  loglike <- log(lambda)+log(P_tz(mu, lambda,x,current_time))+log(end_ntaxa-start_ntaxa)+log(1-u_t(mu, lambda, x,current_time))+(end_ntaxa-start_ntaxa-1)*log(u_t(mu, lambda,x,current_time))-(end_ntaxa-start_ntaxa)*log(u_t(mu, lambda,time1, current_time))
  return(loglike)
}


logcon_pdf_time_interval <- function(mu, lambda, time1, x, current_time) 
{
  
  #check if x is between time1 and current_time
  x <- c(time1, x)  #adding time1 to the duplication time list
  ndup <- length(x)-1
  end_ntaxa <- ndup + 2
  x <- sort(x)
  if(x[1] < time1 | x[ndup] > current_time) print("the duplication times are not between time1 and current_time")
  
 loglike <- log(factorial(end_ntaxa-2))-(end_ntaxa-2)*log(u_t(mu,lambda,time1, current_time))
  for(i in 1:ndup)
  {  
    loglike <- loglike + log(lambda)+log(1-u_t(mu,lambda,x[i+1], current_time))+log(P_tz(mu,lambda,x[i+1], current_time))
  }
  return(loglike)

}

###########################################################
### plot 1st density curve
lik.non.age = x=seq(0.01,9.99,length=100)
for(i in 1:length(x))
{
  lik.non.age[i] = logcon_pdf_time_interval_single(mu, lambda, time1, x[i], current_time, start_ntaxa, end_ntaxa)  #for simulation
}
#par(mfrow=c(1,2))
plot(x,exp(lik.non.age),ylab="Probability Density",xlab="Time", type="l",lty=1,lwd=2,main="(a)")
lines(x,exp(lik.non.age),ylab="Probability Density",xlab="Time", type="l",lty=1,lwd=2,main="(a)")


lik = x=seq(0.01,9.99,length=100)
for(i in 1:length(x))
{
  lik[i] = logcon_pdf_time_interval_single(mu, lambda, time1, x[i], current_time, start_ntaxa, end_ntaxa)  #for simulation
}

#################### generate tree##############
# simulation for non

sim_duplication_time <- function(mu, lambda, time1, current_time, start_ntaxa, end_ntaxa)
{
  ndup <- end_ntaxa - start_ntaxa 
  duptime <- rep(1e-10, ndup+1)
  for(i in 1:ndup)
  {
    starttime <- duptime[i]
    ntaxa <- start_ntaxa + (i-1)
    fn=function(x)
    {
      return(-exp(logcon_pdf_time_interval_single(mu,lambda, time1=starttime,x, current_time, start_ntaxa=ntaxa, end_ntaxa)))
    }
    optim=optim(1, fn, lower=starttime+1e-5, upper = current_time-1e-10, method="L-BFGS-B")
    upper_bound =(-optim$value)
    
    stop = 0
    while(stop==0)
    {
      t <- runif(1, starttime, current_time)
      if(runif(1)*upper_bound < exp(logcon_pdf_time_interval_single(mu, lambda, time1=starttime, t, current_time, start_ntaxa=ntaxa, end_ntaxa=end_ntaxa)))
      {
        time<- t
        stop <- 1
        print(paste(i,time))
      }
    }
    duptime[i+1] <- time
  }
  duptime[-1]
}

set.seed(1)
dup=c(10,20,30,40,60,80,100)
for(i in 1:length(dup))
{
  N=dup[i] #number of duplication times
  m=1 # order of simulated tree among 100 trees 
  ###for neo                                                                                                                                                                                                                                                                                                                                                          
  #mu=0.8
  mu=0.8
  lambda<-0.2
  time1=1e-10
  current_time <- 10
  start_ntaxa=2
  end_ntaxa=start_ntaxa+N  #define number of duplicates
  M=30 # the duplica
  
  #################### generate tree
  
  tree=array(dim = c(M,end_ntaxa-start_ntaxa))
  for(i in 1:M)
  {
    tree[i,] <- sim_duplication_time (mu, lambda, time1, current_time, start_ntaxa, end_ntaxa)
  }
  write.table(tree,paste("tree-non-",N,".txt",sep=""))
}


################################################################################
# 
# Random walk metropolis
#
################################################################################

require(coda)

tess.mcmc <- function(likelihoodFunction,priors,parameters,logTransforms,delta,iterations,burnin=round(iterations/3),thining=1,verbose=FALSE) 
{
  # create a list for the samples
  chain = array(dim = c(floor(iterations/thining)+1,length(parameters))) #reserve memory for the chain, for large chains we might consider writing to a file instead of storing in memory
  # pre-compute current posterior probability
  pp <- likelihoodFunction(parameters)
  
  for ( j in 1:length(parameters) ) {
    pp <- pp + priors[[j]](parameters[j])
  }

  # this is the tuning parameter:  delta <- rep(1,length(parameters))
  accepted <- rep(0,length(parameters))
  tried <- rep(0,length(parameters))
  
  for (i in 1:(burnin+iterations)) 
  {
    # propose new values
    for ( j in 1:length(parameters) ) 
    {
      # increase our tried counter
      tried[j] <- tried[j] + 1
      
      if ( logTransforms[j] == TRUE ) 
      {
        if (parameters[j] == 0) 
        {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }
        eta           <- log(parameters[j]) ### propose a new value for parameter[j]
        if(j==1)
          new_eta = runif(1,max(log(0.6),eta-delta[j]),min(eta+delta[j],log(1)))
        else if(j==2)
          new_eta = runif(1,eta-delta[j],min(eta+delta[j],log(0.5)))
        
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
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) 
        {
          pp <- new_pp
          accepted[j] <- accepted[j] + 1
        } else {
          parameters[j] <- exp(eta)
        }
      } 
      else 
      {
        
      }
      
    }
    
    # sample the parameter
    if (i >= burnin) 
    {
      if ( (i-burnin) %% thining == 0 ) 
      {
        chain[(i-burnin)/thining+1,] <- parameters
      }
    }
    
  }
  
  return(as.mcmc(chain)) #return a mcmc object, used by coda to plot
}


###############construct parameters in tess.mcmc
tree=read.table(paste("tree-non-",N,".txt",sep=""),header=T)
#ave.non=apply(tree,2,mean)
loglike <- function(para)
{
  loglike=logcon_pdf_time_interval(para[1], para[2],time1, as.matrix(tree[m,]), current_time)
  return(loglike)    
}

muf=function(mu)dunif(mu, min=0.6, max=1, log = T)
lambdaf=function(lambda)dunif(lambda, min=0, max=0.5, log = T)
priors=list(muf,lambdaf)

parameters=c(0.8,0.2)
logTransforms=rep(TRUE,length(parameters)) #whether do transform when propose new value
delta=rep(0.1,length(parameters))
iterations=50000

para.chain=tess.mcmc(loglike,priors,parameters,logTransforms,delta,iterations)
summary(para.chain)
write.table(para.chain,file=paste("non_chain_",N,".txt",sep=""),sep="\t")

pdf(paste("non_chain_",N,".pdf",sep=""))
plot(para.chain)
dev.off()

##################################  AIC
##############################################
dup=c(10,20,30,40,60,80,100)
non.tnon.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-non-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=2
    para=c(0.8,0.2)
    non.tnon.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write.table(non.tnon.aic,file="non.tnon.aic.txt")

dup=c(10,20,30,40,60,80,100)
non.tneo.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-neo-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=2
    para=c(0.8,0.2)
    non.tneo.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write.table(non.tneo.aic,file="non.tneo.aic.txt")


dup=c(10,20,30,40,60,80,100)
non.tsub.aic = matrix(nrow=30,ncol=length(dup))
for(j in 1:length(dup))
{
  N=dup[j]
  tree=read.table(paste("tree-sub-",N,".txt",sep=""),header=T)
  
  for(i in 1:30)
  {
    loglike <- function(para)
    {
      loglike=logcon_pdf_time_interval(para[1], para[2],time1, as.matrix(tree[i,]), current_time)
      return(loglike)    
    }
    
    K=2
    para=c(0.8,0.2)
    non.tsub.aic[i,j] = -2 *( loglike(para) ) + 2*K
    
  }
  
}
write.table(non.tsub.aic,file="non.tsub.aic.txt")

