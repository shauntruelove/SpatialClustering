#######################################################
# Code for running stochastic compartmental models
#######################################################



#####################STOCHASTIC FIXED TIMESTEP SIR MODELS#################
#Runs a simple stochastic SIR model.
#Parameters
#  beta - the force of infection
#  gamma - the recovery rate
#  initial.state - the initial state of the system.
#                MUST BE WHOLE NUMBERS!
#  time.step -  what time step should we use
#  freq.dependent - should we use a frequency or densitiy
#             dependent model
#  final.only - speeds things up a little if we are only interested in
#           the final state of the epidemic
#
#Returns
#   a matrix with 4 columns: t, S, I, and R
runStochasticSIR <- function(beta=.001,
                             gamma=.5,
                             initial.state= c(S=999, I=1, R=0),
                             step.size = 1/4,
                             freq.dependent=FALSE,
                             final.only=FALSE) {

  #If the model is frequency dependent we modify beta
  #based on the total populations size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+initial.state["I"]+initial.state["R"], 1)

  #create the parameter vector.
  param <- c(beta=beta/beta.divisor, gamma=gamma)

  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the epidemic ends once there are no mode susceptibles.
  t <- 0
  y <- initial.state


  sir.output <- matrix(ncol=4, nrow=1)
  colnames(sir.output) <- c("time", "S","I","R")
  sir.output[1,] <- c(t,y)

  while (y["I"]>0) {

    t <- t+step.size

    delta <- stochastic.dx.dt(step.size, y, param)

    y <- y+delta

    #trick to speed up the code
    if (!final.only) {
      sir.output<-rbind(sir.output, c(t,y))
    }
  }

  if(final.only) {
    sir.output[1,]<-c(t,y)
  }

  return(sir.output)
}


#The dx.dt method for a stochastic SIR using simple
#rates (e.g., a constant probability of being moved
#from one comparment to the other)
#
#Parameters:
#  step.size -  the size of the step being taken compared
#         to that used to specify the parameters
#
#  y - the current state of the system.
#  param - the parameters of the system
#
#Returns:
#  a list representing the changes in each compartment of the system

stochastic.dx.dt <- function(step.size, y, param) {
  
  #calculate the probability of infection and recovery in this time step
  p.infect  <- 1-exp(-step.size*param["beta"]*y["I"])
  p.recover <- 1-exp(-step.size*param["gamma"])
  
  #Do our random stuff so we know how much to change:
  #The number of incident cases follows a binomial distrib where N=S and p=p.infect
  incident.cases  <- rbinom(1, y["S"], p.infect)

  #the number of recovered cases follows a binomial distrib where N=I and p=p.recover
  recovered.cases <- rbinom(1, y["I"], p.recover)

  #Find the deltas for our compartments
  dS <- -incident.cases
  dI <- incident.cases - recovered.cases
  dR <- recovered.cases

  #Those susceptibles move to the recovered compartment
  return(c(dS,dI,dR))
}




##############STOCHASTIC EVENT DRIVEN SIR MODELS#####################
#Runs an event driven version of the stochastic SIR model.
#Note that we no longer need to specify a time step
#in this model
#
#Parameters
#  beta - the force of infection
#  gamma - the recovery rate
#  initial.state - the initial state of the system.
#                MUST BE WHOLE NUMBERS!
#  freq.dependent - should we use a frequency or densitiy
#             dependent model
#  final.only - speeds things up a little if we are only interested in
#           the final state of the epidemic
#
#Returns
#   a matrix with 4 columns: t, S, I, and R

runEventDrivenSIR <- function(beta=.001,
                              gamma=.5,
                              initial.state= c(S=999, I=1, R=0),
                              freq.dependent=FALSE,
                              final.only=FALSE) {

  #If the model is frequency dependent we modify beta
  #based on the total populations size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+initial.state["I"]+initial.state["R"],
                         1)

  #create the parameter vector.
  param <- c(beta=beta/beta.divisor, gamma=gamma)

  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the epidemic ends once there are no more susceptibles.
  t <- 0
  y <- initial.state

  sir.output <- matrix(ncol=4, nrow=1)
  colnames(sir.output) <- c("time", "S","I","R")
  sir.output[1,] <- c(t,y)

  while (y["I"]>0) {
    d.nxt <- nxt.evt.SIR(t,y,param)
    t <- t+d.nxt[1]
    y <- y+d.nxt[2:4]

    #trick to speed up the code
    if (!final.only) {
      sir.output<-rbind(sir.output, c(t,y))
    }
  }

  if(final.only) {
    sir.output[1,]<-c(t,y)
  }

  return(sir.output)
}



#gets the next event in the event driven model
#Parameters:
#  t - the current system time
#  y - the current state of the system.
#  param - the parameters of the system
#
#Return:
#  a vector of the form {dt,dS,dI,dR)
nxt.evt.SIR <- function (t, y, param) {
  #Get the current event rates
  #note that multiple exponential processes give event rates
  #at the sum of their rates
  infection.rate <- y["S"]*y["I"]*param["beta"]
  recovery.rate <- y["I"]*param["gamma"]

  #DEBUG
  ## if (is.na(infection.rate)| is.na(recovery.rate)) {
##     print("BAD INF/REC")
##     print(y)
##   }

  #now get the event specific next times
  #we need to check if the infection rate
  #is 0 (i.e., there are no more susceptibles)
  if (infection.rate>0) {
    time.to.i <- rexp(1,rate=infection.rate)
    time.to.r <- rexp(1,rate=recovery.rate)
  } else {
    time.to.i <- Inf
    time.to.r <- rexp(1,rate=recovery.rate)
  }

  #DEBUG
 ##  if (is.na(time.to.i) | is.na(time.to.r)) {
##     print("BAD TIMES")
##     print(y)
##     cat(infection.rate ,":", time.to.i,"\n") #DEBUG
##     cat(recovery.rate,"-", time.to.r,"\n") #DEBUG
##   }

  if (time.to.i<time.to.r) {
    rc <- c(time.to.i, -1, 1,0)
  } else {
    rc <- c(time.to.r,0,-1,1)
  }

  return(rc)
}





#  ... -  parameters to pass to runStochasticSIR
#Returns:
#   A vector containing the final sizes of the epidemics
simulateFinalSizes <- function(n,evt.driven=FALSE, ...) {
  final.sizes <- vector(length = n)
  for (i in 1:n) {
    if (evt.driven) {
      epi <- runEventDrivenSIR(final.only=TRUE,...)
    } else {
      epi <- runStochasticSIR(final.only=TRUE,...)
    }

    final.sizes[i] <- epi[nrow(epi),"R"]
  }
  return(final.sizes)
}

