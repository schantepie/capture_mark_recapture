a=commandArgs()[4]

library(R2jags)
load("Mod_sans_trans.RData")

##########################"""age*sex*site*year
if(a==1){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
    for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
        for(sexi in 1:2)  {
           for(sitei in 1:2)  {
              for(anneei in 1:2)  {
                       gamma[agei,sexi,sitei,anneei]~ dnorm(0, 0.0001)
                       S.s2[agei,sexi,sitei,anneei]<-1/(1+exp(-gamma[agei,sexi,sitei,anneei]))
              }
            }
        }
    }

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],sex[i],site[i],annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2,2,2)), r= runif(1, 0, 1),z = init.state,
                         lpsiA = array(rnorm(64),dim=c(2,2,3)), 
                         lpsiB =array(rnorm(64),dim=c(2,2,3)), 
                         lpsiC = array(rnorm(64),dim=c(2,2,3)), 
                         lpsiD = array(rnorm(64),dim=c(2,2,3)))}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}



###########################  age*sex*site
if(a==2){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
      for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
        for(sexi in 1:2)  {
           for(sitei in 1:2)  {
             
                       gamma[agei,sexi,sitei]~ dnorm(0, 0.0001)
                       S.s2[agei,sexi,sitei]<-1/(1+exp(-gamma[agei,sexi,sitei]))
              }
            }
        }
    

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],sex[i],site[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


################################################ age*site
if(a==3){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
           for(sitei in 1:2)  {
                       gamma[agei,sitei]~ dnorm(0, 0.0001)
                       S.s2[agei,sitei]<-1/(1+exp(-gamma[agei,sitei]))
              }
            }

    

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],site[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


################################  age

if(a==4){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
      for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {

       
             
                       gamma[agei]~ dnorm(0, 0.0001)
                       S.s2[agei]<-1/(1+exp(-gamma[agei]))
              }
         

    

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}



######################################  cst

if(a==5){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


                       gamma~ dnorm(0, 0.0001)
                       S.s2<-1/(1+exp(-gamma))

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


######################### age*sex

if(a==6){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
        for(sexi in 1:2)  {
          
                       gamma[agei,sexi]~ dnorm(0, 0.0001)
                       S.s2[agei,sexi]<-1/(1+exp(-gamma[agei,sexi]))
              }
            }
  

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],sex[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


########################### sex*site

if(a==7){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
       for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(sexi in 1:2)  {
              for(sitei in 1:2)  {
                       gamma[sexi,sitei]~ dnorm(0, 0.0001)
                       S.s2[sexi,sitei]<-1/(1+exp(-gamma[sexi,sitei]))
              }
 }


    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[sex[i],site[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

#######################################" sex

if(a==8){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
      for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)



        for(sexi in 1:2)  {
 
                       gamma[sexi]~ dnorm(0, 0.0001)
                       S.s2[sexi]<-1/(1+exp(-gamma[sexi]))
              }
      
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[sex[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


##########################################" site
if(a==9){
sink("ms3-multinomlogit.jags")
cat("
    model {
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }

    r ~ dunif(0, 1)
    p ~ dunif(0, 1)



           for(sitei in 1:2)  {
                       gamma[sitei]~ dnorm(0, 0.0001)
                       S.s2[sitei]<-1/(1+exp(-gamma[sitei]))
              }
 
   
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){

    logit(SS[i,t])<-gamma[site[i]]

    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

####################################"  annee
if(a==10){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
    for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


              for(anneei in 1:2)  {
                       gamma[anneei]~ dnorm(0, 0.0001)
                       S.s2[anneei]<-1/(1+exp(-gamma[anneei]))
              }
 
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


##################################### age * annee
if(a==11){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
      for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
              for(anneei in 1:2)  {
                       gamma[agei,anneei]~ dnorm(0, 0.0001)
                       S.s2[agei,anneei]<-1/(1+exp(-gamma[agei,anneei]))
              }
            }
  
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

################################"  age * site * annee
if(a==12){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


  for(agei in 1:2)  {
           for(sitei in 1:2)  {
              for(anneei in 1:2)  {
                       gamma[agei,sitei,anneei]~ dnorm(0, 0.0001)
                       S.s2[agei,sitei,anneei]<-1/(1+exp(-gamma[agei,sitei,anneei]))
              }
            }
        }


    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[xage[i,t],site[i],annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


###################################### site * year

if(a==13){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)


 
           for(sitei in 1:2)  {
              for(anneei in 1:2)  {
                       gamma[sitei,anneei]~ dnorm(0, 0.0001)
                       S.s2[sitei,anneei]<-1/(1+exp(-gamma[sitei,anneei]))
              }
            }
  

    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[site[i],annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}




###########################""" sex * annee

if(a==14){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
    for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)



        for(sexi in 1:2)  {
              for(anneei in 1:2)  {
                       gamma[sexi,anneei]~ dnorm(0, 0.0001)
                       S.s2[sexi,anneei]<-1/(1+exp(-gamma[sexi,anneei]))
              }
            }
   
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    

    logit(SS[i,t])<-gamma[sex[i],annee[i]]

   
  
    ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}



#############################  (age*sex*year)+site

if(a==15){
sink("ms3-multinomlogit.jags")
cat("
    model {
    for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    
    Sitei[1] <- 0
    Sitei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    
    for(agei in 1:2)  {
    for(sexi in 1:2)  {
    for(anneei in 1:2)  {
    gamma[agei,sexi,anneei]~ dnorm(0, 0.0001)
    S.site1[agei,sexi,anneei]<-1/(1+exp(-gamma[agei,sexi,anneei]))
    S.site2[agei,sexi,anneei]<-1/(1+exp(-gamma[agei,sexi,anneei]-Sitei[2]))
    }
    }
    }
    
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t],sex[i],annee[i]]+Sitei[site[i]]
    
    
   ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1

    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2,2)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.site1","S.site2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


###############################  (age*annee)+site+sex

if(a==16){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    
    
    Sitei[1] ~ dnorm(0, 0.1)I(-10,10)
    Sitei[2] <- 0 
    
    Sexi[1]~ dnorm(0, 0.1)I(-10,10)
    Sexi[2]<- 0
    
    for(agei in 1:2)  {
    for(anneei in 1:2)  {
    gamma[agei,anneei]~ dnorm(0, 0.0001)
    S.s2[agei,anneei]<-1/(1+exp(-gamma[agei,anneei]))
    } 
    }
    
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t],annee[i]]+Sitei[site[i]]+Sexi[sex[i]]
    
    
    
      ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(64),dim=c(2,2)),Sex=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.s2", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

############################ age+site+year+sex


if(a==17){
sink("ms3-multinomlogit.jags")
cat("
    model {
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    
    
    Sitei[1] <- 0 
    Sitei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    Sexi[1] <- 0 
    Sexi[2]~ dnorm(0, 0.1)I(-10,10)
    
    Anneei[1] <- 0 
    Anneei[2]~ dnorm(0, 0.1)I(-10,10)
    
    for(agei in 1:2)  {
    gamma[agei]~ dnorm(0, 0.0001)
    } 
    
    #S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t]]+Sitei[site[i]]+Sexi[sex[i]]+Anneei[annee[i]]
    
    
    
       ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=c(rnorm(1),rnorm(1)),Annee=c(NA,rnorm(1)),Sex=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c("S.Age1.Male.Annee1", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


################################# age+sex+site
if(a==18){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    
    
    Sitei[1] <- 0 
    Sitei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    Sexi[1] <- 0 
    Sexi[2]~ dnorm(0, 0.1)I(-10,10)
    
    
    for(agei in 1:2)  {
    gamma[agei]~ dnorm(0, 0.0001)
    } 
    
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t]]+Sitei[site[i]]+Sexi[sex[i]]
    
    
    
      ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=c(rnorm(1),rnorm(1)),Sex=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

#############################"annee+site+age


if(a==19){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    
    
    Anneei[1] <- 0 
    Anneei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    Sitei[1] <- 0 
    Sitei[2]~ dnorm(0, 0.1)I(-10,10)
    
    
    for(agei in 1:2)  {
    gamma[agei]~ dnorm(0, 0.0001)
    } 
    
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t]]+Sitei[site[i]]+Anneei[annee[i]]
    
    
    
       ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=c(rnorm(1),rnorm(1)),Annee=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

######################################### site+age

if(a==20){
sink("ms3-multinomlogit.jags")
cat("
    model {
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    Sitei[1] <- 0 
    Sitei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    for(agei in 1:2)  {
    gamma[agei]~ dnorm(0, 0.0001)
    } 
    
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t]]+Sitei[site[i]]
    
    
    
      ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=c(rnorm(1),rnorm(1)),Annee=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}

#################################### (age*area)+year+sex
if(a==21){
sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    Sexi[1] <- 0 
    Sexi[2] ~ dnorm(0, 0.1)I(-10,10)
    Anneei[1] <- 0 
    Anneei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    for(agei in 1:2)  {
    for (sitei in 1:2){
    gamma[agei,sitei]~ dnorm(0, 0.0001)
    } 
    }
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t],site[i]]+Sexi[sex[i]]+Anneei[annee[i]]
    
    
    
       ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(4),c(2,2)),Annee=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


################################  (age*area)+sex

if(a==22){

sink("ms3-multinomlogit.jags")
cat("
    model {
    
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
    Sexi[1] <- 0 
    Sexi[2] ~ dnorm(0, 0.1)I(-10,10)
    #     Annee[1] <- 0 
    #     Annee[2] ~ dnorm(0, 0.1)I(-10,10)
    
    for(agei in 1:2)  {
    for (sitei in 1:2){
    gamma[agei,sitei]~ dnorm(0, 0.0001)
    } 
    }
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t],site[i]]+Sexi[sex[i]]#+Anneei[annee[i]]
    
    
    
       ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(4),c(2,2)),Annee=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}


################################### (age*site)+year
if(a==23){
sink("ms3-multinomlogit.jags")
cat("
    model {
     for(age_p in 1:2)  {  
    for(site_p in 1:2)  {
    for (i_p in 1:3){
    lpsiA[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiB[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiC[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    lpsiD[site_p,age_p,i_p] ~ dnorm(0, 0.001)
    }
    
    # Constrain the transitions such that their sum is < 1
    
    psiA[site_p,age_p,1] <- 1 / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,2] <- exp(lpsiA[site_p,age_p,1]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,3] <- exp(lpsiA[site_p,age_p,2]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))
    psiA[site_p,age_p,4] <- exp(lpsiA[site_p,age_p,3]) / (1 + exp(lpsiA[site_p,age_p,1]) + exp(lpsiA[site_p,age_p,2])+ exp(lpsiA[site_p,age_p,3]))  
    
    psiB[site_p,age_p,1] <- exp(lpsiB[site_p,age_p,1]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,2] <- 1 / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,3] <- exp(lpsiB[site_p,age_p,2]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    psiB[site_p,age_p,4] <- exp(lpsiB[site_p,age_p,3]) / (1 + exp(lpsiB[site_p,age_p,1]) + exp(lpsiB[site_p,age_p,2])+ exp(lpsiB[site_p,age_p,3]))
    
    psiC[site_p,age_p,1] <- exp(lpsiC[site_p,age_p,1]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,2] <- exp(lpsiC[site_p,age_p,2]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,3] <- 1 / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    psiC[site_p,age_p,4] <- exp(lpsiC[site_p,age_p,3]) / (1 + exp(lpsiC[site_p,age_p,1]) + exp(lpsiC[site_p,age_p,2])+ exp(lpsiC[site_p,age_p,3]))
    
    psiD[site_p,age_p,1] <- exp(lpsiD[site_p,age_p,1]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,2] <- exp(lpsiD[site_p,age_p,2]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,3] <- exp(lpsiD[site_p,age_p,3]) / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    psiD[site_p,age_p,4] <- 1 / (1 + exp(lpsiD[site_p,age_p,1]) + exp(lpsiD[site_p,age_p,2])+ exp(lpsiD[site_p,age_p,3]))
    
    }
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)
    
#     Sex[1] <- 0 
#     Sex[2] ~ dnorm(0, 0.1)I(-10,10)
        Anneei[1] <- 0 
        Anneei[2] ~ dnorm(0, 0.1)I(-10,10)
    
    for(agei in 1:2)  {
    for (sitei in 1:2){
    gamma[agei,sitei]~ dnorm(0, 0.0001)
    } 
    }
    #     S.Age1.Male.Annee1<-1/(1+exp(-gamma[1]))
    #     S.Age1.Male.Annee2<-1/(1+exp(-gamma[1]-Annee[2]))
    #     S.Age1.Femelle.Annee1<-1/(1+exp(-gamma[1]-Sex[2]))
    #     S.Age1.Femelle.Annee2<-1/(1+exp(-gamma[1]-Sex[2]-Annee[2]))
    #     S.Age2.Male.Annee1<-1/(1+exp(-gamma[2]))
    #     S.Age2.Male.Annee2<-1/(1+exp(-gamma[2]-Annee[2]))
    #     S.Age2.Femelle.Annee1<-1/(1+exp(-gamma[2]-Sex[2]))
    #     S.Age2.Femelle.Annee2<-1/(1+exp(-gamma[2]-Sex[2]-Annee[2]))
    #     
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    
    
    
    logit(SS[i,t])<-gamma[xage[i,t],site[i]]+Anneei[annee[i]]
    
    
    
       ps[1,i,t,1] <- SS[i,t] * psiA[site[i],xage[i,t],1]
    ps[1,i,t,2] <- SS[i,t] * psiA[site[i],xage[i,t],2]
    ps[1,i,t,3] <- SS[i,t] * psiA[site[i],xage[i,t],3]
    ps[1,i,t,4] <- SS[i,t] * psiA[site[i],xage[i,t],4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[site[i],xage[i,t],1]
    ps[2,i,t,2] <- SS[i,t] * psiB[site[i],xage[i,t],2]
    ps[2,i,t,3] <- SS[i,t] * psiB[site[i],xage[i,t],3]
    ps[2,i,t,4] <- SS[i,t] * psiB[site[i],xage[i,t],4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[site[i],xage[i,t],1]
    ps[3,i,t,2] <- SS[i,t] * psiC[site[i],xage[i,t],2]
    ps[3,i,t,3] <- SS[i,t] * psiC[site[i],xage[i,t],3]
    ps[3,i,t,4] <- SS[i,t] * psiC[site[i],xage[i,t],4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[site[i],xage[i,t],1]
    ps[4,i,t,2] <- SS[i,t] * psiD[site[i],xage[i,t],2]
    ps[4,i,t,3] <- SS[i,t] * psiD[site[i],xage[i,t],3]
    ps[4,i,t,4] <- SS[i,t] * psiD[site[i],xage[i,t],4]
    ps[4,i,t,5] <- (1-SS[i,t])*r
    ps[4,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 0
    ps[5,i,t,6] <- 1
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    po[1,i,t,6] <- 1-p
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- p
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 0
    po[2,i,t,6] <- 1-p 
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- p
    po[3,i,t,4] <- 0 
    po[3,i,t,5] <- 0
    po[3,i,t,6] <- 1-p
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- p
    po[4,i,t,5] <- 0  
    po[4,i,t,6] <- 1-p 
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1  
    po[5,i,t,6] <- 0
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 0  
    po[6,i,t,6] <- 1 
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data

known.state <- rCH
known.state[known.state==6] <- NA
for (i in 1:dim(known.state)[1]){
  m <- min(which(!is.na(known.state[i,])))
  known.state[i,m] <- NA
  state2=known.state
  state2[which(is.na(state2))]=0  
  for (j in 1:dim(known.state)[2]){
    if (state2[i,j]== 5) known.state[i,(j+1):(dim(known.state)[2])]= 6 
  }
}

init.state<-known.state 
init.state[which(is.na(init.state))] <- 0
for (i in 1:dim(init.state)[1]) { init.state[i,1:f[i]] <- NA}
states <- max(init.state, na.rm = TRUE)
usable_states <- 1:(states-2)
init.state[which(init.state>0)] <- NA
init.state[which(init.state==0)] <- sample(usable_states,length(which(init.state==0)),replace = TRUE)



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,xage=xage,sex=sex,site=site,annee=annee)#,rel_hs=rel_hs

inits <- function(){list(gamma=array(rnorm(4),c(2,2)),Annee=c(NA,rnorm(1)),Site=c(NA,rnorm(1)), r= runif(1, 0, 1),lpsiA = array(rnorm(64),dim=c(2,2,3)), lpsiB =array(rnorm(64),dim=c(2,2,3)), lpsiC = array(rnorm(64),dim=c(2,2,3)), lpsiD = array(rnorm(64),dim=c(2,2,3)), p= runif(1, 0, 1),z = init.state)}# , lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3) 

#


# Parameters monitored
parameters <-  c( "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1


# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

save.image(paste(a,"survie_sans_HS.Rdata",sep=""))
}
