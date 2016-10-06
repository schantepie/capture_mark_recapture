
library(R2jags)

setwd("/mnt/travail/6_Papiers_presentations/prep_Annechri_outard")
HDVW=read.csv("/mnt/travail/6_Papiers_presentations/prep_Annechri_outard/HDV_by_week_290714.csv", sep=";")

cp=HDVW[1:(nrow(HDVW)-9),2:86]


CH=matrix(0,nrow=nrow(cp),ncol=ncol(cp)+1)
dim(CH)
for (i in 1:dim(cp)[1]){
  for (j in 1:dim(cp)[2]){
    po=strsplit(as.character(cp[i,j]),split="")
    if (po[[1]][1]=="A")     CH[i,j]=1
    if (po[[1]][1]=="B")     CH[i,j]=2
    if (po[[1]][1]=="C")     CH[i,j]=3
    if (po[[1]][1]=="D")     CH[i,j]=4
    if (po[[1]][2]=="1")  {  CH[i,(j+1)]= 5 }
    
  } 
}
rel_hs=HDVW$REL_HS[1:(nrow(HDVW)-9)]

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

rCH <- CH  # Recoded CH
rCH[rCH==0] <- 6



rel_hs2=HDVW$REL_HS[1:(nrow(HDVW)-9)]
rel_hs2=(rel_hs2-mean(rel_hs2))/sd(rel_hs2)


xage <- matrix(NA, ncol = dim(rCH)[2]-1, nrow = dim(rCH)[1])

for (i in 1:dim(rCH)[1]){
  for (t in f[i]:(dim(rCH)[2]-1)){
    xage[i,t] <- 2
    xage[i,f[i]:(f[i]+4)] <- 1   #f[i]:(f[i]+4) signifie 5 semaines
  } #t
} #i

site=as.character(HDVW$SITE[1:180])
site[site=="ORIENTAL"]=1
site[site=="TATA"]=2

sink("ms3-multinomlogit.jags")
cat("
    model {
    
   
      for (i in 1:4){
       psiA[i] <-1
       psiB[i] <-1
       psiC[i] <-1
      psiD[i] <-1
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)

for(site in 1:2){
      for(age in 1:2){
        int[site,age]<-log(mean.phi[site,age]/(1-mean.phi[site,age]))
        mean.phi[site,age]~ dunif(0,1)
        pente[site,age]~ dnorm(0, 0.001)
    }
}
      for (i in 1:nind){
          epsilon[i]~dnorm(0,tau)
      }
      tau<-pow(sigma,-2)
      sigma~dunif(0,5)

    # Define state-transition and observation matrices   
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
  
    logit(SS[i,t])<-int[site[i],xage[i,t]]+(pente[site[i],xage[i,t]]*rel_hs2[i])+epsilon[i]


  
    ps[1,i,t,1] <- SS[i,t] * psiA[1]
    ps[1,i,t,2] <- SS[i,t] * psiA[2]
    ps[1,i,t,3] <- SS[i,t] * psiA[3]
    ps[1,i,t,4] <- SS[i,t] * psiA[4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[1]
    ps[2,i,t,2] <- SS[i,t] * psiB[2]
    ps[2,i,t,3] <- SS[i,t] * psiB[3]
    ps[2,i,t,4] <- SS[i,t] * psiB[4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[1]
    ps[3,i,t,2] <- SS[i,t] * psiC[2]
    ps[3,i,t,3] <- SS[i,t] * psiC[3]
    ps[3,i,t,4] <- SS[i,t] * psiC[4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[1]
    ps[4,i,t,2] <- SS[i,t] * psiD[2]
    ps[4,i,t,3] <- SS[i,t] * psiD[3]
    ps[4,i,t,4] <- SS[i,t] * psiD[4]
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



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,rel_hs2=rel_hs2,xage=xage,site=site)#,rel_hs=rel_hs

inits <- function(){list(mean.phi= array(runif(4, 0, 1),dim=c(2,2)),pente=array(rnorm(64),dim=c(2,2)), r= runif(1, 0, 1),z = init.state,sigma=runif(1,0,10))}  

#, lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3)


# Parameters monitored
parameters <-  c("mean.phi","int","pente", "psiA", "psiB", "psiC", "psiD", "r","p","epsilon")

# MCMC settings
ni <- 110000
nt <- 100
nb <- 10000
nc <- 1

# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image("HB.Rdata")





library(R2jags)

setwd("/mnt/travail/6_Papiers_presentations/prep_Annechri_outard")
HDVW=read.csv("/mnt/travail/6_Papiers_presentations/prep_Annechri_outard/HDV_by_week_290714.csv", sep=";")

cp=HDVW[1:(nrow(HDVW)-9),2:86]


CH=matrix(0,nrow=nrow(cp),ncol=ncol(cp)+1)
dim(CH)
for (i in 1:dim(cp)[1]){
  for (j in 1:dim(cp)[2]){
    po=strsplit(as.character(cp[i,j]),split="")
    if (po[[1]][1]=="A")     CH[i,j]=1
    if (po[[1]][1]=="B")     CH[i,j]=2
    if (po[[1]][1]=="C")     CH[i,j]=3
    if (po[[1]][1]=="D")     CH[i,j]=4
    if (po[[1]][2]=="1")  {  CH[i,(j+1)]= 5 }
    
  } 
}

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

rCH <- CH  # Recoded CH
rCH[rCH==0] <- 6

HS_cont= read.csv("/mnt/travail/6_Papiers_presentations/prep_Annechri_outard/HDVW_matHS_290714_2.csv", sep=";",h=T,row.names=1)


xage <- matrix(NA, ncol = dim(rCH)[2]-1, nrow = dim(rCH)[1])

for (i in 1:dim(rCH)[1]){
  for (t in f[i]:(dim(rCH)[2]-1)){
    xage[i,t] <- 2
    xage[i,f[i]:(f[i]+4)] <- 1   #f[i]:(f[i]+4) signifie 5 semaines
  } #t
} #i


site=as.character(HDVW$SITE[1:180])
site[site=="ORIENTAL"]=1
site[site=="TATA"]=2

sink("ms3-multinomlogit.jags")
cat("
    model {
    
         for (i in 1:4){
       psiA[i] <-1
       psiB[i] <-1
       psiC[i] <-1
      psiD[i] <-1
    }
    r ~ dunif(0, 1)
    p ~ dunif(0, 1)

for(site in 1:2){
    for(age in 1:2)  {
        int[site,age]~ dnorm(0, 0.001)
        pente[site,age]~ dnorm(0, 0.001)

       for(t in 1:(n.occasions-1)){
         mu[site,age,t] ~ dnorm(mu.mu[site,age],mu.tau[site,age])
       }
        mu.mu[site,age] ~ dnorm(0,.001)
        mu.tau[site,age]~ dgamma(.001,.001)
        tau[site,age] ~ dgamma(.001,.001)
  }
}
    
    # Define state-transition and observation matrices   
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
  
    HS_cont[i,f[i]] ~ dnorm(0,.0001)

    for (t in f[i]+1:(n.occasions-1)){
    mu.z[site[i],i,t] <- HS_cont[i,t-1] + mu[site[i],xage[i,t-1],t-1]
    HS_cont[i,t] ~ dnorm(mu.z[site[i],i,t],tau[site[i],xage[i,t-1]])
    }

    for (t in f[i]:(n.occasions-1)){
  
    logit(SS[i,t])<-int[site[i],xage[i,t]]+(pente[site[i],xage[i,t]]*HS_cont[i,t])
  
    ps[1,i,t,1] <- SS[i,t] * psiA[1]
    ps[1,i,t,2] <- SS[i,t] * psiA[2]
    ps[1,i,t,3] <- SS[i,t] * psiA[3]
    ps[1,i,t,4] <- SS[i,t] * psiA[4]
    ps[1,i,t,5] <- (1-SS[i,t])*r
    ps[1,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[2,i,t,1] <- SS[i,t] * psiB[1]
    ps[2,i,t,2] <- SS[i,t] * psiB[2]
    ps[2,i,t,3] <- SS[i,t] * psiB[3]
    ps[2,i,t,4] <- SS[i,t] * psiB[4]
    ps[2,i,t,5] <- (1-SS[i,t])*r
    ps[2,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[3,i,t,1] <- SS[i,t] * psiC[1]
    ps[3,i,t,2] <- SS[i,t] * psiC[2]
    ps[3,i,t,3] <- SS[i,t] * psiC[3]
    ps[3,i,t,4] <- SS[i,t] * psiC[4]
    ps[3,i,t,5] <- (1-SS[i,t])*r
    ps[3,i,t,6] <- (1-SS[i,t])*(1-r)
    ps[4,i,t,1] <- SS[i,t] * psiD[1]
    ps[4,i,t,2] <- SS[i,t] * psiD[2]
    ps[4,i,t,3] <- SS[i,t] * psiD[3]
    ps[4,i,t,4] <- SS[i,t] * psiD[4]
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



jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1],z = known.state,HS_cont=HS_cont,xage=xage,site=site)#,rel_hs=rel_hs

inits <- function(){list(int=array(rnorm(64),dim=c(2,2)),pente=array(rnorm(64),dim=c(2,2)) ,r= runif(1, 0, 1),z = init.state)}  

#, lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3)


# Parameters monitored
parameters <-  c("int","pente", "psiA", "psiB", "psiC", "psiD", "r","p")

# MCMC settings

ni <- 110000
nt <- 100
nb <- 10000
nc <- 1
# Call JAGS from R (BRT 80 min)
outdead<- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
save.image("HC.Rdata")
