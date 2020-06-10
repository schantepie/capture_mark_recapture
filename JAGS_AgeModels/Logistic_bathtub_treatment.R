  # Specify model in BUGS language
  
  setwd("/media/DATA/PROJETS_ANNEXES/EcolLet/analysesMarked/")
  data=read.csv(file = "/media/DATA/PROJETS_ANNEXES/EcolLet/inputMat-Basta.csv",row.names =1)
  
  # for(i in 1:(dim(data)[1])){
  #   x=data[i,which(grepl("X",colnames(data)))] 
  #   reco=which(x==1)
  #   maskAge=rep(0, dim(x)[2])
  #   maskAge[min(reco): dim(x)[2]]=1:(dim(x)[2]- (min(reco)-1))
  #   maskAge=maskAge*x
  #   data[i,which(grepl("X",colnames(data)))]=maskAge
  # }
  
  X1983=rep(0, dim(data)[1])
  rCH=as.matrix(cbind(X1983,data[,which(grepl("X",colnames(data)))]))
  for (i in 1:dim(data)[1]) rCH[i,paste("X",data$Birth[i],sep="")]=1
  
  get.first <- function(x) min(which(x!=0))
  f <- apply(rCH, 1, get.first)
  
  age <- matrix(NA, ncol = dim(rCH)[2], nrow = dim(rCH)[1])
  for (i in 1:dim(rCH)[1]){
    for (t in f[i]:(dim(rCH)[2])){
      age[i,t] <- t-f[i]
    } #t
  } #i
  age= age[,-1]
  
  rCH=as.matrix(data[,which(grepl("X",colnames(data)))])
  f <- apply(rCH, 1, get.first)
  
  
  for (i in 1:dim(rCH)[1]) {
    if (f[i]!=1)    age[i,1:(f[i]-1)]=NA
  }
  
 
  ch=rCH
  i=1
  known.state.cjs <- function(ch){
    state <- ch
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
    }
    state[state==0] <- NA
    return(state)
  }
  
  Treatment=data[,which(grepl("Treatment",colnames(data)))]
  treatment=Treatment$Treatmentc
  treatment[which(Treatment$Treatmentc==1)]=1
  treatment[which(Treatment$Treatmentm2==1)]=2
  treatment[which(Treatment$Treatmentp2==1)]=3
  
  
  
  # Initial values
  cjs.init.z <- function(ch,f){
    for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
    }
    for (i in 1:dim(ch)[1]){
      ch[i,1:f[i]] <- NA
    }
    return(ch)
  }
  
  
  require(jagsUI)
  sink("cjs-age.jags")
  cat("
  model {
  
  # Priors and constraints
  
  mean.p ~ dunif(0, 1)
  
  # Contraints on phi and p
  

for (u in 1:treat){
  b0[u]~ dunif(0, 10) 
  b1[u] ~ dunif(0, 10)
  db0[u] ~ dunif(0, 1) 
  db1[u] ~ dunif(0, 5)
  s[u] ~ dunif(0, 10) 
  c[u] ~ dunif(0, 1) 
}


  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
         phi[i,t] <- exp((exp(1)*log(b0[treatment[i]]*s[treatment[i]]*exp(age[i,t]*b1[treatment[i]])+exp(1)*b1[treatment[i]]))/s[treatment[i]]-(exp(1)*log(b0[treatment[i]]*s[treatment[i]]*exp(age[i,t+1]*b1[treatment[i]])+exp(1)*b1[treatment[i]]))/s[treatment[i]]-(db0[treatment[i]]*(exp(db1[treatment[i]])-1)*exp(-age[i,t+1]*db1[treatment[i]])+c[treatment[i]]*db1[treatment[i]])/db1[treatment[i]])
         p[i,t] <- mean.p
    } #t
  } #i
  

  # Retrieve classic values
  for(age in 0:7){
      for(tre in 1:treat){
    phi.age[age+1,tre]  <-exp((exp(1)*log(b0[tre]*s[tre]*exp(age*b1[tre])+exp(1)*b1[tre]))/s[tre]-(exp(1)*log(b0[tre]*s[tre]*exp((age+1)*b1[tre])+exp(1)*b1[tre]))/s[tre]-(db0[tre]*(exp(db1[tre])-1)*exp(-(age+1)*db1[tre])+c[tre]*db1[tre])/db1[tre])
    instantMorta.age[age+1,tre]  <- (b0[tre]*exp(b1[tre]*age))/((b0[tre]*s[tre]*exp(b1[tre]*age-1))/b1[tre]+1)+c[tre]+db0[tre]*exp(-db1[tre]*age)
    cumulSurv.age[age+1,tre] <- exp(-((exp(1)*(log(b0[tre]*s[tre]*exp((age+1)*b1[tre])+exp(1)*b1[tre])-log(b0[tre]*s[tre]+exp(1)*b1[tre])))/s[tre])-age*c[tre]+(db0[tre]*(exp(-(age+1)*db1[tre])-1))/db1[tre]-c[tre])
  }
}
  
  
    
  # Likelihood
  for (i in 1:nind){
  
  # Define latent state at first capture
  z[i,f[i]] <- 1
  
  for (t in (f[i]+1):n.occasions){
  # State process
  ProbaAlive[i,t] <- phi[i,t-1] * z[i,t-1]
  z[i,t] ~ dbern(ProbaAlive[i,t])
  
  # Observation process
  ProbaSeen[i,t] <- p[i,t-1] * z[i,t]
  y[i,t] ~ dbern(ProbaSeen[i,t])
  
  } #t
  } #i
  }
  ",fill = TRUE)
  sink()
  
  
  
  
  # Bundle data
  bugs.data <- list(y = rCH, f = f, nind = dim(rCH)[1], n.occasions =
                    dim(rCH)[2], z = known.state.cjs(rCH), age = age,
                    treat=length(unique(treatment)),treatment = treatment)
  
  inits <- function(){list(z = cjs.init.z(rCH, f),
                           b0 = runif(3, 0, 10),
                           b1 = runif(3, 0, 10),
                           db0 = runif(3,0, 1) ,
                           db1 = runif(3,0, 5),
                           s = runif(3,0, 10) ,
                           c = runif(3,0, 1) ,
                           mean.p = runif(1, 0, 1))}
  # Parameters monitored
  parameters <- c("b0", "b1","db0","db0","phi.age","instantMorta.age","cumulSurv.age")
  # MCMC settings
  ni <- 110000
  nt <- 100
  nb <- 10000
  nc <- 4
  # Call WinBUGS from R (BRT 3 min)

  cjs <-jags(bugs.data, inits, parameters, "cjs-age.jags", n.chains =nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, n.cores = nc)
  
  # plot(cjs.age$BUGSoutput$summary[4:10,1])
  # 
  # 

   print(cjs, digits = 3)
   plot(cjs$mean$instantMorta.age[,1])
   
   
  #     
  # age=1:7
  # a=0.5
  # b=0.06
  # 
  # exp(-a*exp(b*age))
  #     
  
