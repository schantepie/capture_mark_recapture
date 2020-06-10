  # Specify model in BUGS language
  
  setwd("/mnt/DATA/RESEARCH/project/Celine_flyCatcher")
  data=read.csv(file = "inputMat-Basta.csv",row.names =1)
  
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
  
  
  require(R2jags)
  sink("cjs-age.jags")
  cat("
  model {
  
  # Priors and constraints
  
  b0 ~ dunif(0, 1) # Prior for mean of logit survival
  b1 ~ dunif(0, 5)# Prior for slope parameter
  db0 ~ dunif(0, 1) # Prior for mean of logit survival
  db1 ~ dunif(0, 5)# Prior for slope parameter
  c ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  
  # Contraints on phi and p
  
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
         phi[i,t] <-  exp(((b0/b1) *(exp(age[i,t]*b1)-exp(age[i,t+1]*b1)))-c-((db0/db1) * exp(-age[i,t+1]*db1)*(-1+exp(db1))))
         p[i,t] <- mean.p
    } #t
  } #i
  
  # Retrive vclassic values
  for(aging in 0:7){
    phi.age[aging+1]  <-  exp(((b0/b1) *(exp(aging*b1)-exp((aging+1)*b1)))-c-((db0/db1) * exp(-(aging+1)*db1)*(-1+exp(db1))))
    instantMorta.age[aging+1]  <- (db0*exp(-db1*aging)) + c + (b0*exp(b1*aging))
    cumulSurv.age[aging+1] <- exp(((b0/b1) *(1-exp((aging)*b1)))-(aging*c)-((db0-(db0*exp(-aging*db1)))/db1))
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
                      dim(rCH)[2], z = known.state.cjs(rCH), age = age)
  
  inits <- function(){list(z = cjs.init.z(rCH, f),
                           b0 = runif(1, 0, 1),
                           b1 = runif(1, 0, 1),
                           c = runif(1, 0, 1),
                           mean.p = runif(1, 0, 1))}
  # Parameters monitored
  parameters <- c("b0", "b1","db0", "db1","c","phi.age","instantMorta.age","cumulSurv.age")
  # MCMC settings
  ni <- 2000
  nt <- 1
  nb <- 100
  nc <- 3
  # Call WinBUGS from R (BRT 3 min)
  cjs.age <-jags(bugs.data, inits, parameters, "cjs-age.jags", n.chains =nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  
    
  # plot(cjs.age$BUGSoutput$summary[4:10,1])
  # 
  # 

   print(cjs.age, digits = 3)
  #     

   plot(0.8* exp(- c(0:10)*0.5))
