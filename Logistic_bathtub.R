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
  
  b0 ~ dunif(0, 10) 
  b1 ~ dunif(0, 10)
  db0 ~ dunif(0, 1) 
  db1 ~ dunif(0, 5)
  s ~ dunif(0, 10) 
  c ~ dunif(0, 1) 
  mean.p ~ dunif(0, 1)
  
  # Contraints on phi and p
  
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
         phi[i,t] <-    exp((exp(1)*log(b0*s*exp(age[i,t]*b1)+exp(1)*b1))/s-(exp(1)*log(b0*s*exp(age[i,t+1]*b1)+exp(1)*b1))/s-(db0*(exp(db1)-1)*exp(-age[i,t+1]*db1)+c*db1)/db1)
         p[i,t] <- mean.p
    } #t
  } #i
  

  # Retrive vclassic values
  for(age in 0:7){
    phi.age[age+1]  <-exp((exp(1)*log(b0*s*exp(age*b1)+exp(1)*b1))/s-(exp(1)*log(b0*s*exp((age+1)*b1)+exp(1)*b1))/s-(db0*(exp(db1)-1)*exp(-(age+1)*db1)+c*db1)/db1)
    instantMorta.age[age+1]  <- (b0*exp(b1*age))/((b0*s*exp(b1*age-1))/b1+1)+c+db0*exp(-db1*age)
    cumulSurv.age[age+1] <- exp(-((exp(1)*(log(b0*s*exp((age+1)*b1)+exp(1)*b1)-log(b0*s+exp(1)*b1)))/s)-age*c+(db0*(exp(-(age+1)*db1)-1))/db1-c)
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
                           b0 = runif(1, 0, 10),
                           b1 = runif(1, 0, 10),
                           db0 = runif(1,0, 1) ,
                           db1 = runif(1,0, 5),
                           s = runif(1,0, 10) ,
                           c = runif(1,0, 1) ,
                           mean.p = runif(1, 0, 1))}
  # Parameters monitored
  parameters <- c("b0", "b1","phi.age","instantMorta.age","cumulSurv.age")
  # MCMC settings
  ni <- 2000
  nt <- 1
  nb <- 10
  nc <- 1
  # Call WinBUGS from R (BRT 3 min)
  cjs.age <-jags(bugs.data, inits, parameters, "cjs-age.jags", n.chains =nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  
    
  # plot(cjs.age$BUGSoutput$summary[4:10,1])
  # 
  # 

   print(cjs.age, digits = 3)
   plot(cjs.age$BUGSoutput$summary[grep("phi",rownames(cjs.age$BUGSoutput$summary)),1])
   
   
  #     
  # age=1:7
  # a=0.5
  # b=0.06
  # 
  # exp(-a*exp(b*age))
  #     
  
