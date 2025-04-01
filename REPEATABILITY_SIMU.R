########################
# Set up the environment
########################

rm(list=ls())
library(parallel)


###################################
## Functions
###################################


# Return CI based on Fisher Z transform
#######################################

Z.transform<-function(point,se){
  z_point<-0.5*log((1+point)/(1-point))
  z_se<-sqrt(se^2/(1-point^2)^2)
  z_lower<-z_point-1.96*z_se
  z_upper<-z_point+1.96*z_se
  lower<-(exp(2*z_lower)-1)/(exp(2*z_lower)+1)
  upper<-(exp(2*z_upper)-1)/(exp(2*z_upper)+1)
  
  return(list("z.low"=lower,"z.hig"=upper))
}


# Return median and CI based on DM (K.B)
#######################################

kappa.bayes.bin<-function(data){
  
  n<-nrow(data)
  nrat<-ncol(data)
  npair<-nrat*(nrat-1)/2
  
  number<-1/(2^ncol(data))
  
  l <- rep(list(factor(1:2,levels=c(1,2),labels=c(1,2))), nrat)
  pi.label<-expand.grid(l)
  
  data.f_<-as.data.frame.table(table(data.frame(data))/n)
  names(pi.label)<-names(data.f_)[1:nrat]
  all.pi<-as.data.frame(pi.label)
  
  data.f<-merge(data.f_,as.data.frame(pi.label),all.y=TRUE)
  data.f[is.na(data.f)] <- 0
  
  rat.pair<-t(combn(seq(1,nrat), 2))
  
  fre<-data.f[,nrat+1]
  
  index<-seq(1,length(fre))
  
  ind_1<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==1])
  ind_2<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==2])
  ind_3<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==1 & data.f[,rat.pair[x,2]]==2])
  ind_4<-sapply(1:npair,function(x) index[data.f[,rat.pair[x,1]]==2 & data.f[,rat.pair[x,2]]==1])
  
  ind<-sapply(1:nrat,function(x) index[data.f[,x]==1])
  
  N2 <-30000
  
  alpha<-fre*n+number
  
  w<- lapply(alpha, function(x) rgamma(N2,shape=x,rate=1))
  
  t<-Reduce("+", w) 
  
  z<-sapply(w, function(x) x/t)
  
  if (npair>1){
    
    p11<-rowSums(z[,c(ind_1)])/(npair)
    p22<-rowSums(z[,c(ind_2)])/(npair)
    p12<-rowSums(z[,c(ind_3)])/(npair)
    p21<-rowSums(z[,c(ind_4)])/(npair)
    
    po<-p11+p22
    
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-((marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2+(1-(marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2
    }
    
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
    
  }
  
  if (npair==1){
    p11<-z[,c(ind_1)]/(npair)
    p22<-z[,c(ind_2)]/(npair)
    p12<-z[,c(ind_3)]/(npair)
    p21<-z[,c(ind_4)]/(npair)
    
    po<-p11+p22
    
    marge_<-sapply(1:nrat,function(x) rowSums(z[,c(ind[,x])]))  
    
    pe_pair<-matrix(NA,nrow=N2,ncol=npair)
    
    for (b in 1:npair){
      pe_pair[,b]<-((marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2+(1-(marge_[,rat.pair[b,1]]+marge_[,rat.pair[b,2]])/2)^2
    }
    pe<-rowMeans(pe_pair)
    
    kappa<-(po-pe)/(1-pe)
  }
  
  q.kappa<-quantile(kappa,c(0.025,0.50,0.975),na.rm = TRUE)
  m.kappa<-mean(kappa)
  
  
  
  return(list("ICC"=unname(unlist(q.kappa[2])),
              "LOW"=unname(unlist(q.kappa[1])),
              "HIGH"=unname(unlist(q.kappa[3]))
  ))

}

# Return kappa coefficient with delta CI (K.F, K.FZ)
#####################################################

kappaw<-function(data,weight,chance)
{
  check<-all(sweep(data,MARGIN=1,STATS=data[,1],FUN="=="))
  
  if (check==TRUE){
    return(list("ICC"=1,
                "LOW"=1,
                "HIGH"=1,
                "LOW.Z"=1,
                "HIGH.Z"=1
    ))
  }
  if (check==FALSE){
    n<-nrow(data)
    nrat<-ncol(data)
    npair<-nrat*(nrat-1)/2
    if (min(data==0)){data<-data+1}
    all<-c(data.matrix(data))
    all.f<-factor(all,levels=names(table(all)))
    ncat<-nlevels(all.f)    		#number of categories of the scale
    score <- matrix(1:ncat, nrow = ncat, ncol = ncat, byrow = TRUE)
    
    if (weight=="binary"){
      w<-diag(ncat)
    }
    if (weight=="linear"){
      w<-1-abs(score - t(score))/(ncat-1)
    }
    if (weight=="quadratic"){
      w<-1-(score - t(score))^2/(ncat-1)^2
    }
    
    freq<-t(sapply(1:n,fun2<-function(i){ table(factor(data[i,],levels=names(table(all))))}))
    
    p<-matrix(NA,nrow=ncat,ncol=ncat)
    oprim<-rep(NA,n)
    for (h in 1:n){
      for (i in 1:ncat){
        for (j in 1:ncat){
          if (i==j){
            p[i,j]<-w[i,j]*freq[h,i]*(freq[h,i]-1)/(nrat*(nrat-1))
          }
          if (i!=j){
            p[i,j]<-w[i,j]*freq[h,i]*freq[h,j]/(nrat*(nrat-1))
          }
        }
        
      }
      oprim[h]<-sum(p)
    }
    
    P_o<-mean(oprim)
    
    
    if (chance==1){
      #chance=1
      eprim<-rep(NA,n)
      p_j<-rep(1/ncat,ncat)
      eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
      eprim<-rep(NA,n)
      
      for (h in 1:n){
        for (r in 1:nrat){
          temp<-rep(NA,nrat-1)
          l<-1
          for (s in setdiff(1:nrat,r)){
            temp[l]<-sum(p_j*w[,data[h,s]])
            l<-l+1
          }
          eprim_tmp[h,r]<-sum(temp)
        }
        eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
      }
      
      P_e<-mean(eprim)
      
      d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
      d2<-(P_o*P_e-2*P_e+P_o)^2
      var<-(d1-d2)/(n*(1-P_e)^4)
      kappa<-(P_o-P_e)/(1-P_e)
    }
    
    if (chance==3){
      #chance=3
      eprim<-rep(NA,n)
      p_j<-colMeans(sweep(freq,1,nrat,"/"))
      eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
      eprim<-rep(NA,n)
      
      for (h in 1:n){
        for (r in 1:nrat){
          temp<-rep(NA,nrat-1)
          l<-1
          for (s in seq(1,nrat)[-r]){
            temp[l]<-sum(p_j*w[,data[h,s]])
            l<-l+1
          }
          eprim_tmp[h,r]<-sum(temp)
        }
        eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
      }
      
      P_e<-mean(eprim)
      
      d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
      d2<-(P_o*P_e-2*P_e+P_o)^2
      var<-(d1-d2)/(n*(1-P_e)^4)
      kappa<-(P_o-P_e)/(1-P_e)
    }
    
    
    if (chance==2){
      #chance=2
      eprim<-rep(NA,n)
      M<-replicate(nrat,matrix(0,ncol=ncat,nrow=1), simplify=FALSE)
      M<-lapply(1:nrat,f1<-function(idd){prop.table(table(factor(data[,idd],levels=names(table(all)))))})
      eprim_tmp<-matrix(NA,nrow=n,ncol=nrat)
      eprim<-rep(NA,n)
      
      for (h in 1:n){
        for (r in 1:nrat){
          temp<-rep(NA,nrat-1)
          l<-1
          for (s in seq(1,nrat)[-r]){
            temp[l]<-sum(M[[r]]*w[,data[h,s]])
            l<-l+1
          }
          eprim_tmp[h,r]<-sum(temp)
        }
        eprim[h]<-sum(eprim_tmp[h,])/(nrat*(nrat-1))
      }
      
      P_e<-mean(eprim)
      
      d1<-sum(((1-P_e)*oprim-2*(1-P_o)*eprim)^2)/n
      d2<-(P_o*P_e-2*P_e+P_o)^2
      var<-(d1-d2)/(n*(1-P_e)^4)
      kappa<-(P_o-P_e)/(1-P_e)
    }
    
    z.kap<-Z.transform(kappa,sqrt(var))
    
    return(list("ICC"=kappa,
                "LOW"=kappa-1.96*sqrt(var),
                "HIGH"=kappa+1.96*sqrt(var),
                "LOW.Z"=unlist(unname(z.kap[1])),
                "HIGH.Z"=unlist(unname(z.kap[2]))
    ))
  }
}


# ICC to tetrachoric (Bonnett approximation)
#######################################
tetra_app<-function(phi,pmarge){
  a11<-phi*pmarge*(1-pmarge)+pmarge^2
  a12<-a21<-pmarge-a11
  a22<-1-a11-a12-a21
  w<-a11*a22/(a12*a21)
  c<-(1-(0.5-pmarge)^2)/2 
  rho_app<-cos(pi/(1+w^c))
  return(rho_app)
}

# Tetrachoric to phi (Bonnett approximation)
#######################################

phi_app<-function(tetra,pmarge){
  c<-(1-(0.5-pmarge)^2)/2 
  w<-((pi-acos(tetra))/acos(tetra))^(1/c)
  
  b_<-(2*pmarge*w-2*pmarge+1)/(2*(w-1))
  sq_r<-sqrt(-4*pmarge^2*w+4*pmarge^2+4*pmarge*w-4*pmarge+1)/(2*(w-1))
  a11a<-b_+sq_r
  a11b<-b_-sq_r
  
  a11<-ifelse(a11a>pmarge|a11a<0,ifelse(a11b<pmarge&a11b>0,a11b,0),a11a)
  
  phi<-(a11-pmarge^2)/(pmarge*(1-pmarge))
  
  return(phi)
}



# simulation of data
#######################################

simu<-function(prop,n.rat,ICC.sim,n.sub){

  rho<-tetra_app(phi=ICC.sim,pmarge=prop)

  var.sub<-rho/(1-rho)
  e_ij<-matrix(rnorm(n.sub*n.rat, sd = 1), nrow=n.sub, ncol=n.rat)
  s_i<-matrix(rep(rnorm(n.sub, sd=sqrt(var.sub)),n.rat), nrow=n.sub)
  data.sim<-ifelse((s_i+e_ij)>-qnorm(prop,0,sqrt(var.sub+1)),1,2)
  return(list("data.sim"=data.sim))

}


# Probit - Bayesian (P.BI)
##########################

cat("
model {
  for (i in 1:n) {
    y[i] ~ dbern(p[i])
    p[i] <- phi(mu[i])
    mu[i]<-beta0+S[sub[i]]
  }

  beta0 ~ dnorm(0,1)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(0,tau_S)  
  }
 
 # Priors for variance components
 tau_S<-1/sigma2_S
 sigma2_S<-ICC/(1-ICC)
 ICC~dunif(0,1)
 
 beta0_marge<-beta0/(sqrt(1+sigma2_S))
 pmarge<-phi(beta0_marge)
 
 c<-(1-(0.5-pmarge)^2)/2 
 w<-((3.14159-acos(ICC))/acos(ICC))^(1/c)
  
 b_<-(2*pmarge*w-2*pmarge+1)/(2*(w-1))
 sq_r<-sqrt(-4*pmarge^2*w+4*pmarge^2+4*pmarge*w-4*pmarge+1)/(2*(w-1))
 a11a<-b_+sq_r
 a11b<-b_-sq_r
  
 a11<-ifelse(a11a>pmarge||a11a<0,ifelse(a11b<pmarge&&a11b>0,a11b,0),a11a)
 phi<-(a11-pmarge^2)/(pmarge*(1-pmarge))

    }",file="probit1.txt")


probit.b1<-function(data.sim){

  #Transform the data into long dataset
  n<-nrow(data.sim)
  nrat<-ncol(data.sim)
  
  y<-c(data.sim) 
  sub<-c(rep(seq(1,n),nrat))
  
  data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
  parameters<-c("ICC","phi")
  jags.inits<-list(beta0=0.0,ICC=0.7)
  obj <- jags.model("probit1.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model3.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  model3.sim<-summary(model3.s)
  
  return(list("ICC.LAT"=model3.sim[[2]][1,3],
              "LOW.LAT"=model3.sim[[2]][1,1],
              "HIGH.LAT"=model3.sim[[2]][1,5],
              "ICC"=model3.sim[[2]][2,3],
              "LOW"=model3.sim[[2]][2,1],
              "HIGH"=model3.sim[[2]][2,5])) 
  
}





# Probit - Bayesian (P.BU)
##########################
cat("
model {
  for (i in 1:n) {
    y[i] ~ dbern(p[i])
    p[i] <- phi(mu[i])
    mu[i]<-beta0+S[sub[i]]
  }

  beta0 ~ dnorm(0,1)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(0,tau_S)  
  }
 
 # Priors for variance components
 tau_S<-1/sigma2_S
 sigma2_S<-pow(sigma_S,2)

 sigma_S~dunif(0,30)
 ICC<-sigma2_S/(sigma2_S+1)
 
 beta0_marge<-beta0/(sqrt(1+sigma2_S))
 pmarge<-phi(beta0_marge)
 
 c<-(1-(0.5-pmarge)^2)/2 
 w<-((3.14159-acos(ICC))/acos(ICC))^(1/c)
  
  b_<-(2*pmarge*w-2*pmarge+1)/(2*(w-1))
  sq_r<-sqrt(-4*pmarge^2*w+4*pmarge^2+4*pmarge*w-4*pmarge+1)/(2*(w-1))
  a11a<-b_+sq_r
  a11b<-b_-sq_r
 
 a11<-ifelse(a11a>pmarge||a11a<0,ifelse(a11b<pmarge&&a11b>0,a11b,0),a11a)
 phi<-(a11-pmarge^2)/(pmarge*(1-pmarge))

    }",file="probit3.txt")


probit.b3<-function(data.sim){

  #Transform the data into long dataset
  n<-nrow(data.sim)
  nrat<-ncol(data.sim)
  
  y<-c(data.sim) 
  sub<-c(rep(seq(1,n),nrat))
  
  data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
  parameters<-c("ICC","phi")
  jags.inits<-list(beta0=0.0,sigma_S=4)
  obj <- jags.model("probit3.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model3.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  model3.sim<-summary(model3.s)
  
  return(list("ICC.LAT"=model3.sim[[2]][1,3],
              "LOW.LAT"=model3.sim[[2]][1,1],
              "HIGH.LAT"=model3.sim[[2]][1,5],
              "ICC"=model3.sim[[2]][2,3],
              "LOW"=model3.sim[[2]][2,1],
              "HIGH"=model3.sim[[2]][2,5]))   
  
}







# Probit - Bayesian (P.BT)
##########################
cat("
model {
  for (i in 1:n) {
    y[i] ~ dbern(p[i])
    p[i] <- phi(mu[i])
    mu[i]<-beta0+S[sub[i]]
  }

  beta0 ~ dnorm(0,1)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(0,tau_S)  
  }
 
 # Priors for variance components
 tau_S<-1/sigma2_S
 sigma2_S<-pow(sigma_S,2)
 #sigma_S~dunif(0,100)
 sigma_S~dt(0,1/6.25,3)T(0,)
  
 ICC<-sigma2_S/(sigma2_S+1)
 
 beta0_marge<-beta0/(sqrt(1+sigma2_S))
 pmarge<-phi(beta0_marge)
 
 c<-(1-(0.5-pmarge)^2)/2 
 w<-((3.14159-acos(ICC))/acos(ICC))^(1/c)
  
  b_<-(2*pmarge*w-2*pmarge+1)/(2*(w-1))
  sq_r<-sqrt(-4*pmarge^2*w+4*pmarge^2+4*pmarge*w-4*pmarge+1)/(2*(w-1))
  a11a<-b_+sq_r
  a11b<-b_-sq_r

 a11<-ifelse(a11a>pmarge||a11a<0,ifelse(a11b<pmarge&&a11b>0,a11b,0),a11a)
 phi<-(a11-pmarge^2)/(pmarge*(1-pmarge))

    }",file="probit5.txt")


probit.b5<-function(data.sim){

  #Transform the data into long dataset
  n<-nrow(data.sim)
  nrat<-ncol(data.sim)
  
  y<-c(data.sim) 
  sub<-c(rep(seq(1,n),nrat))
  
  data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
  parameters<-c("ICC","phi")
  jags.inits<-list(beta0=0.0,sigma_S=4)
  obj <- jags.model("probit5.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model3.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  model3.sim<-summary(model3.s)
  
  return(list("ICC.LAT"=model3.sim[[2]][1,3],
              "LOW.LAT"=model3.sim[[2]][1,1],
              "HIGH.LAT"=model3.sim[[2]][1,5],
              "ICC"=model3.sim[[2]][2,3],
              "LOW"=model3.sim[[2]][2,1],
              "HIGH"=model3.sim[[2]][2,5]))   

}


# Probit - Bayesian (P.BC)
##########################

cat("
model {
  for (i in 1:n) {
    y[i] ~ dbern(p[i])
    p[i] <- phi(mu[i])
    mu[i]<-beta0+S[sub[i]]
  }

  beta0 ~ dnorm(0,1)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(0,tau_S)  
  }
 
 # Priors for variance components
 tau_S<-1/sigma2_S
 sigma2_S<-pow(sigma_S,2)
 sigma_S~dt(0,1,1)T(0,)
  
 ICC<-sigma2_S/(sigma2_S+1)
 
 beta0_marge<-beta0/(sqrt(1+sigma2_S))
 pmarge<-phi(beta0_marge)
 
 c<-(1-(0.5-pmarge)^2)/2 
 w<-((3.14159-acos(ICC))/acos(ICC))^(1/c)
  
  b_<-(2*pmarge*w-2*pmarge+1)/(2*(w-1))
  sq_r<-sqrt(-4*pmarge^2*w+4*pmarge^2+4*pmarge*w-4*pmarge+1)/(2*(w-1))
  a11a<-b_+sq_r
  a11b<-b_-sq_r

 a11<-ifelse(a11a>pmarge||a11a<0,ifelse(a11b<pmarge&&a11b>0,a11b,0),a11a)
 phi<-(a11-pmarge^2)/(pmarge*(1-pmarge))

    }",file="probit6.txt")


probit.b6<-function(data.sim){
  check<-all(sweep(data.sim,MARGIN=1,STATS=data.sim[,1],FUN="=="))
  
  #Transform the data into long dataset
  n<-nrow(data.sim)
  nrat<-ncol(data.sim)
  
  y<-c(data.sim) 
  sub<-c(rep(seq(1,n),nrat))
  
  data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
  parameters<-c("ICC","phi")
  jags.inits<-list(beta0=0.0,sigma_S=4)
  obj <- jags.model("probit6.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model3.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  model3.sim<-summary(model3.s)
  
  return(list("ICC.LAT"=model3.sim[[2]][1,3],
              "LOW.LAT"=model3.sim[[2]][1,1],
              "HIGH.LAT"=model3.sim[[2]][1,5],
              "ICC"=model3.sim[[2]][2,3],
              "LOW"=model3.sim[[2]][2,1],
              "HIGH"=model3.sim[[2]][2,5]))   
  

}



# Probit - frequentist (P.F)
##################################


bootp_icc <- function(model, subject, B = 1000, seed = 2025) {
  set.seed(seed)
  sims <- simulate(model, B, seed)
  
  out <- vapply(seq_len(B), function(i) {
    mod <- suppressMessages(lme4::refit(model, newresp = sims[[i]]))
    
    # Check for singular fit and stop execution if detected
    if (isSingular(mod)) {
      warning("Singular fit detected in bootstrap iteration ", i, ". Returning NA.")
      return(c(NA, NA))  # Return NA values if singularity is encountered
    }
    
    prob_icc(mod, 'sub')
  }, numeric(2))  # Ensures numeric output
  
  # If all iterations resulted in NA, return NA
  if (all(is.na(out))) {
    return(NA)
  }
  
  return(apply(out, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}



prob_icc<-function(model,subject){
  sigma2<-lme4::VarCorr(model)[[subject]][1]
  ICC_lat<-sigma2/(sigma2+1)
  mu<-lme4::fixef(model)
  mu_marge<-mu/(sqrt(1+sigma2))
  p_marge<-pnorm(mu_marge)
  
  lower<-c(-mu_marge,-mu_marge)
  upper<-c(rep(Inf, 2))
  
  p11<-mvtnorm::pmvnorm(lower=lower,upper=upper,mean=c(0,0),corr= matrix(c(1,ICC_lat,ICC_lat,1),ncol=2))
  attributes(p11) <- NULL
  phi<-(p11-p_marge^2)/(p_marge*(1-p_marge))
  attributes(phi)<-NULL
  return(c(ICC_lat,phi))
  
}


probit.f<-function(data.sim,sim.num=500,seed.num=123){
  
  #Transform the data into long dataset
  n<-nrow(data.sim)
  nrat<-ncol(data.sim)
  
  y<-c(data.sim) 
  sub<-c(rep(seq(1,n),nrat))
  dat<-data.frame(cbind(y,sub))
  
  l_res<- try(glmer(y~1+(1|sub),data=dat,family='binomial'(link = "probit"),control = glmerControl(optimizer = "bobyqa"),
                    nAGQ = 15))
  
  if (any(inherits(l_res,"try-error")))
  {
    return(list("ICC.LAT"=NA,
                "LOW.LAT"=NA,
                "HIGH.LAT"=NA,
                "ICC"=NA,
                "LOW"=NA,
                "HIGH"=NA
    ))
  }
  
  if (!any(inherits(l_res,"try-error"))) {
    
    probICC<-prob_icc(l_res,'sub')
    
    b.ICC<-try(bootp_icc(l_res, subject, B=sim.num, seed=seed.num))
    
    if (any(inherits(b.ICC,"try-error")))
    {return(list("ICC.LAT"=probICC[1],
                 "LOW.LAT"=NA,
                 "HIGH.LAT"=NA,
                 "ICC"=probICC[2],
                 "LOW"=NA,
                 "HIGH"=NA
    ))}
    
    if (!any(inherits(b.ICC,"try-error"))){
      return(list("ICC.LAT"=probICC[1],
                  "LOW.LAT"=b.ICC[1,1],
                  "HIGH.LAT"=b.ICC[3,1],
                  "ICC"=probICC[2],
                  "LOW"=b.ICC[1,2],
                  "HIGH"=b.ICC[3,2]
      ))}
    
  }
  
}



# Normal approximation - frequentist (N.F)
#########################################

normal.app.f<-function(data.sim){
  
  check<-all(sweep(data.sim,MARGIN=1,STATS=data.sim[,1],FUN="=="))
  
  if (check==TRUE){
    return(list("ICC"=1,
                "LOW"=1,
                "HIGH"=1,
                "LOW.Z"=1,
                "HIGH.Z"=1
    ))
  }
  if (check==FALSE){
    #Normal approximation
    n<-nrow(data.sim)
    nrat<-ncol(data.sim)
    
    y<-c(data.sim) 
    sub<-c(rep(seq(1,n),nrat))
    DAT<-data.frame(cbind(y,sub))
    
    fit.SS3 <- fitVCA(y~sub,Data=DAT)
    terms<-fit.SS3$aov.tab
    
    ICC<-(terms[2,4])/terms[1,4]
    
    inf.SS3 <-vcovVC(fit.SS3, method="scm")
    
    DIR<-matrix(c((1-ICC)/terms[1,4],-ICC/terms[1,4]),nrow=1)
    SE_inter<-sqrt(DIR%*%inf.SS3%*%t(DIR))
    z.res<-Z.transform(ICC,SE_inter)
    
    #Wald confidence interval
    return(list("ICC"=ICC,
                "LOW"=max(0,ICC-1.96*SE_inter),
                "HIGH"=min(ICC+1.96*SE_inter,1),
                "LOW.Z"=unlist(unname(z.res[1])),
                "HIGH.Z"=unlist(unname(z.res[2]))
    ))
    
  }
}


# Normal approximation - Bayesian (N.B)
########################################

cat("
  model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0+S[sub[i]]
  }

  beta0 ~ dnorm(0,0.001)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(0, tau_S)  
  }
 
  
  # Priors for variance components
  tau_S<-1/sigma2_S
  tau<-1/sigma2
  
 sigma2_S<-pow(sigma_S,2)
 sigma2<-pow(sigma,2)
 
 sigma_S~dunif(0,0.5)
 sigma~dunif(0,0.5)

 ICC<-(sigma2_S)/(sigma2_S+sigma2)

    }",file="norm_approx1.txt")


normal.app.b1<-function(data.sim){
  #Jags-normal approximation
  check<-all(sweep(data.sim,MARGIN=1,STATS=data.sim[,1],FUN="=="))
  
  if (check==TRUE){
    return(list("ICC"=1,
                "LOW"=1,
                "HIGH"=1
    ))
  }
  
  
  
  if (check==FALSE){
    
    #Transform the data into long dataset
    n<-nrow(data.sim)
    nrat<-ncol(data.sim)
    
    y<-c(data.sim) 
    sub<-c(rep(seq(1,n),nrat))
    
    data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
    parameters<-c("ICC","beta0","sigma2","sigma2_S")
    jags.inits<-list(beta0=0.0,sigma_S=0.10,sigma=0.05)
    obj <- jags.model("norm_approx1.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
    jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
    model3.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
    model3.sim<-summary(model3.s)
    
    return(list("ICC"=model3.sim[[2]][1,3],
                "LOW"=model3.sim[[2]][1,1],
                "HIGH"=model3.sim[[2]][1,5]
    ))
  }
}



MC_Simulation <- function(nsubject, nrater, prop,ICC.target, nsim=500, seed.s=1234){
  cl<- makeCluster(50, outfile="datsim.txt")
  cat("ICC.target=",ICC.target, "prop=",prop,"nsubject=", nsubject, " nrater=",nrater, " ****Start*****\n")
  st <- Sys.time()
  set.seed(seed.s)
 
   library(lme4)
  library(mvtnorm)
  #g<-lapply(1:nsim,function(i){simu(prop=prop,nr=nrater,n=nsubject,ICC.sim=ICC.target)})
  
  g<-lapply(1:nsim,function(i){simu(prop=prop,n.rat=nrater,n.sub=nsubject,ICC.sim=ICC.target)})
 
  on.exit(stopCluster(cl), add = TRUE)
  
  
  shared_probit1<-"./probit1.txt"
  shared_probit3<-"./probit3.txt"
  shared_probit5<-"./probit5.txt"
  shared_probit6<-"./probit6.txt"
  shared_probit9<-"./norm_approx1.txt"

  clusterExport(cl, list("ICC.target",
                         "prop",
                         "nsubject",
                         "nrater",
                         "nsim",
                         "g",
                         "Z.transform",
                         "simu",
                         "tetra_app",
                         "phi_app",
                         "probit.b1",
                         "probit.b3",
                         "probit.b5",
                         "probit.b6",
                         "normal.app.b1",
                         "kappa.bayes.bin",
                         "kappaw",
                         "shared_probit1",
                         "shared_probit3",
                         "shared_probit5",
                         "shared_probit6",
                         "shared_probit9",
                         "probit.f",
                         "bootp_icc",
                         "prob_icc",
                         "normal.app.f"
  ),
  envir=environment())
  clusterEvalQ(cl,{
   
    library(rjags)
    library(coda)
    library(mvtnorm)
    library(utils)
    library(stats)
    library(Matrix)
    library(MCMCglmm)
    library(rstan)
    library(lme4)
    library(VCA)
    suppressPackageStartupMessages(library(rjags))
    suppressPackageStartupMessages(library(rjags))
    
  })
  
  
  cis<-
    parLapply(cl, 1:nsim, function(x){
      list( "p.b1"= probit.b1(data.sim =2-g[[x]]$data.sim),
            "p.b3"= probit.b3(data.sim =2-g[[x]]$data.sim),
            "p.b5"= probit.b5(data.sim =2-g[[x]]$data.sim),
            "p.b6"= probit.b6(data.sim =2-g[[x]]$data.sim),
            "k.b"= kappa.bayes.bin(data=g[[x]]$data.sim),
            "k.f"= kappaw(data=g[[x]]$data.sim,"binary",3),
            "n.f"= normal.app.f(data.sim =2-g[[x]]$data.sim),
            "n.b1"= normal.app.b1(data.sim =2-g[[x]]$data.sim),
            "p.f"= probit.f(data.sim=2-g[[x]]$data.sim,sim.num=500,seed.num=x),
            "p.simu"=mean(2-g[[x]]$data.sim))
    })
  
  
  
  cat("ICC.target=",ICC.target, "prop=",prop,"nsubject=", nsubject, " nrater=",nrater, " completed in", Sys.time()-st,"\n")
  
  
  
  cov.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("k.f","k.b","p.f","p.b1","p.b3","p.b5","p.b6","n.f","n.b1"), 
           function(y) ifelse(is.na(x[[y]][['LOW']])|is.na(x[[y]][['HIGH']]),0,ifelse(x[[y]][['LOW']]>=ICC.target|x[[y]][['HIGH']]<=ICC.target,0,1)))))
  
  cov.z.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("n.f","k.f"), 
           function(y) ifelse(is.na(x[[y]][['LOW.Z']])|is.na(x[[y]][['HIGH.Z']]),0,ifelse(x[[y]][['LOW.Z']]>=ICC.target|x[[y]][['HIGH.Z']]<=ICC.target,0,1)))))
  
  rho.check<-tetra_app(phi=ICC.target,pmarge=prop)
  
  
  cov.l.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("p.b1","p.b3","p.b5","p.b6"), 
           function(y) ifelse(is.na(x[[y]][['LOW.LAT']])|is.na(x[[y]][['HIGH.LAT']]),0,ifelse(x[[y]][['LOW.LAT']]>=rho.check|x[[y]][['HIGH.LAT']]<=rho.check,0,1)))))
  
  l.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("k.f","k.b","p.f","p.b1","p.b3","p.b5","p.b6","n.f","n.b1"), 
           function(y) x[[y]][['HIGH']]-x[[y]][['LOW']])),na.rm=TRUE)
  
  l.z.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("n.f","k.f"), 
           function(y) x[[y]][['HIGH.Z']]-x[[y]][['LOW.Z']])),na.rm=TRUE)
  
  
  l.l.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("p.b1","p.b3","p.b5","p.b6"), 
           function(y) x[[y]][['HIGH.LAT']]-x[[y]][['LOW.LAT']])),na.rm=TRUE)
  
  m.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("k.f","k.b","p.f","p.b1","p.b3","p.b5","p.b6","n.f","n.b1"), 
           function(y) x[[y]][['ICC']])),na.rm=TRUE)
  
  m.l.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("p.b1","p.b3","p.b5","p.b6","p.b7","p.b8"), 
           function(y) x[[y]][['ICC.LAT']])),na.rm=TRUE)
  
  bias.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("k.f","k.b","p.f","p.b1","p.b3","p.b5","p.b6","n.f","n.b1"), 
           function(y) abs(x[[y]][['ICC']]-ICC.target))),na.rm=TRUE)
  
  bias.l.icc<- rowMeans(sapply(cis, function(x) 
    sapply(c("p.b1","p.b3","p.b5","p.b6"), 
           function(y) abs(x[[y]][['ICC.LAT']]-rho.check))),na.rm=TRUE)
  
  m.p.simu<- mean(sapply(cis, function(x) 
    sapply(c("p.simu"), 
           function(y) x[[y]])),na.rm=TRUE)
  
  
  return(list("ICC.target"=ICC.target,
              "prop" = prop,
              "nsubject" = nsubject,
              "nrater" = nrater,
              "nsim"= nsim,
              "seed.s"= seed.s,
              "cov.icc" =cov.icc,
              "l.icc" =l.icc,
              "m.icc" =m.icc,
              "bias.icc" =bias.icc,
              "cov.z.icc" =cov.z.icc,
              "l.z.icc" =l.z.icc,
              "cov.l.icc" =cov.l.icc,
              "l.l.icc" =l.l.icc,
              "m.l.icc" =m.l.icc,
              "bias.l.icc" =bias.l.icc,
              "m.p.sim" =m.p.simu,
              "CI"  = cis))
}



gd <- expand.grid("ICC.target"=c(0.7,0.8,0.9),
                  "nsubject"  =c(100,80,60,40,30,20),
                  "nrater"  =c(2,3,5,8),
                  "prop"=0.5)

gd[["fname"]] = sapply(1:nrow(gd), function(x){paste0("Sim_",paste0(names(gd[x,]),"_",gd[x,],collapse = "_"), ".RDS")})


sts <- lapply(1:nrow(gd), function(x){
  MC <- MC_Simulation(nsubject        = gd[x,2],
                      nrater        = gd[x,3],
                      ICC.target      = gd[x,1],
                      prop        = gd[x,4],
                      nsim     = 500,
                      seed.s     = x)
  saveRDS(MC, gd[x,5])
})
