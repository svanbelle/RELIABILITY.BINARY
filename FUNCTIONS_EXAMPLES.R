
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



# Return median and CI based on DM
#######################################

kappa.bayes.bin<-function(data){
  
  # check<-all(sweep(data,MARGIN=1,STATS=data[,1],FUN="=="))
  # 
  # if (check==TRUE){
  #   return(list("ICC"=1,
  #               "LOW"=1,
  #               "HIGH"=1
  #   ))
  # }
  # if (check==FALSE){
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
  # }
}

# Return kappa coefficient with delta CI
#######################################

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
  a11<-a22<-phi*pmarge*(1-pmarge)+pmarge^2
  a12<-a21<-pmarge-a11
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

#Probit - Bayesian (half-t prior) - repeatability study
##################################
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


# Probit - Bayesian (half-cauchy 1) - reproducibility study
##################################
cat("
model {
  for (i in 1:n) {
    y[i] ~ dbern(p[i])
    p[i] <- phi(mu[i])
    mu[i]<-S[sub[i]]+R[met[i]]
  }

  beta0 ~ dnorm(0,1)  
 
  # Priors for random effects (standard deviations and variances)
  for (j in 1:ns) {
    S[j] ~ dnorm(beta0,tau_S)  
  }
  
    for (j in 1:nr) {
    R[j] ~ dnorm(0,tau_R)  
  }
 
 # Priors for variance components
 tau_S<-1/sigma2_S
 sigma2_S<-pow(sigma_S,2)
 sigma_S~dt(0,1,1)T(0,)
  
  
  tau_R<-1/sigma2_R
 sigma2_R<-pow(sigma_R,2)
 sigma_R~dt(0,1,1)T(0,)
  
   
 ICC<-sigma2_S/(sigma2_S+sigma2_R+1)
 beta0_marge<-beta0/(sqrt(1+sigma2_S+sigma2_R))
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



# Probit - frequentist - repeatability 
##################################


bootp_icc1<- function(model, subject, B = 1000, seed = 2025) {
  set.seed(seed)
  sims <- simulate(model, B, seed)
  
  out <- vapply(seq_len(B), function(i) {
    mod <- suppressMessages(lme4::refit(model, newresp = sims[[i]]))
    
    # Check for singular fit and stop execution if detected
    if (isSingular(mod)) {
      warning("Singular fit detected in bootstrap iteration ", i, ". Returning NA.")
      return(c(NA, NA))  # Return NA values if singularity is encountered
    }
    
    prob_icc1(mod, 'sub')
  }, numeric(2))  # Ensures numeric output
  
  # If all iterations resulted in NA, return NA
  if (all(is.na(out))) {
    return(NA)
  }
  
  return(apply(out, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}


prob_icc1<-function(model,subject){
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


# Probit - frequentist - reproducibility 
##################################

bootp_icc2 <- function(model, subject, B = 1000, seed = 2025) {
  set.seed(seed)
  sims <- simulate(model, B, seed)
  
  out <- vapply(seq_len(B), function(i) {
    mod <- suppressMessages(lme4::refit(model, newresp = sims[[i]]))
    
    # Check for singular fit and stop execution if detected
    if (isSingular(mod)) {
      warning("Singular fit detected in bootstrap iteration ", i, ". Returning NA.")
      return(c(NA, NA))  # Return NA values if singularity is encountered
    }
    
    prob_icc2(mod, 'sub','met')
  }, numeric(2))  # Ensures numeric output
  
  # If all iterations resulted in NA, return NA
  if (all(is.na(out))) {
    return(NA)
  }
  
  return(apply(out, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}



prob_icc2<-function(model,subject,rater){
  sigma2_S<-lme4::VarCorr(model)[[subject]][1]
  sigma2_R<-lme4::VarCorr(model)[[rater]][1]
  
  ICC_lat<-sigma2_S/(sigma2_S+sigma2_R+1)
  mu<-lme4::fixef(model)
  mu_marge<-mu/(sqrt(1+sigma2_S+sigma2_R))
  p_marge<-pnorm(mu_marge)
  
  lower<-c(-mu_marge,-mu_marge)
  upper<-c(rep(Inf, 2))
  
  p11<-mvtnorm::pmvnorm(lower=lower,upper=upper,mean=c(0,0),corr= matrix(c(1,ICC_lat,ICC_lat,1),ncol=2))
  attributes(p11) <- NULL
  phi<-(p11-p_marge^2)/(p_marge*(1-p_marge))
  phi_app(tetra=ICC_lat,pmarge=p_marge)
  attributes(phi)<-NULL
  return(c(ICC_lat,phi))
  
}

