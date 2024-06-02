#######################################################
########### Simulations for h1 PDF
#######################################################

## h1 PDF
h1<-function(x, a, b, g){
  ifelse((g*(x^(b)))>1,(a*b/x)*((log(g*x^b))^(-a-1))*exp(-(log(g*x^b))^(-a)), 0)
}


## h1 log-likelihood function
ll1 <- function(p,x){
  a<-p[1]
  b<-p[2]
  g<-p[3]
  
  
  if(min(x)>(g^(-1/b))){
    ll_out<- length(x)*(log(a)+log(b)) - sum(log(x)) - (a+1)*sum(log(log(g)+b*log(x))) - sum((log(g)+b*log(x))^(-a))
  }
  if(min(x)<=(g^(-1/b))){
    ll_out<- 0
  }
  return(ll_out)
}

#Estimator for h1 PDF
estimator_h1<-function(x, epslon=0.01, beta_error=0.05, N=5){
  
  m=(N*N*N)
  
  n=ceiling(-(log(beta_error)*m/epslon))
  
  ai=runif(n+1,min=0, max=N)
  bi=runif(n+1,min=0, max=N)
  
  valor_max1=ifelse((min(x)^(-bi[1]))>N, (min(x)^(-bi[1]))+N, N)
  gi<- runif(1,min=(min(x)^(-bi[1])), max=valor_max1)
  for(i in 2:(n+1)){
    valor_maxi=ifelse((min(x)^(-bi[i]))>N, (min(x)^(-bi[i]))+N, N )
    gi[i]=runif(1,min=(min(x)^(-bi[i])), max=valor_maxi)
  }

  
  p=c(ai[1], bi[1], gi[1])
  ll_max=ll1(p, x)
  for(k in 1:n){
    p_test <- c(ai[k+1], bi[k+1], gi[k+1])
    ll_p_test<-ll1(p_test, x)
    if(ll_p_test>ll_max){
      p<-p_test
      ll_max<-ll_p_test
    }
    
  }
  
  return(c(p, ll_max))
}

#extreme value H-function
H<-function(a1,a2, a3, a4, a5){
  integrand <- function(y){exp(-a1*y - (a2*(y^(a3)) + a4)^(a5))}
  integralv<-integrate(integrand, 0, Inf)$value
  return (integralv) 
}


#R estimator for h1 marginals
fit_R_h1=function(x_est, y_est){
  
  
  A1=y_est[1]
  B1=y_est[2] 
  G1=y_est[3]
  
  A2=x_est[1]
  B2=x_est[2] 
  G2=x_est[3]
  
  out<-ifelse(G1^(-1/B1)>=G2^(-1/B2), 
              H(a1=1,a2=(B2/B1), a3=-(1/A1), a4=(log(G2)-(B2/B1)*log(G1)), a5=-A2),
              1-H(a1=1,a2=(B1/B2), a3=-(1/A2), a4=(log(G1)-(B1/B2)*log(G2)), a5=-A1))
  
  return(out)
  
}



#generate random samples from h1
rh1=function(n, a, b,g){
  u <- runif(n)
  # inverse function
  sample= (b^(-1/g))*exp((1/g)*((-log(u))^(-1/a)))  
  return(sample)
}




#parameters matrice
m_par = data.frame(a2=rep(1,10),
                   b2=c(rep(0.5, 3),rep(0.3,3),  rep(1,4)), 
                   g2=c(rep(0.3,3), rep(0.55,3), rep(1,4)),
                   a1=c(rep(0.7,2),1,  .7,.7,1,  .90,.95,.97,1),
                   b1=c(0.5, .75,.75, .55,.55,0.55,  rep(0.5,3), 0.5),
                   g1=c(0.3,0.3, 0.3, rep(.9,3),  .5,.5,.5,.5))
m_par$R=NA
for(k in 1:nrow(m_par)){
  m_par$R[k]=fit_R_h1(c(m_par[k,1], m_par[k,2], m_par[k,3]),
                      c(m_par[k,4], m_par[k,5], m_par[k,6]))
}


m_par$R_mc_mean=NA
m_par$Bias=NA
m_par$RMSE=NA


N_mc=1000
n_sample=30

for(j in 1:nrow(m_par)){

  R_mc=NULL
  for(MC in 1:N_mc){
    y_sim=rh1(n=n_sample,a=m_par$a1[j],b=m_par$b1[j],g=m_par$g1[j])
    x_sim=rh1(n=n_sample,a=m_par$a2[j],b=m_par$b2[j],g=m_par$g2[j])
    
    x_sim_fit=estimator_h1(x=x_sim, N=5)
    y_sim_fit=estimator_h1(x=y_sim, N=5)

    R_mc[MC]=fit_R_h1(x_sim_fit,y_sim_fit)
  }
  m_par$R_mc_mean[j]=mean(R_mc)
  m_par$Bias[j]= m_par$R_mc_mean[j]- m_par$R[j]
  m_par$RMSE[j]= sqrt( (mean(m_par$R_mc_mean[j] -m_par$R[j]))^2 + var(R_mc))
}


#######################################################
########### Simulations for h2 PDF
#######################################################

## h2 PDF
h2<-function(x, a, b, g){
  y<- g*x^(b)
  ifelse(x>0 & x<g^(-1/b),(a*b/x)*((-log(g*(x^b)))^(a-1))*exp(-(-log(g*(x^b)))^(a)), 0)
}


## log-likelihood for h2 
ll2<-function(p, x){
  a<-p[1]
  b<-p[2]
  g<-p[3]
  if(max(x)<(g^(-1/b))){
    ll_out<- length(x)*(log(a)+log(b)) - sum(log(x)) + (a-1)*sum(log(-log(g)-b*log(x))) - sum((-log(g)-b*log(x))^a)
  }
  if(max(x)>=(g^(-1/b))){
    ll_out<- -Inf
  }
  return(ll_out)
}

## Estimator for h2 PDF
estimator_h2<-function(x, epslon=0.01, beta_error=0.05, N=10){
  
  m=(N*N*N)
  
  n=ceiling(-(log(beta_error)*m/epslon))
  
  ai=runif(n+1,min=0, max=N)
  bi=runif(n+1,min=0, max=N)
  
  gi<- runif(1,min=0, max=min(N,(min(x)^(-bi[1]))))
  for(i in 2:(n+1)){
    gi[i]=runif(1,min=0, max=min(N,min(x)^(-bi[i])))
  }
  
  p=c(ai[1], bi[1], gi[1])
  
  ll_max=ll2(p, x)
  for(k in 1:n){
    p_test <- c(ai[k+1], bi[k+1], gi[k+1])
    ll_p_test<-ll2(p_test, x)
    if(ll_p_test>ll_max){
      p<-p_test
      ll_max<-ll_p_test
    }
    
  }
  
  return(c(p, ll_max))
}

#Extreme value H-function
H<-function(a1,a2, a3, a4, a5){
  integrand <- function(y){exp(-a1*y - (a2*(y^(a3)) + a4)^(a5))}
  integralv<-integrate(integrand, 0, Inf)$value
  return (integralv) 
}


#Estimator of R for h2 marginals
fit_R=function(x_est, y_est){
  
  
  A1=y_est[1]
  B1=y_est[2] 
  G1=y_est[3]
  
  A2=x_est[1]
  B2=x_est[2] 
  G2=x_est[3]
  
  out<-ifelse(G1^(-1/B1)<=G2^(-1/B2), 
              H(a1=1,a2=(B2/B1), a3=(1/A1), a4=(-log(G2)+(B2/B1)*log(G1)), a5=A2),
              1-H(a1=1,a2=(B1/B2), a3=(1/A2), a4=(-log(G1)+(B1/B2)*log(G2)), a5=A1))
  
  return(out)
  
}



#generate random samples from h2
rh2=function(n, a, b,g){
  u <- runif(n)
  sample= (b^(-1/g))*exp((-1/g)*((-log(u))^(1/a)))  
  return(sample)
}

#parameter matrice
m_par = data.frame(a2=rep(1,10),
                   b2=c(rep(0.5, 3),rep(0.3,3),  rep(1,4)), 
                   g2=c(rep(0.3,3), rep(0.55,3), rep(1,4)),
                   a1=c(rep(0.7,2),1,  .7,.7,1,  .90,.95,.97,1),
                   b1=c(0.5, .75,.75, .55,.55,0.55,  rep(0.5,3), 0.5),
                   g1=c(0.3,0.3, 0.3, rep(.9,3),  .5,.5,.5,.5))
m_par$R=NA
for(k in 1:nrow(m_par)){
  m_par$R[k]=fit_R(c(m_par[k,1], m_par[k,2], m_par[k,3]),
                   c(m_par[k,4], m_par[k,5], m_par[k,6]))
}


m_par$R_mc_mean=NA
m_par$Bias=NA
m_par$RMSE=NA


N_mc=1000
n_sample=200

for(j in 1:nrow(m_par)){
  
  R_mc=NULL
  for(MC in 1:N_mc){
    y_sim=rh2(n=n_sample,a=m_par$a1[j],b=m_par$b1[j],g=m_par$g1[j])
    x_sim=rh2(n=n_sample,a=m_par$a2[j],b=m_par$b2[j],g=m_par$g2[j])
    
    x_sim_fit=estimator_h2(x=x_sim, N=10)
    y_sim_fit=estimator_h2(x=y_sim, N=10)

    R_mc[MC]=fit_R(x_sim_fit,y_sim_fit)
  }
  m_par$R_mc_mean[j]=mean(R_mc)

  m_par$Bias[j]= m_par$R_mc_mean[j]- m_par$R[j]
  m_par$RMSE[j]= sqrt( (mean(m_par$R_mc_mean[j] -m_par$R[j]))^2 + var(R_mc))
}

#######################################################
########### Simulations for h5 PDF
#######################################################

## h5 PDF 
h5<-function(x, b, g){
  ifelse(x>0 ,(b/g)*(x^(-b-1))*exp(-(g*x^(b))^(-1)), 0)
}


## log-likelihood for h5
ll5<-function(p, x){
  b<-p[1]
  g<-p[2]
  if(min(x)>0){
    ll_out<- length(x)*(log(b)-log(g)) -(b+1)*sum(log(x)) - sum((g*(x^b))^(-1))
  }
  if(min(x)<=0){
    ll_out<- 0
  }
  return(ll_out)
}

#Estimator for h5 PDF
estimator_h5<-function(x,epslon=0.01, beta_error=0.05, N=10){
  
  m=(N*N)
  n=ceiling(-(log(beta_error)*m/epslon))
  
  bi=runif(n+1,min=0, max=N)
  gi<- runif(n+1,min=0, max=N)
  
  
  p<-c(bi[1], gi[1])
  ll_max=ll5(p, x)
  for(k in 1:n){
    p_test <- c(bi[k+1], gi[k+1])
    ll_p_test<-ll5(p_test, x)
    if(ll_p_test>ll_max){
      p<-p_test
      ll_max<-ll_p_test
    }
    
  }

  return(c(p, ll_max))
}

#extreme value H-function
H<-function(a1,a2, a3, a4, a5){
  integrand <- function(y){exp(-a1*y - (a2*(y^(a3)) + a4)^(a5))}
  integralv<-integrate(integrand, 0, Inf)$value
  return (integralv) 
  }


#Estimator of R for h5 marginals
fit_R_h5=function(x_est, y_est){
  
  
  B1=y_est[1] 
  G1=y_est[2]
  
  B2=x_est[1] 
  G2=x_est[2]
  
  out<-(1/G1)* H(a1=1/G1,a2=1/G2, a3=B2/B1, a4=0, a5=1)
  
  return(out)
  
}

#generate random samples from h5 PDF
rh5=function(n, b,g){
  u <- runif(n)
  sample= (-b*log(u))^(-1/g)
  return(sample)
}

#parameter matrices
m_par = data.frame(
  b2=c(rep(0.3, 3),rep(0.5,3),  rep(1,4)), 
  g2=c(rep(0.5,3), rep(0.7,3), rep(1,4)),
  b1=c(rep(c(0.2,0.6,0.9),3), 1),
  g1=c(rep(c(0.3,0.5,1),3), 1))
m_par$R=NA
for(k in 1:nrow(m_par)){
  m_par$R[k]=fit_R_h5(c(m_par[k,1], m_par[k,2]),
                      c(m_par[k,3], m_par[k,4]))
}


m_par$R_mc_mean=NA
m_par$Bias=NA
m_par$RMSE=NA

N_mc=1000
n_sample=30

for(j in 1:nrow(m_par)){
  
  R_mc=NULL
  for(MC in 1:N_mc){
    y_sim=rh5(n=n_sample,b=m_par$b1[j],g=m_par$g1[j])
    x_sim=rh5(n=n_sample,b=m_par$b2[j],g=m_par$g2[j])
    
    x_sim_fit=estimator_h5(x=x_sim, N=10)
    y_sim_fit=estimator_h5(x=y_sim, N=10)

    R_mc[MC]=fit_R_h5(x_sim_fit,y_sim_fit)
  }
  m_par$R_mc_mean[j]=mean(R_mc)
  m_par$Bias[j]= m_par$R_mc_mean[j]- m_par$R[j]
  m_par$RMSE[j]= sqrt( (mean(m_par$R_mc_mean[j] -m_par$R[j]))^2 + var(R_mc))
}

