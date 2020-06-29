
Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 1:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      #ifelse(A==0,NA,A^(1/(1-q)))
      A^(1/(1-q))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

Diversity_Tsallis <- function(x,q){
  qD = Diversity_profile(x, q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}

Diversity_Tsallis_MLE <- function(x,q){
  qD = Diversity_profile_MLE(x,q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}

bootstrap_forq = function(data,B,q,conf,FUNNAME){
  data <- data[data!=0]
  n <- sum(data)
  f1 = sum(data==1); f2 = sum(data==2)
  f0 = ceiling(ifelse( f2>0, (n-1)*f1^2/n/2/f2, (n-1)*f1*(f1-1)/2/n ))
  C_hat = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
  lamda_hat = (1-C_hat)/sum((data/n)*(1-data/n)^n) 
  pi_hat = (data/n)*(1-lamda_hat*((1-data/n)^n)) 
  p_hat = c( pi_hat, rep( (1-C_hat)/f0, f0 )) 
  random = rmultinom( B, n, p_hat ) 
  #Bt_estimate <- sapply(c(1:B),function(i) FUNNAME(random[,i],q))
  Bt_estimate <- apply(random,MARGIN = 2,function(i) FUNNAME(i,q))
  estimate <- FUNNAME(data,q)
  #Interval_mean = apply( Bt_estimate, 1, mean)
  Interval_mean = rowMeans(Bt_estimate)
  Interval_sd = apply(Bt_estimate, 1, sd)
  Interval_quantileL = apply( Bt_estimate, 1, quantile, p=(1-conf)/2)
  Interval_quantileU = apply( Bt_estimate, 1, quantile, p=1-(1-conf)/2)
  Upper_bound = estimate+Interval_quantileU-Interval_mean 
  Lower_bound = estimate+Interval_quantileL-Interval_mean
  result <- cbind("estimate"=estimate,"sd"=Interval_sd,"LCL"=Lower_bound,"UCL"=Upper_bound)
  result 
}

MakeTable_Proposeprofile_nose = function(data,q){
  Diversity = Diversity_profile(data,q)
  #Entropy = Diversity_Tsallis(Diversity,q)
  output = data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity)
  output[,c(1,3)] = round(output[,c(1,3)],9)
  
  return(output)
}

MakeTable_Proposeprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq(data, B, q, conf, Diversity_profile)
  #Entropy = bootstrap_forq(data, B, q, conf, Diversity_Tsallis)
  output = data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4])
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)
  
  return(output)
}

MakeTable_Empericalprofile_nose = function(data,q){
  Diversity = Diversity_profile_MLE(data,q)
  #Entropy = Diversity_Tsallis_MLE(Diversity,q)
  output = data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity)
  output[,c(1,3)] = round(output[,c(1,3)],9)
  
  return(output)
}

MakeTable_Empericalprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq( data, B,q,conf,Diversity_profile_MLE)
  #Entropy = bootstrap_forq( data, B,q,conf,Diversity_Tsallis_MLE)
  # tmp <- Diversity_Tsallis(Diversity[,1],q)
  # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4])
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)
  
  return(output)
}

plot_diversity_profile_nose_yhc=function(output){
  ggplot(output, aes(x=Order.q, y=Diversity, colour=Community,lty = method))+
  geom_line(size=1.2) +    
  labs(x="Order q", y="Diversity")+
  theme(text=element_text(size=18))
}

plot_diversity_profile_yhc=function(output){
  ggplot(output, aes(x=Order.q, y=Diversity, colour=Community, lty = method))+
  labs(x="Order q", y="Diversity")+
  theme(text=element_text(size=18))+
  geom_ribbon(data = output %>% filter(method=="Asymptotic"),
              aes(ymin=LCL, ymax=UCL, fill=Community, colour=NULL), alpha=0.4, linetype=0)+
  geom_line(size=1.1)
  
}

#=======incidence profile==========#
#cpp function in Diversity_profile.inc
cppFunction(
  "NumericVector qDFUN(NumericVector q,NumericVector Xi,const int n){
    const int length = q.size();
    const int Sobs = Xi.size();
    NumericVector Q(length);
    NumericVector delta(n);
    NumericVector temp(Sobs);
    for(int k=0;k<=(n-1);k++){
    for(int i = 0;i<Sobs;i++){
    temp[i] = (Xi[i]/n)*exp(Rf_lchoose(n-Xi[i],k)-Rf_lchoose(n-1,k));
    }
    delta[k] = sum(temp);
    }
    
    for(int i=0;i<length;i++){
    float temp = 0;
    for(int k=0;k<=(n-1);k++){
    temp = temp + (Rf_choose(q[i]-1,k)*pow(-1,k)*delta[k]);
    }
    Q[i] = temp;
    }
    return Q;
  }")
#cpp function in Diversity_profile_MLE.inc
cppFunction(
  "NumericVector qD_MLE(NumericVector q,NumericVector ai){
    const int length = q.size();
    const int S = ai.size();
    NumericVector Q(length);
    NumericVector temp(S);
    for(int j = 0; j<length;j++){
    for(int i = 0 ; i<S;i++){
    temp[i] = pow(ai[i],q[j]);
    }
    Q[j] = pow(sum(temp),1/(1-q[j]));
    }
    return Q;
}")
AA.inc <- function(data){
  T = data[1]
  U <- sum(data[-1])
  data = data[-1]
  Yi = data[data!=0]
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  if(Q2>0 & Q1>0){
    A <- 2*Q2/((T-1)*Q1+2*Q2)
  }
  else if(Q2==0 & Q1>1){
    A <- 2/((T-1)*(Q1-1)+2)
  }
  else{
    A <- 1
  }
  return(A)
}

Diversity_profile.inc <- function(data,q){
  T = data[1]
  Yi = data[-1]
  Yi <- Yi[Yi!=0]
  U <- sum(Yi)
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  Sobs <- length(Yi)
  A <- AA.inc(data)
  Q0hat <- ifelse(Q2 == 0, (T - 1) / T * Q1 * (Q1 - 1) / 2, (T - 1) / T * Q1 ^ 2/ 2 / Q2)
  B <- sapply(q,function(q) ifelse(A==1,0,(Q1/T)*(1-A)^(-T+1)*(A^(q-1)-sum(sapply(c(0:(T-1)),function(r) choose(q-1,r)*(A-1)^r)))))
  qD <- (U/T)^(q/(q-1))*(qDFUN(q,Yi,T) + B)^(1/(1-q))
  qD[which(q==0)] = Sobs+Q0hat
  yi <- Yi[Yi>=1 & Yi<=(T-1)]
  delta <- function(i){
    (yi[i]/T)*sum(1/c(yi[i]:(T-1)))
  }
  if(sum(q %in% 1)>0){
    C_ <- ifelse(A==1,0,(Q1/T)*(1-A)^(-T+1)*(-log(A)-sum(sapply(c(1:(T-1)),function(r) (1-A)^r/r))))
    qD[which(q==1)] <- exp((T/U)*( sum(sapply(c(1:length(yi)),function(i) delta(i))) + C_)+log(U/T))
  }
  return(qD)
}

Diversity_profile_MLE.inc <- function(data,q){
  Yi = data[-1]
  U = sum(Yi)
  Yi <- Yi[Yi!=0]
  ai <- Yi/U
  qD = qD_MLE(q,ai)
  qD[which(q==1)] <- exp(-sum(ai*log(ai)))
  return(qD)
}

Diversity_Tsallis.inc <- function(qD,q){
  #qD = Diversity_profile.inc(data,q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}

Diversity_Tsallis_MLE.inc <- function(qD,q){
  #qD = Diversity_profile_MLE.inc(data,q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}

bootstrap_forq.inc = function(data,B,q,conf,FUNNAME){
  T <- data[1]
  Yi <- data[-1]
  Yi <- Yi[Yi>0]
  Sobs <- sum(Yi > 0)   
  Q1 <- sum(Yi == 1) 	
  Q2 <- sum(Yi == 2) 
  Q0.hat <- ifelse(Q2 == 0, (T - 1) / T * Q1 * (Q1 - 1) / 2, (T - 1) / T * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
  A <- ifelse(Q1>0, T*Q0.hat/(T*Q0.hat+Q1), 1)
  a <- Q1/T*A
  b <- sum(Yi / T * (1 - Yi / T) ^ T)
  w <- a / b  		
  Prob.hat <- Yi / T * (1 - w * (1 - Yi / T) ^ T)					
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))  	
  p_hat =c(Prob.hat, Prob.hat.Unse)
  random = t(sapply(1:length(p_hat), function(i){rbinom(B,T,p_hat[i])}))
  random = rbind(rep(T,B),random)
  Bt_estimate <- apply(random,MARGIN = 2,function(i) FUNNAME(i,q))
  estimate <- FUNNAME(data,q)
  Interval_mean = rowMeans(Bt_estimate)
  Interval_sd = apply(Bt_estimate, 1, sd)
  Interval_quantileL = apply( Bt_estimate, 1, quantile, p=(1-conf)/2)
  Interval_quantileU = apply( Bt_estimate, 1, quantile, p=1-(1-conf)/2)
  Upper_bound = estimate+Interval_quantileU-Interval_mean 
  Lower_bound = estimate+Interval_quantileL-Interval_mean
  result <- cbind("estimate"=estimate,"sd"=Interval_sd,"LCL"=Lower_bound,"UCL"=Upper_bound)
  result 
}

MakeTable_EmpericalDiversityprofile.inc_nose = function(data,q){
  Diversity = Diversity_profile_MLE.inc(data,q)
  #Entropy = Diversity_Tsallis_MLE.inc(Diversity,q)
  output = data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity)
  return(output)
}

MakeTable_Proposeprofile.inc_nose = function(data,q){
  Diversity = Diversity_profile.inc(data,q)
  #Entropy = Diversity_Tsallis.inc(Diversity,q)
  output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity))
  return(output)
}

MakeTable_EmpericalDiversityprofile.inc = function(data, B, q,conf){
  Diversity = bootstrap_forq.inc( data, B,q,conf,Diversity_profile_MLE.inc)
  #Entropy = bootstrap_forq.inc( data, B,q,conf,Diversity_Tsallis_MLE.inc)
  #tmp <- Diversity_Tsallis_MLE.inc(Diversity[,1],q)
  #Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4])
  return(output)
}

MakeTable_Proposeprofile.inc = function(data, B, q,conf){
  Diversity = bootstrap_forq.inc(data, B, q,conf,Diversity_profile.inc)
  #Entropy = bootstrap_forq.inc( data, B,q,conf,Diversity_Tsallis.inc)
  #tmp <- Diversity_Tsallis.inc(Diversity[,1],q)
  #Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4])
  return(output)
}

plot_diversity_profile.inc_nose_yhc <- function(data){
  ggplot(data, aes(x=Order.q, y=Diversity, colour=Community, lty = method))+
  geom_line(size=1.2) +    
  labs(x="Order q", y="Diversity")+
  theme(text=element_text(size=18))

}

plot_diversity_profile.inc_yhc <- function(data){
  ggplot(data, aes(x=Order.q, y=Diversity, colour=Community, lty = method))+
  labs(x="Order q", y="Diversity")+
  theme(text=element_text(size=18))+
  geom_ribbon(aes(ymin=LCL, ymax=UCL, fill=Community, colour=NULL), alpha=0.4)+
  geom_line(size=1.2) 
}

#=======Sample Completeness Curve=========#

sample_coverage = function(freq, q, datatype = c("abundance","incidence_freq")){
  
  if(datatype=="abundance"){
    freq = freq[freq>0]
    n = sum(freq)
    f1 = sum(freq==1)
    f2 = sum(freq==2)
    A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
    
    c_hat = function(q){
      if (q==0){
        
        S_obs = length(freq)
        # f0_hat = if ( f2 == 0 ){( (n-1)/n ) * ( f1*(f1-1)/2 )} else {( (n-1)/n ) * ( (f1^2) / (2*f2) )}
        # f0_hat_star = ceiling(f0_hat)
        # c_hat = S_obs / (S_obs + f0_hat_star)
        f0_hat = ifelse ( f2 == 0 ,( (n-1)/n ) * ( f1*(f1-1)/2 ), ( (n-1)/n ) * ( (f1^2) / (2*f2) ))
        c_hat = S_obs / (S_obs + f0_hat)
        return(c_hat)
        
      } else if (q==1){  
        
        c_hat = 1 - (f1/n)*(1-A)
        return(c_hat)
        
      } else if (q==2){
        
        x = freq[freq>=2]
        c_hat = 1 - (f1/n)*( (A*(1-A))/sum( x*(x-1) / (n*(n-1)) ) )
        return(c_hat)
        
      } else {
        
        r <- 0:(n-1)
        sort.data = sort(unique(freq))
        tab = table(freq)
        term = sapply(sort.data,function(z){
          k=0:(n-z)
          sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
        })
        lambda_hat =  sum(tab*term) + ifelse(f1==0|A==1,0,f1/n*(1-A)^(1-n)*(A^(q-1)-sum(choose(q-1,r)*(A-1)^r)))
        c_hat = 1 - ((f1/n)*(A^(q-1))*(1-A)/lambda_hat)
        return(c_hat)
        
      }
    }
  } else {
    
    t = freq[1]
    freq = freq[-1]; freq = freq[freq>0]
    u = sum(freq)
    Q1 = sum(freq==1)
    Q2 = sum(freq==2)
    B = ifelse(Q2>0,2*Q2/((t-1)*Q1+2*Q2),ifelse(Q1>0,2/((t-1)*(Q1-1)+2),1))
    
    c_hat = function(q){
      if (q==0){
        
        S_obs = length(freq)
        # Chao2 = S_obs + ceiling(if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )})
        Q0_hat = ifelse( Q2 == 0,( (t-1)/t ) * ( Q1*(Q1-1)/2 ), ( (t-1)/t ) * ( (Q1^2) / (2*Q2) ))
        c_hat = S_obs / (S_obs + Q0_hat)
        return(c_hat)
        
      } else if (q==1){  
        
        c_hat = 1 - (Q1/u)*(1-B)
        return(c_hat)
        
      } else if (q==2){
        
        x = freq[freq>=2]
        c_hat = 1 - (t-1)*Q1*( (B*(1-B))/sum( x*(x-1) ) )
        return(c_hat)
        
      } else {
        
        r <- 0:(t-1)
        sort.data = sort(unique(freq))
        tab = table(freq)
        term = sapply(sort.data,function(z){
          k=0:(t-z)
          sum(choose(k-q,k)*exp(lchoose(t-k-1,z-1)-lchoose(t,z)))
        })
        phi_hat = sum(tab*term) + ifelse(Q1==0|B==1,0,Q1/t*(1-B)^(1-t)*(B^(q-1)-sum(choose(q-1,r)*(B-1)^r)))
        c_hat = 1 - ((Q1/t)*(B^(q-1))*(1-B)/phi_hat)
        return(c_hat)
      }
    }
    
  }
  
  sapply(q,c_hat)
  
}

bootstrap_sample = function(freq, B, datatype = c("abundance","incidence_freq")){
  
  if(datatype=="abundance"){
    
    freq = freq[freq>0]
    n = sum(freq)
    f1 = sum(freq == 1)
    f2 = sum(freq == 2)
    
    S_obs = length(freq)
    f0_hat = if ( f2 == 0 ){( (n-1)/n ) * ( f1*(f1-1)/2 )} else {( (n-1)/n ) * ( (f1^2) / (2*f2) )}
    f0_hat_star = ceiling(f0_hat)
    S_hat_Chao1 = S_obs + f0_hat_star
    
    c_hat = if ( f2 != 0 ){ 1 - (f1/n)*((n-1)*f1/(((n-1)*f1)+2*f2))
      
    } else if (f1 != 0) { 1 - (f1/n)*((n-1)*(f1-1)/(((n-1)*(f1-1))+2*f2)) } else { 1 }
    
    lambda_hat = (1-c_hat) / sum((freq/n)*(1-freq/n)^n )
    p_i_hat_obs = (freq/n) * (1-lambda_hat* (1-freq/n)^n ) 
    p_i_hat_unobs = rep( (1-c_hat)/ f0_hat_star, f0_hat_star )
    
    bootstrap_population = c(p_i_hat_obs,p_i_hat_unobs)
    bootstrap_sample = rmultinom(n=B, size=n, prob=bootstrap_population)
    return(bootstrap_sample)
    
  } else {
    
    t = freq[1]
    freq = freq[-1]; freq = freq[freq>0]
    u = sum(freq)
    Q1 = sum(freq==1)
    Q2 = sum(freq==2)
    
    S_obs = length(freq)
    Q_0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )}
    Q_0_hat_star = ceiling(Q_0_hat)
    
    c_hat = if ( Q2 > 0 ){ 1 - (Q1/u)*((t-1)*Q1/(((t-1)*Q1)+2*Q2))
      
    } else { 1 - (Q1/u)*((t-1)*(Q1-1)/(((t-1)*(Q1-1))+2)) }
    
    tau_hat = (u/t) * (1-c_hat) / sum((freq/t)*(1-freq/t)^t )
    pi_i_hat_obs = (freq/t) * (1-tau_hat* (1-freq/t)^t ) 
    pi_i_hat_unobs = rep( (u/t) * (1-c_hat)/ Q_0_hat_star, Q_0_hat_star )
    
    bootstrap_population = c(1,pi_i_hat_obs,pi_i_hat_unobs)
    bootstrap_sample = sapply(1:length(bootstrap_population), function(i) rbinom(n=B, size=t, prob=bootstrap_population[i]))
    bootstrap_sample = if(B==1) {as.matrix(bootstrap_sample)} else {t(bootstrap_sample)}
    return(bootstrap_sample)
    
  }
  
  
  
}

sc_profile.nose = function(freq, q, datatype = c("abundance","incidence_freq")) {
  
  data.frame(Order.q=q, Estimate=sample_coverage(freq, q, datatype))
  
}

sc_profile = function(freq, q, B, conf, datatype = c("abundance","incidence_freq")) {
  
  bootstrap_samples = bootstrap_sample(freq, B, datatype)
  sc_bs = sapply(1:B, function(i) sample_coverage(bootstrap_samples[,i], q, datatype))
  
  LCL = sample_coverage(freq, q, datatype) - qnorm(1-(1-conf)/2)*apply(sc_bs, 1, sd); LCL[LCL<0]=0
  UCL = sample_coverage(freq, q, datatype) + qnorm(1-(1-conf)/2)*apply(sc_bs, 1, sd); UCL[UCL>1]=1
  
  data.frame(Order.q=q, Estimate=sample_coverage(freq, q, datatype), s.e.=apply(sc_bs, 1,sd), LCL=LCL, UCL=UCL)
  
}

plot_sc_profile_nose <- function(data){

      ggplot(data, aes(x=Order.q, y=Estimate, colour=Community))+
      geom_line(size=1.2) +    
      labs(x="Order q", y="sample completeness")+
      theme(text=element_text(size=18), legend.position="bottom")
  
}

plot_sc_profile <- function(data){
  
    ggplot(data, aes(x=Order.q, y=Estimate, colour=Community))+
      labs(x="Order q", y="sample completeness")+
      theme(text=element_text(size=18))+
      geom_ribbon(aes(ymin=LCL, ymax=UCL, fill=Community, colour=NULL), alpha=0.4,show.legend = FALSE)+
      geom_line(size=1.1,show.legend = FALSE)
  
}

#191018 yhc add ggiNEXT
ggiNEXT.iNEXT <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE){
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="order") color.var <- "site"
  if(facet.var=="site") color.var <- "order"
  
  options(warn = -1)
  z <- fortify(x, type=type)
  z$order <- factor(paste0("q = ",z$order),levels = paste0("q = ",unique(z$order)))
  options(warn = 0)
  if(ncol(z) ==7) {se <- FALSE}
  datatype <- unique(z$datatype)
  if(color.var=="none"){
    if(levels(factor(z$order))>1 & "site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep="-")
      
    }else if("site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }else if(levels(factor(z$order))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="order"){     
    z$col <- z$shape <- factor(z$order)
  }else if(color.var=="site"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }else if(color.var=="both"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
  }
  zz=z
  z$method[z$method=="observed"]="Interpolated"
  z$method[z$method=="interpolated"]="Interpolated"
  z$method[z$method=="extrapolated"]="Extrapolated"
  z$lty <- z$lty <- factor(z$method, levels=unique(c("Interpolated", "Extrapolated")))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method=="observed"),]
  
  g <- ggplot(z, aes_string(x="x", y="y", colour="col")) + theme_bw() +
    geom_point(aes_string(shape="shape"), size=3, data=data.sub)
  
  
  g <- g + geom_line(aes_string(linetype="lty"), lwd=1.1) +
    guides(linetype=guide_legend(title="Method"),
           colour=guide_legend(title="Guides"), 
           fill=guide_legend(title="Guides"), 
           shape=guide_legend(title="Guides")) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          legend.key.width = unit(1.2,"cm")) 
  
  if(type==2L) {
    g <- g + labs(x="Number of sampling units", y="Sample coverage")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Sample coverage")
  }
  else if(type==3L) {
    g <- g + labs(x="Sample coverage", y="Species diversity")
  }
  else {
    g <- g + labs(x="Number of sampling units", y="Species diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Species diversity")
  }
  
  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.4)
  
  
  if(facet.var=="order"){
    if(length(levels(factor(z$order))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")      
    }else{
      g <- g + facet_wrap(~order, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$order))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="site"){
    if(!"site"%in%names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }else{
      g <- g + facet_wrap(~site, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$order)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="both"){
    if(length(levels(factor(z$order))) == 1 | !"site"%in%names(z)){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }else{
      g <- g + facet_wrap(site~order) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$site))), byrow=TRUE),
                         fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(grey){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      guides(linetype=guide_legend(title="Method"), 
             colour=guide_legend(title="Guides"), 
             fill=guide_legend(title="Guides"), 
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank())
  }
  # g <- g + theme(legend.box = "vertical")
  return(g)
  
}


#' Fortify method for classes from the iNEXT package.
#'
#' @name fortify.iNEXT
#' @param model \code{iNEXT} to convert into a dataframe.
#' @param data not used by this method
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param ... not used by this method
#' @export
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggplot2::fortify(out1, type=1)

fortify.iNEXT <- function(model, data = model$iNextEst, type = 1, ...) {
  datatype <- ifelse(names(model$DataInfo)[2]=="n","abundance","incidence")
  z <- data
  if(class(z) == "list"){
    z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
    rownames(z) <- NULL
  }else{
    z$site <- ""
  }
  
  if(ncol(z)==6) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if(ncol(z)>6) {
    se <- TRUE
  }
  
  if(type==1L) {
    z$x <- z[,1]
    z$y <- z$qD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }else if(type==2L){
    if(length(unique(z$order))>1){
      z <- subset(z, order==unique(z$order)[1])
    }
    z$x <- z[,1]
    z$y <- z$SC
    if(se){
      z$y.lwr <- z[,8]
      z$y.upr <- z[,9]
    }
  }else if(type==3L){
    z$x <- z$SC
    z$y <- z$qD
    if(se){
      z$y.lwr <- z[,5]
      z$y.upr <- z[,6]
    }
  }
  z$datatype <- datatype
  z$plottype <- type
  if(se){
    data <- z[,c("datatype","plottype","site","method","order","x","y","y.lwr","y.upr")]
  }else{
    data <- z[,c("datatype","plottype","site","method","order","x","y")]
  }
  data
}

