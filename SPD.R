###################################################################################################
###Title: Spectral Density Functions.
###Author: Tetsuya Takabatake.
###Last updated: 11 December, 2023.
###################################################################################################
###Numerical computation of an inverse function.
Newton_method = function(func,x0,lower,upper,eps_Newton=1e-1000,max_iter_step=1e+5,increasing=T){
  x1 = x0
  x2 = x1 - func(x1)/func(x1,dv=1)
  cou = 0
  if(increasing==T) sign_Newton = 1 else sign_Newton = -1
  while( abs(x2-x1) >= eps_Newton ){
    if( x2 <= lower ) x2 = x1 - sign_Newton*(x1-lower)/2
    if( x2 >= upper ) x2 = x1 + sign_Newton*(upper-x1)/2
    if( cou >= max_iter_step || abs(x2-x1) <= eps_Newton ) break 
    x1 = x2
    x2 = x1 - func(x1)/func(x1,dv=1)
    cou = cou+1
  }
  return(x2)
}
##################################################
##################################################
###1.SPDs for AR processes.
###1.1.Case of the AR1 process.
spd_AR1 = function(x,theta,dv=0){ ###dv:the order of the derivative. 
  sig=theta[1]; rho=theta[2] 
  if(dv == 0) spd = (sig^2/(2*pi))*(1-2*rho*cos(x)+rho^2)^{-1}
  if(dv == 1){
    spd = matrix(0,nrow=length(theta),ncol=length(x))
    for( j in 1:length(x) ){
      spd[1,j] = (sig/pi)*(1-2*rho*cos(x[j])+rho^2)^{-1}
      spd[2,j] = (sig^2/(2*pi))*(1-2*rho*cos(x[j])+rho^2)^{-2}*{-1*(-2*cos(x[j])+2*rho)}
    }
  }
  return(spd)
}

log_spd_AR1 = function(x,theta,dv){ ###dv:the order of the derivative. 
  sig=theta[1]; rho=theta[2]
  if(dv == 0) log_spd = 2*log(sig) -log(2*pi) -log(1-2*rho*cos(x)+rho^2)
  if(dv == 1){
    log_spd = matrix(0,nrow=length(theta),ncol=length(x))
    for( j in 1:length(x) ){
      log_spd[1,j] = 2/sig
      log_spd[2,j] = -(-2*cos(x[j])+2*rho)/(1-2*rho*cos(x[j])+rho^2)
    }
  }
  return(log_spd)
}

inv_spd_AR1 = function(x,theta,dv=1){
  sig=theta[1]; rho=theta[2]
  if( dv == 0) sig^{-2}*(2*pi)*(1-2*rho*cos(x)+rho^2)
  if( dv == 1 ){
    inv_spd = matrix(0,nrow=length(theta),ncol=length(x))
    for( j in 1:length(x) ){
      inv_spd[1,j] = -2*sig^{-3}*(2*pi)*(1-2*rho*cos(x[j])+rho^2)
      inv_spd[2,j] = (sig^2/(2*pi))^{-1}*(-2*cos(x[j])+2*rho)
    }
  }
}

var_AR1 = function(N,theta){
  sig=theta[1]; rho=theta[2]
  acvf = rep(0,N)
  acvf[1] = sig^2/(1-rho^2)
  for(tau in 1:(N-1)){ acvf[tau+1]=rho^{tau}*acvf[1] }
  return(toeplitz(acvf))
}
##############################
phi_AR1_func = function(rho,dv=0){
  phi_AR1 = rep(0,length(rho))
  for( j in 1:length(rho) ){
    if( dv == 0 ){
      if( rho[j]>=0 && rho[j]<1 ) phi_AR1[j] = 1/(1-rho[j])-1
      if( rho[j]>-1 && rho[j]<0 ) phi_AR1[j] = 1/(-1-rho[j])+1
    }
    if( dv == 1 ){
      if( rho[j]>=0 && rho[j]<1 ) phi_AR1[j] = (1-rho[j])^{-2} 
      if( rho[j]>-1 && rho[j]<0 ) phi_AR1[j] = (-1-rho[j])^{-2} 
    }
  }
  return( phi_AR1 )
} 

inv_phi_AR1_func = function(phi){
  phi_AR1_func_zero = function(rho,dv=0){
    if(dv==0) phi_AR1_zero = phi_AR1_func(rho)-phi
    if(dv==1) phi_AR1_zero = phi_AR1_func(rho,dv=1)
    return( phi_AR1_zero )
  } 
  inv_phi_AR1 = Newton_method(phi_AR1_func_zero,x0=0,lower=-1,upper=1) 
  return( inv_phi_AR1 )
} 

spd_AR1_rep = function(x,theta_rep){
  sig = exp( theta_rep[1] ) 
  rho = inv_phi_AR1_func( theta_rep[2] ) 
  theta = c(sig,rho)
  spd = spd_AR1(x,theta)
  return(spd)
}

inv_phi_AR1_func_dv = function(phi){
  rho = inv_phi_AR1_func(phi) 
  inv_phi_dv = 1/phi_AR1_func(rho,dv=1)
  return( inv_phi_dv )
} 

log_spd_AR1_rep = function(x,theta_rep,dv){ ###dv:the order of the derivative. 
  
  sig = exp( theta_rep[1] ) 
  rho = inv_phi_AR1_func( theta_rep[2] ) 
  theta = c(sig,rho)
  
  if(dv == 0) log_spd_rep = log_spd_AR1(x,theta,dv=0)
  
  if(dv == 1){
    dim_theta = length(theta) 
    log_spd_rep = matrix(0,nrow=dim_theta,ncol=length(x))
    inv_phi_dv = inv_phi_AR1_func_dv(theta_rep[2])
    for( j in 1:length(x) ){
      log_spd_rep[1,j] = 2
      log_spd_rep[2,j] = -1*(-2*cos(x[j])*inv_phi_dv+2*rho*inv_phi_dv)/(1-2*rho*cos(x[j])+rho^2)
    }
  }
  return(log_spd_rep)
}
##############################
###1.2.Case of the ARMA(p,q) process.
###Compute the transfer function of the ARMA(p,q) process.
TrasFuncARMA = function(x,paraAR,paraMA){
  
  z = exp(-1i*x)
  
  TrasFuncAR.z = 1
  if( max(abs(paraAR)) > 0 ){
    p = length(paraAR)
    for( p1 in 1:p ) TrasFuncAR.z = TrasFuncAR.z - paraAR[p1]*z^{p1}
  }
  
  TrasFuncMA.z = 1
  if( max(abs(paraMA)) > 0 ){
    q = length(paraMA)
    for( q1 in 1:q ) TrasFuncMA.z = TrasFuncMA.z + paraMA[q1]*z^{q1}
  }
  
  TrasFuncARMA.z = TrasFuncMA.z/TrasFuncAR.z
  
  return(TrasFuncARMA.z)
}

###Compute the SPD for the ARMA(p,q) process.
spd_ARMA = function(x,theta,ARMA_order=c(1,0)){ 
  
  p = ARMA_order[1]
  q = ARMA_order[2]
  sig = theta[1]
  paraAR = theta[2:(p+1)]
  if( q != 0 ) paraMA = theta[(p+2):(p+q+1)] else paraMA = 0 
  
  spd_WN = sig^2/(2*pi)
  PowerTrasFuncARMA = abs(TrasFuncARMA(x,paraAR,paraMA))^2
  spd = spd_WN*PowerTrasFuncARMA
  
  return( spd )
}
####################
###1.2.1.SPD for the AR(p) process.
spd_AR = function(x,paraAR,sig) spd_ARMA(x,paraAR,paraMA=0,sig)
###1.2.2.SPD for the MA(q) process.
spd_MA = function(x,paraMA,sig) spd_ARMA(x,paraAR=0,paraMA,sig)
###1.2.3.SPD for the ARMA process with reparameterized parameter.
spd_ARMA_rep = function(x,theta_rep,ARMA_order=c(1,0)){
  sig = exp(theta_rep[1])
  theta = c(sig,theta_rep[-c(1)])
  spd = spd_ARMA(x,theta,ARMA_order)
  return( spd )
}
####################
log_spd_ARMA = function(x,theta,ARMA_order,dv){ ###dv:the order of the derivative. 
  
  p = ARMA_order[1]
  q = ARMA_order[2]
  sig = theta[1]
  paraAR = theta[2:(p+1)]
  if( q != 0 ) paraMA = theta[(p+2):(p+q+1)] else paraMA = 0 
  
  TrasFuncAR = function(x) TrasFuncARMA(x,paraAR,paraMA=0)
  TrasFuncMA = function(x) TrasFuncARMA(x,paraAR=0,paraMA)
  
  if(dv == 0){
    log_spd = 2*log(sig) + log(abs(TrasFuncMA)^2) - log(abs(TrasFuncAR)^2)
  } 
  
  if(dv == 1){
    dim_theta = length(theta)
    log_spd = matrix(0,nrow=dim_theta,ncol=length(x))
    for( j in 1:length(x) ){
      log_spd[1,j] = 2/sig
      for( p0 in 2:(p+1) ){
        log_spd[p0,j] = -abs(TrasFuncAR(x[j]))^{-2}*( 
          2*Re(TrasFuncAR(x[j]))*(-cos(-p0*x[j])) + 2*Im(TrasFuncAR(x[j]))*(-sin(-p0*x[j]))
          )
        }
      if( q > 0 ){
        for( q0 in (p+2):(p+q+1) ){
          log_spd[q0,j] = abs(TrasFuncMA(x[j]))^{-2}*(
            2*Re(TrasFuncMA(x[j]))*cos(-q0*x[j]) + 2*Im(TrasFuncMA(x[j]))*sin(-q0*x[j])
            )
        }
      }
    }
  }
  return( log_spd )
}


log_spd_ARMA_rep = function(x,theta_rep,ARMA_order,dv){ ###dv:the order of the derivative. 
  sig = exp( theta_rep[1] )
  theta = c(sig,theta_rep[-c(1)])
  log_spd_rep = log_spd_ARMA(x,theta,ARMA_order,dv)
  return( log_spd_rep )
}
##################################################
##################################################
###4.Brune's Model.
###4.1.SPD for Brune's Model wrt theta.
spd_route_func = function(x,Q){
  spd = exp(-abs(x)/Q)
  return(spd)
}

spd_Brune = function(x,theta,para_known){
  dim_theta = length(theta)
  sig = theta[1]
  omega_c = theta[2]
  m = para_known[1]
  spd = sig^2*abs(x)^2*( 1 + (abs(x)/omega_c)^m )^{-2}
  
  if(dim_theta==3){
    Q = theta[3] ###theta_rep[3]=log(Q).
    spd_route = spd_route_func(x,Q)
    spd = spd*spd_route
  }
  
  return(spd)
}
##################################################
###4.2.SPD for Brune's Model wrt theta_rep.
phi_Brune_func = function(omega_c,dv=0){
  if( dv==0 ) phi_Brune = log(omega_c/(pi-omega_c)) 
  if( dv==1 ){
    phi_Brune = omega_c^{-1} + (pi-omega_c)^{-1}
  }
  return( phi_Brune )
}

inv_phi_Brune_func = function(phi,dv=0){
  if( dv==0 ) inv_phi_Brune = pi*exp(phi)/(1+exp(phi)) 
  if( dv==1 ){
    inv_phi_Brune = pi*exp(phi)/(1+exp(phi))^2 
  } 
  return(inv_phi_Brune)
} 

spd_Brune_rep = function(x,theta_rep,para_known){
  dim_theta = length(theta_rep)
  sig = exp( theta_rep[1] ) 
  omega_c = inv_phi_Brune_func(theta_rep[2]) 
  if(dim_theta==3){
    Q = exp( theta_rep[3] ) 
    theta = c(theta,Q)
  }
  spd = spd_Brune(x,theta,para_known)
  return(spd)
}
##################################################
###4.3.Partial Derivatives of SPD for Brune's Model wrt theta_rep.
log_spd_Brune = function(x,theta,para_known,dv){ ###dv:the order of the derivative. 
  
  dim_theta = length(theta)
  sig = theta[1]
  omega_c = theta[2]
  if(dim_theta==3) Q = theta[3]
  m = para_known[1]

  if(dv == 0){
    log_spd = 2*log(sig) +2*log(x) -2*log(1+(x/omega_c)^m)
    if(dim_theta==3) log_spd = log_spd -x/Q
  } 
  
  if(dv == 1){
    log_spd = matrix(0,nrow=dim_theta,ncol=length(x))
    for( j in 1:length(x) ){
      log_spd[1,j] = 2/sig
      log_spd[2,j] = -2*{1+(x[j]/omega_c)^m}^{-1}*{x[j]^m*(-m*omega_c^{-m-1})}
      if(dim_theta==3) log_spd[3,j] = x[j]/Q^2
    }
  }
  return(log_spd)
}


log_spd_Brune_rep = function(x,theta_rep,para_known,dv){ ###dv:the order of the derivative. 
  
  dim_theta = length(theta_rep)
  sig = exp( theta_rep[1] ) 
  omega_c = inv_phi_Brune_func(theta_rep[2]) 
  theta = c(sig,omega_c)
  if(dim_theta==3){
    Q = exp( theta_rep[3] )
    theta = c(theta,Q)
  } 
  m = para_known[1]
  
  if(dv == 0) log_spd_rep = log_spd_Brune(x,theta,para_known,dv=0)
    
  if(dv == 1){
    log_spd_rep = matrix(0,nrow=dim_theta,ncol=length(x))
    for( j in 1:length(x) ){
      log_spd_rep[1,j] = 2
      log_spd_rep[2,j] = -2*{1+(x[j]/omega_c)^m}^{-1}*{ x[j]^m*(-m*omega_c^{-m-1})*inv_phi_Brune_func(theta_rep[2],dv=1) }
      if(dim_theta==3) log_spd_rep[3,j] = x[j]*exp(-theta_rep[3])
    }
  }
  return(log_spd_rep)
}



