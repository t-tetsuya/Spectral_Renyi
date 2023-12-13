###################################################################################################
###Title: Simulation of ARMA Process and Spectral Simulation.
###Author: Tetsuya Takabatake.
###Last updated: 11 December, 2023.
###################################################################################################
source("SPD.R")
##################################################
sim_AR1 = function(N,theta){
  AR = rep(0,N+1)
  sig = theta[1]; phi = theta[2] 
  AR[1] = rnorm(1) ###Initial Value.
  WN = rnorm(N)
  for(t in 1:N){
    AR[t+1]=phi%*%AR[t]+sig*WN[t]
  }
  return(AR[2:(N+1)])
}

sim_ARMA = function(N,theta,ARMA_order=c(1,0),x0){
  ###Memo: x0 is a vector of initial values of (X_{t-p},...,X_{t-1}).  
  ###Memo: This function generates a sample from an AR(1) process if ARMA_order=c(1,0).
  p = ARMA_order[1]
  q = ARMA_order[2]
  sig = theta[1]
  paraAR = theta[2:(p+1)]
  if( q != 0 ) paraMA = theta[(p+2):(p+q+1)] 
  
  WN = sig*rnorm(q+N) ###White Noise.
  ARMA = rep(0,N+p) 
  ARMA[1:p] = x0 
  for(n in 1:N){
    if( q != 0 ){
      ARMA[p+n] = paraAR%*%rev(ARMA[n:(p+n-1)]) + c(1,paraMA)%*%rev(WN[n:(q+n)])
    }else{
      ARMA[p+n] = paraAR%*%rev(ARMA[n:(p+n-1)]) + WN[n]
    }
  }
  
  return(ARMA[(p+1):(p+N)])
}

spectral_sim <- function(n,N_app,theta,spd){ ###Should take T<=N_app.
  
  a = rep(0,2*N_app) 
  U = rnorm(N_app); V = rnorm(N_app-1)
  
  for(k in 1:(N_app-1)){
    tk = pi*k/N_app 
    spd_tk = spd(tk,theta)
    
    a[k] = (1/2)*(U[k]+1i*V[k])*sqrt(2*pi*spd_tk/N_app)
    a[2*N_app-k+1] = (1/2)*(U[k]-1i*V[k])*sqrt(2*pi*spd_tk/N_app)
  } 
  
  a[N_app] = U[N_app]
  X = fft(a)
  
  return(Re(X)[1:n])
}



