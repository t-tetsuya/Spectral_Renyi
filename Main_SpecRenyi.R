###################################################################################################
###Title: Numerical Study of the Spectral Renyi Divergence (ARMA, Brune).
###Author: Tetsuya Takabatake.
###Last updated: 11 December, 2023.
###################################################################################################
###0.Preliminary Results.
source("SPD.R")
source("SimTS.R")
source("SpecDiv.R") 
##################################################
###0.Select the Setup of the Numerical Study.
###0.1.Choose Model = "ARMA" or "Brune".
###0.1.1.Set Model = "ARMA" if the data generating process is an ARMA process. 
#Model = "ARMA" 
###0.1.2.Set Model = "Brune" if the data generating process is Brune's model. 
Model = "Brune" 
###0.2.Choose Trend = TRUE or FALSE. 
###Remark: Set Trend = T if some trigonometric trends are contaminated in observations.
Trend = T　
###0.3.Select Other Settings.
###Select the number of replications of the simulation.
n_rep = 1 
###Select spectral divergences used for estimating parameters.
###Remark: Choose elements of the vector from "SpecRenyi" (Spectral Renyi Divergence) and "PIS" (Pinsker-Itakura-Saito Divergence). 
method.div_vec = c("PIS","SpecRenyi")
###Select the value of alpha used in spectral Renyi divergence. 
###Remark: Set alpha = 1 if the PIS divergence is used.
alp_vec = c(1,0.9)
###Select the value of alpha used in the spectral Renyi divergence. 
###Remark: Set alpha = 1 if the PIS divergence is used.
alp_vec = c(1,0.9)
###Select whether smoothing methods for the spectral density estimates are used or not used. 
###Remark: Select TRUE (resp.FALSE) if the smoothing method is used (resp. not used). 
smooth.In_vec = c(F,T)
##################################################
###1.Setup of the Model Parameters (ARMA, Brune). 
n = 2^10 ###2^10=1024.
delta = 1
##############################
###1.1.ARMA process.
if( Model == "ARMA" ){
  ###1.2.1.Known Parameters.
  p = 2 ###Order of the AR Components.
  q = 0 ###Order of the MA Components. (Remark: q must be less than p.)
  ###1.2.2.Unknown Parameters.
  ###Parameter of the MA Coefficients.
  sig = 2
  ###Parameter of the AR Coefficients.
  ###Case of AR1 Process.
  if( p == 1 ){
    paraAR = c(-0.1)
  }
  ###Case of the AR2 Process
  if( p == 2 ){
    paraAR = c(0.9,-0.2) 
  }
  ###Case of AR3 Process.
  if( p == 3 ){
    paraAR = c(0.2,0.3,-0.4) 
  }
  ####################
  ###1.2.3.Preliminary Calculations of theta.
  if( q == 0 ){
    theta = c(sig,paraAR)
    if( p == 1 ){ 
      theta_rep = c(log(theta[1]),phi_AR1_func(theta[2]))
    }
    if( p > 1 ){
      theta_rep = c(log(sig),theta[-c(1)])
    }
  }
  if( q > 1 ){
    paraMA = c(0.1)
    theta=c(sig,paraAR,paraMA)
    theta_rep = c(log(sig),theta[-c(-1)])
  }
  ARMA_order = c(p,q)
}
##############################
###1.2.Brune model.
if( Model == "Brune" ){
  ###Select Route = T if a route effect is taking into account.
  Route = T 
  ###1.2.1.Unknown Parameters.
  sig = 3
  omega_c = 2
  Q = 1/2
  ###1.2.2.Known Parameter.
  m = 2 
  ###1.2.3.Preliminary calculations of theta and theta_rep.
  theta = c(sig,omega_c,Q)
  theta_rep = c(log(sig),phi_Brune_func(omega_c),log(Q))
  para_known = c(m)
}
##################################################
###2.Select the Initial Value of the Gradient Descent.
###2.1.ARMA Process.
if( Model == "ARMA" ){
  ###2.1.1.Select the Initial Value of theta (ARMA).
  if( q == 0 ){
    ###Case of the AR1 Process.
    if( p == 1 ){
      theta_ini = c(10,-0.5)
    }
    ###Case of the AR2 Process.
    if( p == 2 ){
      theta_ini = c(0.5,0,0)
    }
    ###Case of the AR3 Process.
    if( p == 3 ){
      theta_ini = c(0.1,-0.3,-0.2,0.1)
    }
  }else{
    theta_ini = c(0.1,paraAR,paraMA)
  }
  ###2.1.2.Compute the Initial Value of theta_rep (ARMA).
  if( p == 1 ) theta_rep_ini = c(log(theta_ini[1]),phi_AR1_func(theta_ini[2]))
  if( p > 1 ){
    theta_rep_ini = c(log(theta_ini[1]),theta_ini[-c(1)])
  }
}
##############################
###2.2.Brune Model.
if( Model == "Brune" ){
  ###2.2.1.Select the Initial Value of theta (Brune Model).
  theta_ini = c(0.5,0.5,1) 
  ###2.2.2.Compute the Initial Value of theta_rep (Brune Model).
  theta_rep_ini = c(log(theta_ini[1]),phi_Brune_func(theta_ini[2]),log(theta_ini[3]))
  ###2.2.3.Compute the Unknown Parameter and Initial Values when Route == FALSE.
  if( Route == FALSE ){
    theta = theta[1:2]
    theta_rep = theta_rep[1:2]
    theta_ini = theta_ini[1:2]
    theta_rep_ini = theta_rep_ini[1:2]
  }
}
####################
result_theta_hat = matrix(0,nrow=n_rep*length(method.div_vec),ncol=1+length(theta))

if( Model == "ARMA" ){
  if( p == 1 && q == 0 ) colnames(result_theta_hat) = c("alpha","sigma","rho")
  if( p == 2 && q == 0 ) colnames(result_theta_hat) = c("alpha","sigma","rho1","rho2")
} 
if( Model == "Brune" ){
  if( Route == TRUE ) colnames(result_theta_hat) = c("alpha","sigma","omega_c","Q")
  if( Route == FALSE ) colnames(result_theta_hat) = c("alpha","sigma","omega_c")
} 

row_names_temp = NULL
for( cou_method in 1:length(method.div_vec) ){
  method.div = method.div_vec[cou_method]
  row_names_temp = c(row_names_temp,paste(method.div,seq(1,n_rep),sep=""))
}
rownames(result_theta_hat) = row_names_temp
##################################################
###3.Estimation Step.
for( cou_rep in 1:n_rep ){
  ###3.1.Generate the data (AR1,WN,RespAR1,Noisy)
  if( Model == "ARMA" ){
    ###Select sim_spec = TRUE if the spectral simulation method is used to generate observations.
    sim_spec = F
    ###3.1.1.Preliminary Calculations of theta and thera_rep.
    spd_func = function(x,theta,dv=0) spd_ARMA(x,theta,ARMA_order)
    log_spd_func = function(x,theta_rep,dv) log_spd_ARMA(x,theta,ARMA_order,dv)
    if( sim_spec == FALSE ){
      data0 = sim_ARMA(n,theta,ARMA_order,x0=rnorm(p)) 
    }else{
      data0 = spectral_sim(n,N_app=n*2^7,theta=theta,spd=spd_func)
    }
    ###Memo:x0 is an initial value of the ARMA process.
    ylim.spd = c(0,200)
    ###3.1.2.Define the reparametrized spectral density function of the AR1 model.
    if( p == 1 ){
      spd_func_rep = function(x,theta_rep,dv=0) spd_AR1_rep(x,theta_rep)
      log_spd_func_rep = function(x,theta_rep,dv) log_spd_AR1_rep(x,theta_rep,dv)
    }
    
    if( p > 1 ){
      spd_func_rep = function(x,theta_rep,dv=0) spd_ARMA_rep(x,theta_rep,ARMA_order)
      log_spd_func_rep = function(x,theta_rep,dv) log_spd_ARMA_rep(x,theta_rep,ARMA_order,dv)
    }
  } 
  
  if( Model == "Brune" ){
    ###Select sim_spec = TRUE if the spectral simulation method is used to generate observations.
    sim_spec = T
    ###1.1.Generate the data of Brune's model.
    spd_func = function(x,theta,dv=0) spd_Brune(x,theta,para_known)
    log_spd_func = function(x,theta_rep,dv) log_spd_Brune(x,theta_rep,para_known,dv)
    data0 = spectral_sim(n,N_app=n*2^7,theta=theta,spd=spd_func)
    ylim.spd = c(0,2)
    ###1.2.Define the reparametrized spectral density function of Brune's model.
    spd_func_rep = function(x,theta_rep,dv=0) spd_Brune_rep(x,theta_rep,para_known)
    log_spd_func_rep = function(x,theta_rep,dv) log_spd_Brune_rep(x,theta_rep,para_known,dv)
  }
  ##################################################
  ###2.Compute the Gradient Descent.
  ###2.1.Preliminaries.
  ord.rou = 4
  if( Trend == FALSE ) data = data0
  ###Case1: ARMA Process.
  if( Model == "ARMA" && Trend == TRUE ){
    A1 = 1*pi/4
    z1 = 20
    A2 = pi/8
    z2 = 20   
    trend1 = z1*sqrt(2*pi/n)*sin(A1*seq(1,n))
    trend2 = z2*sqrt(2*pi/n)*sin(A2*seq(1,n))
    data = trend1 + trend2 + data0
  }
  ###Case2: Brune's Model.
  if( Model == "Brune" && Trend == TRUE ){
    A1 = 1*pi/4
    z1 = 3
    A2 = pi/8
    z2 = 3
    trend1 = z1*sqrt(2*pi/n)*sin(A1*seq(1,n))
    trend2 = z2*sqrt(2*pi/n)*sin(A2*seq(1,n))
    data = trend1 + trend2 + data0
  }
  ##################################################
  ###2.2.Computations of the Spectral Divergences (ARMA Process, Brune's model).
  for( cou_method in 1:length(method.div_vec) ){
    method.div = method.div_vec[cou_method]
    smooth.In = smooth.In_vec[cou_method]
    alp = alp_vec[cou_method]
    
    reparametrization = F
    if( reparametrization == TRUE ){
      result_GD = GD_Div(data,theta_rep_ini,spd_func_rep,log_spd_func_rep,inv_spd=0,method.div,smooth.In,sim_spec,alp)
      num_iter = ncol(result_GD)
      theta_rep_hat = result_GD[1:length(theta),num_iter]
      
      if( Model == "ARMA" ){
        sig_hat = exp(theta_rep_hat[1])
        if( p == 1 ){
          rho_hat = inv_phi_AR1_func(theta_rep_hat[2])
          theta_hat = c(sig_hat,rho_hat)
        }
        if( p > 1 ) theta_hat = c(sig_hat,theta_rep_hat[-c(1)])
      }
      
      if( Model == "Brune" ){
        sig_hat = exp(theta_rep_hat[1])
        omega_c_hat = inv_phi_Brune_func(theta_rep_hat[2])
        theta_hat = c(sig_hat,omega_c_hat)
        if( Route == T ){
          Q_hat = exp(theta_rep_hat[3])
          theta_hat = c(theta_hat,Q_hat)
        }
      }
    }else{
      result_GD = GD_Div(data,theta_ini,spd_func,log_spd_func,inv_spd=0,method.div,smooth.In,sim_spec,alp)
      num_iter = ncol(result_GD)
      theta_hat = result_GD[1:length(theta),num_iter]
    }
    
    cou_row = n_rep*( cou_method-1 ) + cou_rep
    result_theta_hat[cou_row,] = c(alp,theta_hat)
    print(paste("cou=",cou_rep,"_method=",method.div,"_alp=",alp,"_smooth.In=",smooth.In,"_theta_hat",seq(1,length(theta)),"=",round(theta_hat,ord.rou),sep=""))
  }
  print(paste("TrueValue_theta",seq(1,length(theta)),"=",theta,sep=""))
}
##################################################
result_theta_hat







