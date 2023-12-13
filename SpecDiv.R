###################################################################################################
###Title: Gradient Descent for Spectral Divergences (Spectral Renyi and Pinsker-Itakura-Saito).
###Author: Tetsuya Takabatake.
###Last updated: 11 December, 2023.
###################################################################################################
###0.Define Preliminary Functions.
norm = function(x) sqrt(sum(x^2)) 
pardv.z = function(f,x,z,h=1e-4){
  len.x = length(x); len.z = length(z)
  pardv.z.res = matrix(0,nrow=len.z,ncol=len.x)
  for( j in 1:len.x ){
    for( k in 1:len.z ){
      hk = rep(0,len.z)
      hk[k] = h
      pardv.z.res[k,j] = {f(x[j],z+hk)-f(x[j],z)}/h
    } 
  } 
  return(pardv.z.res)
} 
##################################################
###1.Define Spectral Divergences.
###1.1.Preliminary Functions.
m_func_PIS = function(x) x^{-1}-1+log(x) 
m_func_SpecRenyi = function(x,alp) (alp-1)^{-1}*{alp*log(x) -log(alp*x+(1-alp))}
###1.2.Define Several Divergences.
Div_func = function(theta,data,spd,log_spd=0,inv_spd=0,dv=0,method.div,smooth.In,alp=0){
  ###Choose method.div from "PIS" and "SpecRenyi".
  ###Remark: Choose alp\in(0,1) if the Spectral Renyi Divergence is used.
  n = length(data)
  lamb = (2*pi/n)*seq(1,floor(n/2))
  if( smooth.In == F ){
    In = {(2*pi*n)^{-1}*abs(fft(data))^2}[2:(floor(n/2)+1)]
  }else{
    In = (2*pi)^{-1}*spec.pgram(ts(data),spans=c(3,5),plot=F)$spec
  }
  lamb = lamb[-floor(n/2)]
  In = In[-floor(n/2)]
  spd.para = spd(lamb,theta,dv=0)
  
  if(dv == 0){
    if( method.div == "PIS" ) Div = n^{-1}*2*sum( m_func_PIS(spd.para/In) )
    if( method.div == "SpecRenyi" ) Div = n^{-1}*2*sum( m_func_SpecRenyi(spd.para/In,alp) )
  }
  
  if(dv == 1){
    Div = rep(0,length(theta))
    if( method.div == "PIS" || method.div == "SpecRenyi" ){
      log_spd_dv = log_spd(lamb,theta,dv=1)
    }
    for(p in 1:length(theta)){
      if( method.div == "PIS" ){
        Div[p] = n^{-1}*2*sum( (1-In/spd.para)*log_spd_dv[p,] ) 
      } 
      if( method.div == "SpecRenyi" ){
        Div[p] = {alp/(alp-1)}*n^{-1}*2*sum( { 1-spd.para/{alp*spd.para+(1-alp)*In} }*log_spd_dv[p,] ) 
      }
    }
  }
  return(Div)
}
##################################################
###2.Gradient Descent for Divergences.
GD_Div = function(data,theta_ini,spd,log_spd,inv_spd,method.div,smooth.In,sim_spec=F,alp=0,max_iter_step=10000){ ###GD:Gradient Descent
  dim_para = length(theta_ini)
  res_GD = matrix(0,nrow=2*dim_para,ncol=max_iter_step)
  
  LearningRate = 0.005
  eps.GD = 1e-3
  cou_iter = 0
  
  theta_m = theta_ini
  Div_dv = Div_func(theta_m,data,spd,log_spd,inv_spd,dv=1,method.div,smooth.In,alp)
  while( norm(Div_dv) >= eps.GD  ){
    if( cou_iter >= max_iter_step ) break 
    Div_dv = Div_func(theta_m,data,spd,log_spd,inv_spd,dv=1,method.div,smooth.In,alp)
    theta_m = theta_m - LearningRate*Div_dv
    cou_iter = cou_iter+1
    res_GD[1:dim_para,cou_iter] = theta_m
    res_GD[dim_para+1:dim_para,cou_iter] = Div_dv
  }
  
  res_GD = res_GD[,1:cou_iter]
  return(res_GD)
}
##################################################

