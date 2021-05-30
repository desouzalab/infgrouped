

# mult_simul -------------------------------------------------------------------

#' Simulated Data Multivariated Case
#'
#' @description The function stimulates group data for the multivariate case.
#'
#' @param mm The vector of mu that we want to simulate from.
#' @param ss The covariance matrix for the simulation.
#' @param n_data_sets Integer, representa o numero de data sets a serem
#' simulados.
#' @param breaks_x Vector of values, represents the intervals that x belongs to.
#' @param breaks_y Vector of values, represents the intervals that y belongs to.
#' @param lower_x Vector of values, the lower bounds of x.
#' @param lower_y Vector of values, the lower bounds of y.
#' @param upper_x Vector of values, the upper bounds of x.
#' @param upper_y Vector of values, the upper bounds of y.
#'
#' @return Return a list where each element is a simulation of multivariate data.
#' @export
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#'
#' simulateddata
#'



mult_simul <- function(mm = c(68,68),
                       ss = base::matrix(c(3,2,2,6),2,2) ,
                       ssdata = array(base::rep(0,1000*2*30),
                                      c(1000,2,30)),
                       n_data_sets = 30,
                       breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
                       breaks_y = c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,
                                    71.2,72.2,Inf),
                       lower_x = base::rep(c(-Inf,64,65,66,67,68,69,
                                             70,71,72),10),
                       lower_y = c(base::rep(-Inf,10),
                                   base::rep(64.2,10),
                                   base::rep(65.2,10),
                                   base::rep(66.2,10),
                                   base::rep(67.2,10),
                                   base::rep(68.2,10),
                                   base::rep(69.2,10),
                                   base::rep(70.2,10),
                                   base::rep(71.2,10),
                                   base::rep(72.2,10)),
                       upper_x = base::rep(c(64,65,66,67,68,69,70,
                                             71,72,Inf),10),
                       upper_y = c(base::rep(64.2,10),
                                   base::rep(65.2,10),
                                   base::rep(66.2,10),
                                   base::rep(67.2,10),
                                   base::rep(68.2,10),
                                   base::rep(69.2,10),
                                   base::rep(70.2,10),
                                   base::rep(71.2,10),
                                   base::rep(72,2,10),
                                   base::rep(Inf,10))
){
  
  ssdata = base::array(base::rep(0,1000*2*n_data_sets),
                       c(1000,2,n_data_sets))
  
  for(i in 1: n_data_sets){
    
    ssdata[,,i]<- MASS::mvrnorm(n=1000,mm,ss)
  }
  
  x<- base::matrix(base::rep(0,1000*n_data_sets),ncol=n_data_sets)
  y<- base::matrix(base::rep(0,1000*n_data_sets),ncol=n_data_sets)
  
  for(i in 1:n_data_sets){
    x[,i]<- ssdata[,1,i]
    y[,i]<- ssdata[,2,i]
  }
  
  Freqtable1<- base::array(base::rep(0,10*10*n_data_sets),
                           c(10,10,n_data_sets))
  
  for(i in 1:n_data_sets){
    x_cut<- base::cut(x[,i],breaks= breaks_x)
    y_cut<- base::cut(y[,i],breaks=breaks_y)
    Freqtable1[,,i]<- base::table(x_cut,y_cut)
    
  }
  
  
  simulateddata<- base::array(base::rep(0,n_data_sets*5*100),
                              c(100,5,n_data_sets))
  
  
  for(i in 1:n_data_sets){
    simulateddata[,1,i]<- lower_x
    simulateddata[,3,i]<- lower_y
    simulateddata[,2,i]<- upper_x
    simulateddata[,4,i]<- upper_y
    simulateddata[,5,i]<- c(Freqtable1[,,i])
  }
  
  return(simulateddata)
}





################################################################################



# mexi -------------------------------------------------------------------------

#' ME-XI For First Moment
#'
#' @description Find the first moment of the truncated multivariate
#'  normal distributions for the rectangles which will be used to calculate
#'  the estimate of the covariance matrix.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param mu Mean vector.
#' @param sigma Covariance matrix
#' @return Returns a matrix with the references of the first moments.
#' @examples
#'
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#'
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' out<- em.estimate::mexi(data = simulateddata[,,1],
#'            mu = mu2,
#'            sigma = sigma2)
#'
#'
#' out
#'



mexi<- function(data,mu,sigma){
  n<- base::nrow(data)
  d<- (base::ncol(data)-1)/2

  ind<- base::ncol(data)-1
  lowerb<- NULL
  upperb<- NULL
  for (i in 1:ind){
    if (i%%2!=0) {lowerb<- base::append(lowerb,i)}
    else {upperb<- base::append(upperb,i)}

  }
  data1<- as.matrix(data)



  ex<- base::matrix(base::rep(0,n*d),
                    ncol=d,
                    nrow=n)

  for(i in 1:base::nrow(data)){
    mnts<- tmvtnorm::mtmvnorm(mean=mu,
                              sigma=sigma,
                    lower=data1[i,c(lowerb)],
                    upper=data1[i,c(upperb)])

    ex[i,]<- mnts$tmean

  }
  return(ex)

}

################################################################################

# mcovxi -----------------------------------------------------------------------

#' Mcov - XI E-Step For Variance Estimate
#'
#' @description The function to find the covariance matrices using the
#' current theta of the truncated multivariate normal distributions for the
#' rectangles which will be used to calculate the estimate of the covariance
#' matrix.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param mu Mean vector.
#' @param sigma Covariance matrix
#' @return returns a list containing in each item a variation and estimated
#' covariance matrix.
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' out<- em.estimate::mcovxi(data = simulateddata[,,1],
#'              mu = mu2,
#'              sigma = sigma2)
#'
#' out
#'
#'



mcovxi<- function(data,mu,sigma){
  n<- base::nrow(data)
  d<- (base::ncol(data)-1)/2
  ind<- base::ncol(data)-1
  lowerb<- NULL
  upperb<- NULL

  for (i in 1:ind){
    if (i%%2!=0) {lowerb<- base::append(lowerb,i)}
    else {upperb<- base::append(upperb,i)}

  }
  data1<- as.matrix(data)

  s1<- base::array(base::rep(0,n*d*d),
                   c(d,d,n))

  for(i in 1:base::nrow(data)){
    mnts<- tmvtnorm::mtmvnorm(mean=mu, sigma=sigma,
                    lower=data1[i,c(lowerb)],
                    upper=data1[i,c(upperb)])
    s1[,,i]<- mnts$tvar
    }
  return(s1)
}

################################################################################

# MU:E-step mem ----------------------------------------------------------------

#' Estimate MU verctor for E-step
#'
#' @description The function to calculate the mean vectors.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param exi the first moments matrix.
#' @return returns a vector with the estimated means.
#'
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' exi<- em.estimate::mexi(data=simulateddata[,,1],
#'            mu=mu2,
#'            sigma=sigma2)
#'
#' out<-em.estimate:: mem(data=simulateddata[,,1],
#'           exi=exi)
#'
#' out
#'



mem<- function(data,exi){
  F<- base::ncol(data)
  Freq<- data[,F]
  A<- exi*data[,F]
  mupd<- base::apply(A,2,sum)/base::sum(data[,F])

  return(mupd)

}

################################################################################

# E(XX) ------------------------------------------------------------------------

#' Estimate expected value
#'
#' @description The function to calculate the second moments, for the for the
#'  multivariate case.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param exi The first moments matrix.
#' @param ss1 The array of covariances.
#' @return returns a list where each element is an estimate for the second
#'  moment.
#'
#' @examples 
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#'
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' exi<- em.estimate::mexi(data=simulateddata[,,1],
#'            mu=mu2,
#'            sigma=sigma2)
#'
#' ss1<- em.estimate::mcovxi(data=simulateddata[,,1],
#'              mu=mu2,
#'              sigma=sigma2)
#'
#' out<- em.estimate::exxest(data=simulateddata[,,1],
#'              exi=exi,
#'              ss1=ss1)
#'
#' out
#'
#'


exxest<- function(data,exi,ss1){
  n<- base::nrow(data)
  d<- (base::ncol(data)-1)/2
  exxNew<- base::array(base::rep(0,n*d*d),c(d,d,n))
  for(i in 1:n){
    exxNew[,,i]<- ss1[,,i]+exi[i,]%*%t(exi[i,])
  }
  return(exxNew)

}

################################################################################

# est_sig  ---------------------------------------------------------------------

#' Estimate Sigma in Multivariate Case
#'
#' @description The function to calculate the estimate of the covariances,
#' also for multivariate case.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param exx The list of second moments matrices.
#' @param ex Matrix with the references of the first moments.
#' @param updmu The estimate mu.
#' @return returns a matrix with the variance and covariance estimates.
#'
#' @examples 
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' exi<- em.estimate::mexi(data=simulateddata[,,1],
#'            mu=mu2,
#'            sigma=sigma2)
#'
#' ss1<- em.estimate::mcovxi(data=simulateddata[,,1],
#'              mu=mu2,
#'              sigma=sigma2)
#'
#' exx<- em.estimate::exxest(data=simulateddata[,,1],
#'              exi=exi,
#'              ss1=ss1)
#'
#' updmu <- em.estimate::mem(data=simulateddata[,,1],
#'              exi=exi)
#'
#' out<- em.estimate::est_sig(data=simulateddata[,,1],
#'             exx=exx,
#'             ex=exi,
#'             updmu = updmu)
#'
#'
#' out
#'

est_sig<- function(data,
                 exx,
                 ex,
                 updmu){
  n<- base::nrow(data)
  d<- (base::ncol(data)-1)/2
  s2<- base::array(base::rep(0,n*d*d),c(d,d,n))

  for(i in 1:n){
    s2[,,i]<- (exx[,,i]-updmu%*%t(ex[i,])-ex[i,]%*%t(updmu)+
                 (updmu%*%t(updmu)))*data[i,base::ncol(data)]

  }

  SigmaNew<- base::apply(s2,c(1,2),sum)/base::sum(data[,base::ncol(data)])

  if (base::ncol(SigmaNew)>2){

    for (i in 1:base::ncol(SigmaNew)){
      if (SigmaNew[base::upper.tri(SigmaNew)][i]!=
          SigmaNew[base::lower.tri(SigmaNew)][i]){
        SigmaNew[base::lower.tri(SigmaNew)][i]<-
          SigmaNew[base::upper.tri(SigmaNew)][i]
      }

    }

  }

  return(SigmaNew)

}


################################################################################

# embivg -----------------------------------------------------------------------

#' Estimate Mean Vector and Covariance Matrix
#'
#' @description The EM function, to find the estimate of the mean vector and
#'  covariance matrix.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param mu_init The initial values of mu vector.
#' @param sigma_init The initial covariance matrix.
#' @param maxit Integer, number of max iterations.
#' @param tol1 Stopping criteria for means.
#' @param tol2 Stopping criteria for variances.
#' @return Returns a list with the mu estimate and the variance matrix.
#' @export
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' out<- em.estimate::embivg(data=simulateddata[,,1],
#'              mu_init=mu2,
#'              sigma_init=sigma2,
#'              maxit=1000,
#'              tol1=1e-4,
#'              tol2=1e-3)
#'
#' out
#'

embivg<- function(data,
                  mu_init,
                  sigma_init,
                  maxit=1000,
                  tol1=1e-4,
                  tol2=1e-3){
  flag<- 0
  M_cur<- mu_init
  S_cur<- sigma_init

  for (i in 1:maxit){
    MU_cur<- M_cur
    SS_cur<- S_cur
    MuNew<- em.estimate::mem(data=data,
                exi=em.estimate::mexi(data=data,
                         mu= MU_cur,
                         sigma=SS_cur))
    SigmaNew<- em.estimate::est_sig(data=data,
                     exx=em.estimate::exxest(data=data,
                                exi=em.estimate::mexi(data=data,
                                         mu=MU_cur,
                                         sigma=SS_cur),
                                ss1=em.estimate::mcovxi(data=data,
                                           mu=MU_cur,
                                           sigma=SS_cur)),
                     ex=em.estimate::mexi(data=data,
                             mu=MU_cur,
                             sigma=SS_cur),
                     updmu=em.estimate::mem(data=data,
                               exi=em.estimate::mexi(data=data,
                                        mu=MU_cur,
                                        sigma=SS_cur)))

    Mu_new<- MuNew
    Sigma_new<- SigmaNew
    diff1<- MU_cur-Mu_new
    diff2<- SS_cur-Sigma_new
    D1<- base::abs(base::mean(diff1))
    D2<- base::abs(base::mean(diff2))

    if(D1<tol1 & D2<tol2){
      flag<- 1 ;break
      }

    M_cur<- Mu_new
    S_cur<- Sigma_new

  }

  if(!flag) warning("Didn't Converge \n")

  update<- base::list("mu_estimate" = M_cur,
              "sigma_estimate" = S_cur)

  return(update)

}

################################################################################

# caso trivariate DEVO COLOCAR O exEMPLO? --------------------------------------
# mm<- c(68,68,68) ## The vector of mu that we want to simulate from
# ss<- matrix(c(3,2,1.5,2,4,2.5,1.5,2.5,5),3,3) ## The covariance matrix for the simulation
#
# ssdata<- array(rep(0,3000*3),c(3000,3))
# ssdata<- mvrnorm(n=3000,mm,ss)
#
# x<- rep(0,3000)
# y<- rep(0,3000)
# z<- rep(0,3000)
#
# x<- ssdata[,1]
# y<- ssdata[,2]
# z<- ssdata[,3]
#
# Freqtable1<- array(rep(0,8*8*8),c(8,8,8))
# x_cut<- cut(x,breaks=c(-Inf,65,66,67,68,69,70,71,Inf))
# y_cut<- cut(y,breaks=c(-Inf,63,65,67,69,71,73,75,Inf))
# z_cut<- cut(z,breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf))
#
# Freqtable1<- table(x_cut,y_cut,z_cut)
#
#
# simulateddata3<- array(rep(0,7*512),c(512,7))
# lower_x<- rep(c(-Inf,65,66,67,68,69,70,71),64)
# ly<- c(rep(-Inf,8),rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8))
# lower_y<- rep(ly,8)
# lower_z<- c(rep(-Inf,64),rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),rep(71.5,64),rep(73,64))
#
# upper_x<- rep(c(65,66,67,68,69,70,71,Inf),64)
# uy<- c(rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8),rep(Inf,8))
# upper_y<- rep(uy,8)
# upper_z<- c(rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),
#             rep(71.5,64),rep(73,64),rep(Inf,64))
#
# simulateddata3[,1]<- lower_x
# simulateddata3[,3]<- lower_y
# simulateddata3[,5]<- lower_z
# simulateddata3[,2]<- upper_x
# simulateddata3[,4]<- upper_y
# simulateddata3[,6]<- upper_z
# simulateddata3[,7]<- c(Freqtable1)
#
#
# mu2<- c(67,67,67)# case trivariate??
# sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3) # case trivariate??
#
# out13<- embivg(data=simulateddata3,
#                mu_init=mu2,
#                sigma_init=sigma2,
#                maxit=1000,
#                tol1=1e-4,
#                tol2=1e-3) ## run the embivg funtion to see the
# ## estimate of the parameters using the EM algorithm
# out13 ## see the output
# unlist(out13) ## unlist the output
################################################################################






# zmcem_mult -----------------------------------------------------------------------

#' Generate Samples
#'
#' @description This function generate samples of data for the MCEM algorithm.
#'
#' @param data Data in the form of multivariate grouped data.
#' @param mu The mean vector values of mu.
#' @param sigma The covariance matrix.
#' @return Returns an array with a generated sample
#'
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu <- c(67,67)
#' sigma <- base::matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#'
#' out<- em.estimate::zmcem_mult(data=simulateddata[,,1],
#'                  mu=mu,
#'                  sigma=sigma)
#'
#' out
#'


zmcem_mult<- function(data,
                  mu,
                  sigma){

  k<- 1000
  n<- base::nrow(data)
  d<- (base::ncol(data)-1)/2

  ind<- base::ncol(data)-1
  lowerb<- NULL
  upperb<- NULL
  for (i in 1:ind){
    if (i%%2!=0) {lowerb<- base::append(lowerb,i)}
    else {upperb<- base::append(upperb,i)}

  }
  data1<- as.matrix(data)

  m<- base::array(base::rep(0,k*d*n),c(k,d,n))

  for(i in 1:n){
    m[,,i]<- tmvtnorm::rtmvnorm(k,
                      mean=mu,
                      sigma=sigma,
                      lower=data1[i,c(lowerb)],
                      upper=data1[i,c(upperb)],
                      algorithm="gibbs")
  }

  return(m)
}

################################################################################

# mu_zmcem_mult ----------------------------------------------------------------

#' Estimate Mu for Sample
#'
#' @description This is the function for calculation mean vector from the
#' generated samples.
#'
#' @param data Data in the form of multivariated grouped data.
#' @param generated_sample Array, which is the simulated (generated) data.
#' @return Return the vector of means.
#'
#' @examples
#'
#'
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu <- c(67,67)
#' sigma <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#'
#' generated_sample<- em.estimate::zmcem_mult(data=simulateddata[,,1],
#'                               mu=mu,
#'                               sigma=sigma)
#'
#' out<- em.estimate::mu_zmcem_mult(data=simulateddata[,,1],
#'                      generated_sample=generated_sample)
#'
#'
#' out
#'


mu_zmcem_mult<- function(data,
              generated_sample){

  n<- base::nrow(data)
  F<- base::ncol(data)
  T<- base::sum(data[,F])
  d<- (base::ncol(data)-1)/2
  mm<- base::matrix(base::rep(0,d*n),
                    ncol=d)
  for(i in 1:n){

    mm[i,]<- base::apply(generated_sample[,,i],2,mean)
  }

  Ym<- mm*data[,F]
  sm<- base::rep(0,d)
  for(i in 1:d){
    sm[i]<- base::sum(Ym[,i])
  }

  mu<- (sm/T)

  return(mu)
}


################################################################################

# sigma_zmcem_mult  -------------------------------------------------------------

#' Estimate Covariance Matrix for Sample
#'
#' @description This function calculate the covariance matrix.
#'
#' @param data Data in the form of multivariated grouped data.
#' @param generated_sample Array, which is the simulated (generated) data.
#' @param mu_vector The mean vector from the Mu.
#' @return Return the covariance matrix.
#'
#' @examples
#'
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu <- c(67,67)
#' sigma <- base::matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#'
#' generated_sample<- em.estimate::zmcem_mult(data=simulateddata[,,1],
#'                               mu=mu,
#'                               sigma=sigma)
#'
#' mu_vector<- em.estimate::mu_zmcem_mult(data=simulateddata[,,1],
#'                           generated_sample=generated_sample)
#'
#'
#' out<- em.estimate::sigma_zmcem_mult(data=simulateddata[,,1],
#'                        generated_sample = generated_sample,
#'                        mu_vector = mu_vector)
#'
#' out
#'



sigma_zmcem_mult<- function(data,
                            generated_sample,
                            mu_vector){
  n<- base::nrow(data)
  F<- base::ncol(data)
  total<- base::sum(data[,F])
  k=1000
  d<- (base::ncol(data)-1)/2
  myarray<- base::array(base::rep(0,d*d*k*n),c(d,d,k,n))

  for(i in 1:n){
    for(j in 1:k){

      myarray[,,j,i]<-
        (generated_sample[j,,i]-mu_vector)%*%t(generated_sample[j,,i]-mu_vector)
    }

  }

  Newarray<- base::array(base::rep(0,d*d*n),c(d,d,n))

  for(i in 1:n){

    Newarray[,,i]<- base::apply(myarray[,,,i],c(1,2),mean)
  }

  NNarray<- base::array(base::rep(0,d*d*n),c(d,d,n))

  for(i in 1:n){

    NNarray[,,i]<- Newarray[,,i]*data[i,F]
  }

  S1<- base::apply(NNarray,c(1,2),sum)
  mcov<- S1/total

  return(mcov)
}


################################################################################

# mcem_mult --------------------------------------------------------------------

#' Estimate mean and Variance by MCEM multivariate
#'
#' @description This is the function of M-step, to iterate over the results
#' for mean and variance.
#'
#' @param data Data in the form of multivariated grouped data.
#' @param mu_init The initialized mean vector.
#' @param sigma_init The initialize covariance matrix.
#' @param maxit Integer, number of max iterations .
#' @param tol1 The stopping criterias for means vectors.
#' @param tol2 The stopping criterias for covariance matrix.
#' @return Returns a list with the mu estimate and the variance matrix.
#' @export
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' mu2<- c(67,67)
#' sigma2<- base::matrix(c(3.1,2.2,2.2,4.3),2,2)
#'
#'
#' out<- em.estimate::mcem_mult(data=simulateddata[,,1],
#'               mu_init=mu2,
#'               sigma_init=sigma2,
#'               maxit=1000,
#'               tol1=1e-4,
#'               tol2=1e-3)
#'
#' out
#'


mcem_mult<- function(data,
                   mu_init,
                   sigma_init,
                   maxit=1000,
                   tol1=1e-3,
                   tol2=1e-3){

  flag<- 0
  M_cur<- mu_init
  S_cur<- sigma_init
  iter<- base::rep(0,maxit)

  for (i in 1:maxit){

    MU_cur<- M_cur
    SS_cur<- S_cur

    MuNew<- em.estimate::mu_zmcem_mult(data=data,
                          generated_sample=
                            em.estimate::zmcem_mult(data=data,
                                        mu=MU_cur,
                                        sigma=SS_cur))

    SigmaNew<- em.estimate::sigma_zmcem_mult(data=data,
                                generated_sample=
                                  em.estimate::zmcem_mult(data=data,
                                             mu=MU_cur,
                                             sigma=SS_cur),
                                mu_vector=
                                  em.estimate::mu_zmcem_mult(data=data,
                                      generated_sample=
                                        em.estimate::zmcem_mult(data=data,
                                                            mu=MU_cur,
                                                            sigma=SS_cur)))

    Mu_new<- MuNew
    Sigma_new<- SigmaNew

    diff1<- MU_cur-Mu_new
    diff2<- SS_cur-Sigma_new

    D1<- base::abs(base::mean(diff1))
    D2<- base::abs(base::mean(diff2))

    if(D1<tol1 & D2<tol2){flag<- 1 ;break}

    M_cur<- Mu_new

    S_cur<- Sigma_new
  }

  if(!flag) warning("Didn't Converge \n")

  updatedest = base::list("mu_estimate" = M_cur,
                          "sigma_estimate" = S_cur)
  return(updatedest)

}


################################################################################







# mult_llik --------------------------------------------------------------------

#' Defining the Log-Likelihood Function
#'
#' @description In this function we will calculate the log-likelihood function.
#'
#' @param data Data in the form of multivariated grouped data.
#' @param theta Is a the vector of all parameters, the first d elements
#'  of it are the means, d+1 are v varinces and d+v+1 are c covariances.
#'  In the form:
#'  theta=[mu_1,...,mu_d,
#'  Var_1,...,Var_v,
#'  Cov_1(1,(1+1)), ...,Cov_d(1,d),
#'  Cov_d+1(2,(2+1)), ...,Cov_c(d-1,d)].
#'
#' @return Return the optimum piont or the maximum to log-likelihood
#'
#' @examples
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#' theta<- c(65,65,3,6,2.12132)
#'
#' out<- em.estimate::mult_llik(data = simulateddata[,,1],
#'                 theta = theta)
#'
#' out
#'

mult_llik<- function(data,
                     theta){

  d<- (base::ncol(data)-1)/2
  mu<- theta[c(1:d)]
  sig<- base::matrix(rep(0,d*d),ncol=d)
  subtheta<- theta[-c(1:d)]
  diag(sig)<- subtheta[1:d]
  subtheta1<- subtheta[-c(1:d)]
  sig[base::lower.tri(sig)] <- sig[base::upper.tri(sig)]<- subtheta1
  ind<- base::ncol(data)-1
  lowerb<- NULL
  upperb<- NULL

  for (i in 1:ind){
    if (i%%2!=0) {lowerb<- base::append(lowerb,i)}
    else {upperb<- base::append(upperb,i)}

  }

  data1<- as.matrix(data)
  b<- 0
  for(i in 1:base::nrow(data1)){

    a<- mvtnorm::pmvnorm(lower=data1[i,c(lowerb)],
                upper=data1[i,c(upperb)],
                mean=mu,sigma=sig)

    if (a[1]==0) {a<- 1e-03}
    else {a<- a[1]}
    b<- b+data[i,base::ncol(data)]*base::log(a)

  }

  return(-b)
}


################################################################################

# mle_exact_mult ---------------------------------------------------------------

#' Calculate MLEeexact for Case Multivariate
#'
#' @description Calculate the exact MLE of the parameters (mean, variance and
#' covariance).
#'
#' @param data Data in the form of multivariated grouped data.
#' @param theta Is a the vector of all parameters, the first d elements
#'  of it are the means, d+1 are v varinces and d+v+1 are c covariances.
#'  In the form:
#'  theta=[mu_1,...,mu_d,
#'  Var_1,...,Var_v,
#'  Cov_1(1,(1+1)), ...,Cov_d(1,d),
#'  Cov_d+1(2,(2+1)), ...,Cov_c(d-1,d)].
#'
#' @return Return a vector, for the exact MLE of the parameters for mean,
#'  variance and covariance.
#' @export
#' @examples
#'
#'
#' library(MASS)
#' library(em.estimate)
#' library(tmvtnorm)
#'
#' set.seed(12345)
#'
#' simulateddata = em.estimate::mult_simul(mm = c(68,68),
#'                           ss = base::matrix(c(3,2,2,6),2,2) ,
#'                           n_data_sets = 1,
#'                           breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
#'                           breaks_y = c(-Inf,64.2,65.2,66.2,
#'                                        67.2,68.2,69.2,70.2,
#'                                        71.2,72.2,Inf),
#'                           lower_x = base::rep(c(-Inf,64,65,66,67,
#'                                           68,69,70,71,72),10),
#'                           lower_y = c(base::rep(-Inf,10),
#'                                       base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72.2,10)),
#'                           upper_x = base::rep(c(64,65,66,67,68,69,70,
#'                                                 71,72,Inf),10),
#'                           upper_y = c(base::rep(64.2,10),
#'                                       base::rep(65.2,10),
#'                                       base::rep(66.2,10),
#'                                       base::rep(67.2,10),
#'                                       base::rep(68.2,10),
#'                                       base::rep(69.2,10),
#'                                       base::rep(70.2,10),
#'                                       base::rep(71.2,10),
#'                                       base::rep(72,2,10),
#'                                       base::rep(Inf,10))
#' )
#'
#'
#' theta<- c(65,65,3,6,2.12132)
#'
#' out = em.estimate::mle_exact_mult(theta<- c(65,65,3,6,2.12132),
#'                      data=simulateddata[,,1])
#'
#' out
#'
#'


mle_exact_mult<- function(data,
                          theta){

  MLEexact<- stats::optim(theta<- theta,
                   fn=mult_llik,
                   data=data,
                   method="Nelder-Mead")


  return(MLEexact$par)

}

################################################################################



