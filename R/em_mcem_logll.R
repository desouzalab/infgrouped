
### E-Step For Mean Estimate ---------------------------------------------------

#' E-Step For Mean Estimate
#'
#' @param theta Initial theta value, vector containing mean and sd.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param freq frequency over the intervals, values in vector.
#' @return Return M which are the estimates of mean in E-step.
#' @examples
#'
#'  library(infgrouped)
#'  library(stats)
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  munew <- infgrouped:::mest(theta=c(67,2),
#'                bl=simdataaaa$simul_data[,1,1],
#'                bu=simdataaaa$simul_data[,2,1],
#'                freq=simdataaaa$simul_data[,3,1])
#'
#'
#'  munew



mest<- function(theta,bl,bu,freq){

  Aj<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))

  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
  }

  dinom<- NULL
  for(i in 1:base::length(bl)){
    dinom[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    if(dinom[i]==0) {
      Aj[i]<- theta[1]-theta[2]*((stats::dnorm(bstar[i])-
                                    stats::dnorm(astar[i]))/0.0001)
      }
    else {
      Aj[i]<- theta[1]-theta[2]*((stats::dnorm(bstar[i])-
                                    stats::dnorm(astar[i]))/dinom[i])
      }
  }

  M <- base::sum(Aj*freq)/base::sum(freq)
  return(M)

}
################################################################################

### E-Step For Variance Estimate -----------------------------------------------

#' E-Step For Variance Estimate
#'
#' @param theta Initial theta value, vector containing mean and sd.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param freq Frequency over the intervals, values in vector.
#' @param muupdate Is the updated estimates of mu.
#' @return Return ss which are the estimates of sigma (variance) in E-step.
#' @examples
#'
#'  library(infgrouped)
#'  library(stats)
#'
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                         n=50,
#'                         nclass = 10,
#'                         mean = 68,
#'                         sd = 1.80,
#'                         fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  ssnew <- infgrouped:::ssest(theta=c(67,2),
#'               bl=simdataaaa$simul_data[,1,1],
#'               bu=simdataaaa$simul_data[,2,1],
#'               muupdate=mest(theta=c(67,2),
#'                             bl=simdataaaa$simul_data[,1,1],
#'                             bu=simdataaaa$simul_data[,2,1],
#'                             freq=simdataaaa$simul_data[,3,1]),
#'               freq=simdataaaa$simul_data[,3,1])
#'
#'  ssnew


ssest<- function(theta,bl,bu,muupdate,freq){

  Bj<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))

  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]

  }

  astar[1] <- (-1000)
  bstar[base::length(bl)]<- (1000)


  for(i in 1:base::length(bl)){
    dinom[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    if(dinom[i]==0) {
      Bj[i]<- theta[2]^2*(1-(bstar[i]*stats::dnorm(bstar[i])-
                               astar[i]*stats::dnorm(astar[i]))/0.0001)
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((stats::dnorm(bstar[i])-
                                          stats::dnorm(astar[i]))/0.0001))
      }
    else{
      Bj[i]<- theta[2]^2*(1-(bstar[i]*stats::dnorm(bstar[i])-
                               astar[i]*stats::dnorm(astar[i]))/dinom[i])
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((stats::dnorm(bstar[i])-
                                          stats::dnorm(astar[i]))/dinom[i]))}

  }

  ss<- base::sum(Bj*freq)/base::sum(freq)
  return(ss)
}
################################################################################

### M-step maximization step of the EM algorithm -------------------------------

#' Maximization step of the EM algorithm
#'
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param freq Frequency over the intervals, values in vector.
#' @param theta_init The initial values of the parameters, for mu and sigma.
#'  Vector with two values.
#' @param maxit The maximum number of iteration of the EM algorithm.
#' @param tol1 A number, the stopping criteria for updating mu.
#' @param tol2 A number, the stopping criteria for updating sigma.
#' @return This is the maximization step of the EM algorithm (M-step) that has
#'  defined it using the function E-Step for mean estimate and variance
#'  estimate. Return a list has as arguments "mu_estimate" for the average and
#'  "sigma_estimate" for the variance.
#' @export
#' @examples
#'
#'  library(infgrouped)
#'
#'
#'  output2 <- base::list()
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                        n=50,
#'                        nclass = 10,
#'                        mean = 68,
#'                        sd = 1.80,
#'                        fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  output2 <- infgrouped::em(bl=simdataaaa$simul_data[,1,1],
#'                bu=simdataaaa$simul_data[,2,1],
#'                freq=simdataaaa$simul_data[,3,1],
#'                theta_init=c(67,2),
#'                maxit = 1000,
#'                tol1=1e-3,
#'                tol2=1e-4)
#'  output2


em<- function(bl,bu,freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]

  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    munew<- infgrouped:::mest(theta=cur,bl,bu,freq)
    ssnew<- infgrouped:::ssest(theta=cur,bl,bu,
                  muupdate=infgrouped:::mest(theta=cur,bl,bu,freq) ,freq)

    Mu_new <- munew
    S_new <- base::sqrt(ssnew)
    new_step<- c(Mu_new,S_new)

    if(base::abs(cur[1]-new_step[1])<tol1 &
       base::abs(cur[2]-new_step[2])<tol2){
      flag<-1 ;break
      }
    Mu_cur<- Mu_new
    S_cur<- S_new
  }

  if(!flag){
    warning("Didn't Converge \n")
  }

  updateres<- base::list("mu_estimate" = Mu_cur,
                         "sigma_estimate" = (S_cur)^2)

  return(updateres)
}

################################################################################

### Standard Errors For EM Algorithm -------------------------------------------

#' Standard Errors For EM Algorithm
#'
#' @param thetaupd Output of the EM algorithm function.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param data Matrix of the raw data used in the estimation in the EM
#' algorithm.
#' @return Returns a matrix with the standard errors of the errors of the EM
#' algorithm estimates.
#' @export
#' @examples
#'
#'set.seed(1245)
#'#simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'simdata50m15 <- univ_simul(ncol_matrix=500,
#'                         n=50,
#'                         nclass = 15,
#'                         mean = 68,
#'                         sd = 1.80,
#'                         fr_breaks=c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,
#'                                     73,74.5,76,77.5,79))
#'
#'
#'  output50m15em<- matrix(rep(0,2*500),
#'                         ncol=2)
#'
#'  for(i in 1:500){
#'
#'    output2 <- infgrouped::em(bl=simdata50m15$simul_data[,1,i],
#'                 bu=simdata50m15$simul_data[,2,i],
#'                 freq=simdata50m15$simul_data[,3,i],
#'                 theta_init=c(67,2),
#'                 maxit = 100,
#'                 tol1=1e-3,
#'                 tol2=1e-4)
#'
#'  output50m15em[i,1] <- output2$mu_estimate
#'  output50m15em[i,2] <- output2$sigma_estimate
#'
#'  }
#'
#'output50m15em
#'
#'se_em_50<- matrix(rep(0,2*nrow(output50m15em)),ncol=2)
#'
#'for(i in 1:nrow(output50m15em)){
#'  se_em_50[i,]<- infgrouped::muvarstd(thetaupd=c(output50m15em[i,1],
#'                                     sqrt(output50m15EM[i,2])),
#'                          bl=simdata50m15$simul_data[,1,i],
#'                          bu=simdata50m15$simul_data[,2,i],
#'                          data=simdata50m15$simul_data[,,i])
#'
#'}
#'
#'ci_mu_50<- matrix(rep(0,2*nrow(output50m15em)),
#'                  nrow=nrow(output50m15em),ncol=2)
#'colnames(ci_mu_50)<- c("lower bound","upper bound")
#'
#'ci_sigma2_50<- matrix(rep(0,2*nrow(output50m15em)),
#'                      nrow=nrow(output50m15em),ncol=2)
#'colnames(ci_sigma2_50)<- c("lower bound","upper bound")
#'
#'
#'for(i in 1:nrow(output50m15em)){
#'  ci_mu_50[i,1]<- output50m15em[i,1]-(se_em_50[i,1]*qnorm(0.975))
#'  ci_mu_50[i,2]<- output50m15em[i,1]+(se_em_50[i,1]*qnorm(0.975))
#'
#'}
#'
#'
#'for(i in 1:nrow(output50m15em)){
#'  ci_sigma2_50[i,1]<- output50m15em[i,2]-(se_em_50[i,2]*qnorm(0.975))
#'  ci_sigma2_50[i,2]<- output50m15em[i,2]+(se_em_50[i,2]*qnorm(0.975))
#'
#'}
#'
#'
#'true_mu<- 68
#'true_sigma2<- (1.80)^2
#'
#'
#'empirical_conf_mu_50<- 0
#'empirical_conf_sigma2_50<- 0
#'
#'for(i in 1:nrow(output50m15em)){
#'
#'  if(true_mu >= ci_mu_50[i,1] & true_mu <= ci_mu_50[i,2]){
#'
#'    empirical_conf_mu_50<- empirical_conf_mu_50+1
#'  }
#'
#'  if(true_sigma2>=ci_sigma2_50[i,1] & true_sigma2<=ci_sigma2_50[i,2]){
#'
#'    empirical_conf_sigma2_50<- empirical_conf_sigma2_50+1
#'
#'  }
#'
#'}
#'
#'emp_conf_mu_prop<- c(empirical_conf_mu_50/500)
#'emp_conf_sigma2_prop<- c(empirical_conf_sigma2_50/500)
#'
#'n<- c(50)
#'emp_ci_mu_prop<- emp_conf_mu_prop*100
#'emp_ci_sigma2_prop<- emp_conf_sigma2_prop*100
#'
#'sd_mu <- c(sd(output50m15em[,1]))
#'ave_mu <- c(mean(output50m15em[,1]))
#'
#'se_mu_hat<- c(mean(se_em_50[,1]))
#'sd_var<- c(sd(output50m15em[,2]))
#'
#'
#'ave_var<- c(mean(output50m15em[,2]))
#'se_sigma2_hat<- c(mean(se_em_50[,2]))
#'
#'
#'
#'std_mu_estimates<- cbind(n,ave_mu,sd_mu,se_mu_hat,emp_ci_mu_prop)
#'std_sigma2_estimates<- cbind(n,ave_var,sd_var,se_sigma2_hat,emp_ci_sigma2_prop)
#'
#'std_sigma2_estimates
#'std_mu_estimates
#'
#'


muvarstd<- function(thetaupd,
                    bl,
                    bu,
                    data){
  Wj<- base::rep(0,base::length(bl))
  WjN<- base::rep(0,base::length(bl))

  astar<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))
  Mstdj<- base::rep(0,base::length(bl))
  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-thetaupd[1])/thetaupd[2]
    astar[i]<- (bl[i]-thetaupd[1])/thetaupd[2]

  }
  dinom1<- NULL
  astar[1]<- -1000
  bstar[base::length(bl)]<- 1000
  for(i in 1:base::length(bl)){

    dinom1[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    if(dinom1[i]==0) {
      WjN[i]<- (-thetaupd[2])*
        ((stats::dnorm(bstar[i])-stats::dnorm(astar[i]))/0.0001)
    }
    else {
      WjN[i]<- (-thetaupd[2])*
        ((stats::dnorm(bstar[i])-stats::dnorm(astar[i]))/dinom1[i])
    }
  }

  Ej2<- base::rep(0,base::length(bl))

  dinom<- NULL

  astar[1]<- -1000
  bstar[base::length(bl)]<- 1000

  for(i in 1:base::length(bl)){

    dinom[i]<- (stats::pnorm(bstar[i])-stats::pnorm(astar[i]))
    if(dinom[i]==0) {
      Ej2[i]<- (bstar[i]*
                  stats::dnorm(bstar[i])-astar[i]*
                  stats::dnorm(astar[i]))/0.0001
    }
    else{
      Ej2[i]<- (bstar[i]*
                  stats::dnorm(bstar[i])-astar[i]*
                  stats::dnorm(astar[i]))/dinom[i]
    }
  }
  SMU<- WjN/(thetaupd[2]^2)
  SVAR<- (-1/(2*(thetaupd[2]^2)))*Ej2
  Sj1<- base::rbind(SMU,SVAR)
  IME<- base::array(base::rep(0,2*2*(base::nrow(data))),
                    c(2,2,base::nrow(data)))
  for(i in 1:base::nrow(data)){
    IME[,,i]<- Sj1[,i]%*%t(Sj1[,i])*data[i,3]
  }

  InfME<- base::apply(IME,c(1,2),sum)

  var_EM_est<- base::solve(InfME)
  std_EM_est<- base::sqrt(base::diag(var_EM_est))

  return(std_EM_est)

}

################################################################################






### func univ_simul (gerando simdata2) -----------------------------------------

#' Generating data from a gaussian
#'
#' @param ncol_matrix A number, number of repetitions of the random simulated
#'  samples.
#' @param n A number, number of observations in each sample.
#' @param nclass A number, number of classes on the contingency table.
#' @param mean A number, parameter mu value to generate data.
#' @param sd A number, parameter sigma (standard error) value to generate data.
#' @param fr_breaks A vector, vector containing the values of the referenced
#'  as counts.
#' @return Returns a list with two objects, the first "simul_data" returns
#'  a matrix in which the first column the lower limits of the intervals,
#'  in the second column the values of the upper limits of the intervals and
#'  the third column the counts of the observations. The second object "med"
#'  returns a matrix, in which the first column shows the sampled value and
#'  the second column the observation count.
#' @export
#' @examples
#'
#'  library(infgrouped)
#'  library(stats)
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                        n=50,
#'                        nclass = 10,
#'                        mean = 68,
#'                        sd = 1.80,
#'                        fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  simdataaaa
#'


univ_simul <- function(ncol_matrix=30,
                       n=50,
                       nclass = 10,
                       mean = 68,
                       sd = 1.80,
                       fr_breaks=c(62,64,66,68,70,72,74,76,78)){

  sim2<- base::matrix(rep(0,n*ncol_matrix),
                ncol=ncol_matrix)
  for(i in 1:base::ncol(sim2)){
    sim2[,i]<- stats::rnorm(n = n,
                            mean = mean,
                            sd = sd)
  }

  ### ### ### ###
  Fr<- base::matrix(0,
                    ncol= ncol_matrix,
                    nrow = length(fr_breaks)+1)


  for(i in 1:base::ncol(sim2)){
    Fr[,i]<- base::table(base::cut(sim2[,i],
                       breaks=c(-Inf,fr_breaks,Inf)))

  }
  ### ### ### ###

  simdata2<- base::array(base::rep(0,10*3*ncol_matrix),
                   c(nclass,3,ncol_matrix))

  med2<- base::array(base::rep(0,10*2*ncol_matrix),
               c(nclass,2,ncol_matrix))

  for(i in 1:base::ncol(sim2)){
    simdata2[,1,i]<- c(-Inf,fr_breaks)
    simdata2[,2,i]<- c(fr_breaks,Inf)
    simdata2[,3,i]<- Fr[,i]

    med2[,1,i]<- (simdata2[,1,i]+simdata2[,2,i])/2
    med2[1,1,i]<- base::min(fr_breaks)-1
    med2[nclass,1,i]<- base::max(fr_breaks)+1
    med2[,2,i]<- Fr[,i]
  }

  final_list <- base::list("simul_data" = simdata2,
                    "med" = med2)

  return(final_list)

}
################################################################################


### Galton Data ----------------------------------------------------------------

#' Galton Data
#'
#' @return Returns a list with three formats of Galton data, the first is
#' the parents data in the frequency table, the second is the children data
#' in the frequency table and the third is the raw format of the data.
#' @export
#' @examples
#'
#'
#' galton_data_formats = galton_data()
#'
#'


library(HistData)
library(utils)


galton_data<- function(){

  library(HistData)

  utils::data(Galton)
  data_galton_parent <- Galton
  data_galton_parent <- data_galton_parent$parent

  freq <- base::table(base::cut(data_galton_parent,
                                breaks=c(-Inf,
                                         unique(data_galton_parent),
                                         Inf)))
  freq <- as.data.frame(freq)

  intervalo <- base::unique(data_galton_parent)
  intervalo <- intervalo[base::order(intervalo)]

  breaks <- c(-Inf,
              intervalo,
              Inf)

  TL <- breaks[-13]
  TU <- breaks[-1]

  data_galton_parent <- base::data.frame("TL_parent" = TL,
                                         "TU_parent" = TU,
                                         "freq" = freq$Freq)


  ########## ########## ##########

  data_galton_children <- Galton
  data_galton_children <- data_galton_children$child

  freq <- base::table(base::cut(
    data_galton_children,breaks=c(-Inf,
                                  unique(data_galton_children),
                                  Inf)))

  freq <- as.data.frame(freq)


  intervalo <- base::unique(data_galton_children)
  intervalo <- intervalo[base::order(intervalo)]


  breaks <- c(-Inf,
              intervalo,
              Inf)

  TL <- breaks[-16]
  TU <- breaks[-1]

  data_galton_children <- base::data.frame("TL_parent" = TL,
                                           "TU_parent" = TU,
                                           "freq" = freq$Freq)


  final_list <- base::list("data_galton_children" = data_galton_children,
                           "data_galton_parent" = data_galton_parent,
                           "data_galton" = Galton)

  detach("package:HistData", unload = TRUE)
  return(final_list)

}

################################################################################






### E-Step for Mu estimation: Simulating Z's  ----------------------------------
###Calculation of updates for Mu###

#' E-Step for Mu estimation, Simulating Z's
#'
#' @param theta The current state of the parameters, vector containing
#'  mean and sd.
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies.
#' @return Simulate data about each specific interval and assign it to the
#'  simulate matrix.
#' @examples
#'
#'  library(truncnorm)
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  zm_matriz <- infgrouped:::zmcem(theta=c(67,2),
#'  data = simdataaaa$simul_data[,,1])
#'  zm_matriz

zmcem<- function(theta,data){

  k<- 1000
  sim<- base::matrix(rep(0,k*base::nrow(data)),ncol=base::nrow(data))

  for(i in 1:base::nrow(data)) {
    sim[,i]<- truncnorm::rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2])
  }

  return(sim)
}

################################################################################

### E-Step for Mu --------------------------------------------------------------

#' E-Step estimate the parameter of mu.
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies. The current state of the parameters, vector containing
#'  mean and sd.
#' @param simz The matrix, which is the samples simulated over the intervals,
#'  of the zmcem function.
#' @return Return mu which are the estimates of mean in E-step.
#' @examples
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  zm_matriz <- infgrouped::zmcem(theta=c(67,2),
#'                                  data = simdataaaa$simul_data[,,1])
#'
#'  mu_estimate <- infgrouped:::mumcem(data = simdataaaa$simul_data[,,1],
#'                            simz = zm_matriz)
#'  mu_estimate


mumcem<- function(data,simz){

  n<- base::sum(data[,3])
  Z<- base::colMeans(simz)
  numerator<- base::rep(0,base::nrow(data))
  for (i in 1:base::nrow(data)) {
    numerator[i]<- data[i,3]*Z[i]
  }

  MuN<- (1/n)*base::sum(numerator)
  return(MuN)
}

################################################################################

### Calculate estimate of sigma ------------------------------------------------

#' E-Step estimate the parameter of sigma.
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies. The current state of the parameters, vector containing
#'  mean and sd.
#' @param simzz The matrix, which is the samples simulated over the intervals,
#'  of the zmcem function.
#' @param mupd Is the updated estimates of mu.
#' @return Return sigma which are the estimates of variance in E-step.
#' @examples
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  zm_matriz <- infgrouped:::zmcem(theta=c(67,2),
#'                                  data = simdataaaa$simul_data[,,1])
#'
#'  mu_estimate <- infgrouped:::mumcem(data = simdataaaa$simul_data[,,1],
#'                                     simz = zm_matriz)
#'
#'  sigma_estimate <- infgrouped:::sigmamcem(data = simdataaaa$simul_data[,,1],
#'                             simz = zm_matriz,
#'                             mupd = mu_estimate)
#'  sigma_estimate



sigmamcem<- function(data,simzz,mupd){

  n<- base::sum(data[,3])
  ZZ<- simzz
  NewZ<- (ZZ-mupd)^2
  SZNEW<- base::colMeans(NewZ)
  numerator<- base::rep(0,nrow(data))

  for (i in 1:base::nrow(data)){
    numerator[i]<- data[i,3]*SZNEW[i]
  }

  sigmaNN<- (1/n)*base::sum(numerator)
  sig<- base::sqrt(sigmaNN)

  return(sig)
}

################################################################################

### MONTE CARLO EM -------------------------------------------------------------

#' This is the maximization step of the mcem algorithm (M-step).
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies.
#' @param theta_init The initial values of the parameters, for mu and sigma.
#'  Vector with two values.
#' @param maxit The maximum number of iteration of the EM algorithm.
#' @param tol1 A number, the stopping criteria for updating mu.
#' @param tol2 A number, the stopping criteria for updating sigma.
#' @return This is the maximization step of the mcem algorithm (M-step) that has
#'  defined it using the function E-Step for mean estimate and variance
#'  estimate. Returns the estimates for mean (mu) and sigma (variance).
#' @export
#' @examples
#'
#'  library(infgrouped)
#'
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  outputmcem2 <- base::list()
#'
#'  outputmcem2[[1]] <- infgrouped::mcem(data=simdataaaa$simul_data[,,1],
#'                            theta_init=c(67,2),
#'                            maxit = 1000,tol1=1e-2,tol2=1e-3)
#'
#'  outputmcem2


mcem<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){

  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]
  iter<- base::rep(0,maxit)
  Svec<- base::rep(0,maxit)
  Mvec<- base::rep(0,maxit)
  mydat<- data

  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    munew<- infgrouped:::mumcem(data=mydat,
                            simz=infgrouped:::zmcem(theta=cur,
                                                data=mydat))

    Snew<- infgrouped:::sigmamcem(data=mydat,simzz=zmcem(theta=cur,data=mydat),
                     mupd=infgrouped:::mumcem(data=mydat,
                                          simz=infgrouped:::zmcem(theta=cur,
                                                              data=mydat)))

    Mu_new<- munew
    S_new<- Snew
    new_step<- c(Mu_new,S_new)

    if(base::abs(cur[1]-new_step[1])<tol1 & base::abs(cur[2]-new_step[2])<tol2){
      flag<-1 ;break
    }

    Mu_cur<- Mu_new
    S_cur<- S_new
    iter[i]<- i
    Svec[i]<- S_new
    Mvec[i]<- Mu_new

  }

  if(!flag) warning("Didn't Converge \n")
  update <- base::list("mu_estimate" = Mu_cur,
                       "sigma_estimate" = (S_cur)^2)

  return(update)
}
################################################################################


### Standard Errors For MCEM Algorithm -----------------------------------------

#' Standard Errors For EM Algorithm
#'
#' @param thetaupd Output of the EM algorithm function.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param data Matrix of the raw data used in the estimation in the EM
#' algorithm.
#' @return Returns a matrix with the standard errors of the errors of the MCEM
#' algorithm estimates. First column is the standard deviation and the
#' second column is the variance.
#' @export
#' @examples
#'
#'
#'outputmcem50m15<- matrix(rep(0,2*500),ncol=2)
#'colnames(outputmcem50m15)<- c("mean","var")
#'
#'set.seed(1245)
#'simdata50m15 <- univ_simul(ncol_matrix=500,
#'                        n=50,
#'                        nclass = 15,
#'                        mean = 68,
#'                        sd = 1.80,
#'                        fr_breaks=c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,
#'                                    73,74.5,76,77.5,79))
#'
#'
#'for(i in 1:500){
#'  outputmcem2 <- infgrouped::mcem(data=simdata50m15$simul_data[,,i],
#'                             theta_init=c(67,2),
#'                             maxit = 1000,
#'                             tol1=1e-2,
#'                             tol2=1e-3)
#'
#'  outputmcem50m15[i,1] <- outputmcem2$mu_estimate
#'  outputmcem50m15[i,2] <- outputmcem2$sigma_estimate
#'
#'
#'}
#'
#'
#'theta50<-matrix(rep(0,2*nrow(output50m15em)),ncol=2)
#'for(i in 1:nrow(output50m15em)){
#'  theta50[i,]<- c(output50m15em[i,1],
#'                  sqrt(output50m15em[i,2]))
#'}
#'
#'
#'se_mcem50m15<- matrix(rep(0,2*nrow(output50m15em)),ncol=2)
#'for(i in 1:nrow(output50m15em)){
#'  se_mcem50m15[i,]<-  zsimmcem(theta50[i,],
#'                               data=simdata50m15$simul_data[,,i])
#'}
#'
#'# remove cases with not converge
#'errors_lines_c1 = which(is.na(se_mcem50m15[,1]))
#'errors_lines_c2 = which(is.na(se_mcem50m15[,2]))
#'errors_lines = c(errors_lines_c1,errors_lines_c2)
#'se_mcem50m15 = se_mcem50m15[-errors_lines, ]
#'outputmcem50m15 = outputmcem50m15[-errors_lines, ]
#'
#'trueval<- c(68,(1.8)^2)
#'
#'ci_mu_50_mcem<- matrix(rep(0,2*nrow(outputmcem50m15)),
#'                       nrow=nrow(outputmcem50m15),
#'                       ncol=2)
#'
#'colnames(ci_mu_50_mcem)<- c("lower bound",
#'                            "upper bound")
#'
#'ci_sigma2_50_mcem<- matrix(rep(0,2*nrow(outputmcem50m15)),
#'                           nrow=nrow(outputmcem50m15),
#'                           ncol=2)
#'
#'colnames(ci_sigma2_50_mcem)<- c("lower bound",
#'                                "upper bound")
#'
#'
#'for(i in 1:nrow(outputmcem50m15)){
#'
#'  ci_mu_50_mcem[i,1]<- outputmcem50m15[i,1]-(se_mcem50m15[i,1]*qnorm(0.975))
#'  ci_mu_50_mcem[i,2]<- outputmcem50m15[i,1]+(se_mcem50m15[i,1]*qnorm(0.975))
#'
#'}
#'
#'for(i in 1:nrow(outputmcem50m15)){
#'
#'  ci_sigma2_50_mcem[i,1]<- outputmcem50m15[i,2]-(se_mcem50m15[i,2]*qnorm(0.975))
#'  ci_sigma2_50_mcem[i,2]<- outputmcem50m15[i,2]+(se_mcem50m15[i,2]*qnorm(0.975))
#'
#'}
#'
#'true_mu<- 68
#'true_sigma2<- (1.8)^2
#'
#'
#'empirical_conf_mu_50_mcem<- 0
#'empirical_conf_sigma2_50_mcem<- 0
#'
#'for(i in 1:nrow(outputmcem50m15)){
#'
#'  if(true_mu >= ci_mu_50_mcem[i,1] &
#'     true_mu <= ci_mu_50_mcem[i,2]){
#'
#'    empirical_conf_mu_50_mcem<- empirical_conf_mu_50_mcem+1
#'
#'  }
#'
#'  if(true_sigma2 >= ci_sigma2_50_mcem[i,1] &
#'     true_sigma2 <= ci_sigma2_50_mcem[i,2]){
#'
#'    empirical_conf_sigma2_50_mcem<- empirical_conf_sigma2_50_mcem+1
#'
#'  }
#'
#'}
#'
#'
#'conf_mu_prop_mcem<- c(empirical_conf_mu_50_mcem/500)
#'conf_sigma2_prop_mcem<- c(empirical_conf_sigma2_50_mcem/500)
#'
#'
#'n<- c(50)
#'ci_mu_prop_mcem<- conf_mu_prop_mcem*100
#'ci_sigma2_prop_mcem<- conf_sigma2_prop_mcem*100
#'
#'
#'sd_mu_mcem<- c(sd(outputmcem50m15[,1]))
#'ave_mu_mcem<- c(mean(outputmcem50m15[,1]))
#'
#'
#'se_mu_hat_mcem<- c(mean(se_mcem50m15[,1]))
#'sd_var_mcem<- c(sd(outputmcem50m15[,2]))
#'
#'
#'ave_var_mcem<- c(mean(outputmcem50m15[,2]))
#'se_sigma2_hat_mcem<-  c(mean(se_mcem50m15[,2]))
#'
#'
#'
#'std_mu_estimates_mcem<- cbind(n,ave_mu_mcem,sd_mu_mcem,
#'                              se_mu_hat_mcem,ci_mu_prop_mcem)
#'std_sigma2_estimates_mcem<- cbind(n,ave_var_mcem,sd_var_mcem,
#'                                  se_sigma2_hat_mcem,ci_sigma2_prop_mcem)
#'
#'std_mu_estimates_mcem
#'std_sigma2_estimates_mcem
#'




zsimmcem<- function(theta,
                    data){

  k<- 1000
  sim<- base::matrix(rep(0, k*base::nrow(data)),
                     ncol=base::nrow(data))

  for(i in 1:base::nrow(data)) {
    sim[,i]<- truncnorm::rtruncnorm(k,
                                    a=data[i,1],
                                    b=data[i,2],
                                    mean=theta[1],
                                    sd=theta[2]) }

  wj<- base::matrix(base::rep(0,k*nrow(data)),
                    ncol=base::nrow(data))
  wj2<- base::matrix(base::rep(0,k*nrow(data)),
                     ncol=base::nrow(data))

  for(i in 1:base::nrow(data)){

    wj[,i]<- sim[,i]-theta[1]
    wj2[,i]<- (sim[,i]-theta[1])^2

  }

  d2<- base::array(base::rep(0,2*2*base::nrow(data)),
                   c(2,2,base::nrow(data)))
  for (i in 1:base::nrow(data)){

    d2[1,1,i]<- 1/theta[2]^2
    d2[1,2,i]<- d2[2,1,i]<- (1/theta[2]^4)*base::mean(wj[,i])
    d2[2,2,i]<- (-1/(2*theta[2]^4))+((1/theta[2]^6)*base::mean(wj2[,i]))

  }

  d1mu<- wj/(theta[2]^2)
  d1sigma2<- (-1/(2*theta[2]^2))+(wj2/(2*theta[2]^4))

  dd1<- base::array(base::rep(0,2*k*base::nrow(data)),c(k,2,base::nrow(data)))
  for(i in 1:base::nrow(data)){

    dd1[,1,i]<- d1mu[,i]
    dd1[,2,i]<- d1sigma2[,i]

  }

  d1<- base::array(base::rep(0,2*1*base::nrow(data)),
                   c(2,1,base::nrow(data)))
  for(i in 1:base::nrow(data)){

    d1[1,1,i]<- ((base::mean(wj[,i]))/(theta[2]^2))
    d1[2,1,i]<- (-1/(2*theta[2]^2))+((base::mean(wj2[,i]))/(2*theta[2]^4))

  }

  diff1<- base::array(base::rep(0,2*k*base::nrow(data)),
                      c(k,2,base::nrow(data)))
  for(i in 1:base::nrow(data)){

    diff1[,1,i]<- dd1[,1,i]-d1[1,1,i]
    diff1[,2,i]<- dd1[,2,i]-d1[2,1,i]

  }

  diff1_2<- base::array(base::rep(0,2*2*k*base::nrow(data)),
                        c(2,2,k,base::nrow(data)))
  for(i in 1:base::nrow(data)){

    for(j in 1:k){

      diff1_2[,,j,i]<- diff1[j,,i]%*%t(diff1[j,,i])

    }
  }

  mydif<- base::apply(diff1_2,c(1,2,4),mean)
  mydifmat<- base::array(rep(0,2*2*base::nrow(data)),
                         c(2,2,base::nrow(data)))

  for(i in 1:base::nrow(data)){

    mydifmat[,,i]<- (d2[,,i]-mydif[,,i])*data[i,base::ncol(data)]

  }

  myinf_init<- base::apply(mydifmat,c(1,2),sum)
  myinf_final<- base::solve(myinf_init)

  se_MCEM_est<- base::sqrt(base::diag(myinf_final))
  base::names(se_MCEM_est)<- c("std_mu","std_sigma2")
  return(se_MCEM_est)

}


################################################################################





### LOG L FUNCTION -------------------------------------------------------------

#' Likelihood Estimation.
#'
#' @param tl Lower end of the intervals.
#' @param freq Frequency on each interval.
#' @param theta The arguments including mu and sigma.
#' @return Return the result of the log likelihood for the mean (mu) and sigma
#'  (standard deviation) estimates.
#' @examples
#'
#'  library(stats)
#'
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'
#'  tl <- simdataaaa$simul_data[,1,1]
#'  freq <- simdataaaa$simul_data[,3,1]
#'  res <- stats::optim(c((67/2),(1/2)),fn=Logll,tl=tl,freq=freq,
#'  method="L-BFGS-B",lower=c(0.02,0.2),upper=c(180,2))
#'
#'  estimate<- res$par
#'  mu_estimate<- base::round(res$par[1]/res$par[2],4)
#'  sigma_estimate<- base::round(1/res$par[2],4)
#'  sigma_squared_estimate = sigma_estimate^2
#'
#'  mu_estimate
#'  sigma_squared_estimate


Logll <- function(tl,freq,theta){
  m<- base::length(tl)

  if( (stats::pnorm(theta[2]*tl[2]-theta[1])) < 1e-16  ){
    a <- -1e+6
  } #end if

  else{
    a <- freq[1]*base::log(stats::pnorm(theta[2]*tl[2]- theta[1]))
  } #end else

  if( (1-stats::pnorm(theta[2]*tl[m]-theta[1])) < 1e-16  ) {
    b <- -1e+6
  }#end if

  else{
    b <- freq[m]*base::log(1-stats::pnorm(theta[2]*tl[m]-theta[1]))
  }#end else

  c<-0
  for(i in 2:(m-1)){

    if ( (stats::pnorm(theta[2]*tl[i+1]-theta[1]) -
          stats::pnorm(theta[2]*tl[i]-theta[1])) < 1e-16 ){
      c <- c -1e+6
    }#end if

    else{
        c <- c + freq[i]*
          (base::log( stats::pnorm(theta[2]*tl[i+1]-theta[1]) -
                  stats::pnorm(theta[2]*tl[i]-theta[1])))
    }#end else
  }#end for

  L <- -(a+b+c)

  return(L)
}
################################################################################


### exact_mle FUNCTION -------------------------------------------------------------

#' Maximizing the Likelihood Estimation 1.
#'
#' @param tl Lower end of the intervals.
#' @param freq Frequency on each interval.
#' @param intial_mu Initial numeric value for mu (mean).
#' @param intial_sigma Initial numeric value for sigma.
#' @return Return the result of the log likelihood for the mean (mu) and sigma
#'  (variance) estimates.
#' @export
#' @examples
#'  simdataaaa <- infgrouped::univ_simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'
#'  tl<- simdataaaa$simul_data[,1,1]
#'  freq<- simdataaaa$simul_data[,3,1]
#'
#'  estimates_simul <- infgrouped::exact_mle(tl=tl,
#'                              freq=freq,
#'                              intial_mu = (67/2),
#'                              intial_sigma = (1/2))
#'
#'  estimates_simul

exact_mle <- function(tl = c(-Inf,62,64,66,68,70,72,74,76,78),
                      freq = c(0,1,9,16,15,8,1,0,0,0),
                      intial_mu = (67/2),
                      intial_sigma = (1/2)){

  res <- stats::optim(c(intial_mu,intial_sigma),
                      fn=Logll,
                      tl=tl,
                      freq=freq,
                      method="L-BFGS-B",
                      lower=c(0.02,0.2),
                      upper=c(180,2))

  mu_estimate<- base::round(res$par[1]/res$par[2],4)
  sigma_estimate<- base::round(1/res$par[2],4)
  sigma_squared_estimate = sigma_estimate^2

  final_list <- base::list("mu_estimate" = mu_estimate,
                           "sigma_estimate" = sigma_squared_estimate)

  return(final_list)

}
################################################################################









