#' Using imported population and control parameters, run a population viability analysis (PVA) on the target species.
#' @import ggplot2
#' @import tidyr
#' @param params A list of initialized population and control parameters to inform the PVA. Parameters should be provided in the form of a named list. We suggest filling in a parameter template, which can be created and loaded using the \code{pva_template()} and \code{load_pva_parameters()} functions.
#' @param custom_inits (Optional, invoked by the `rank_uncertainty` function) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pvA_params}. Should be a named list_ Can be be outputs of the \code{init()} function from \code{PVAInvas}.
#' @param sens_percent (Optional, invoked by the `rank_uncertainty` function) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param create_plot (Optional) Should a ggplot heatmap object also be provided? Default: false. Will return a list with the final entry being the created plot.
#' @return pva
#' * A named list of PVA outputs, including calculated parameters (from init())
#'   and outputs of the PVA_
#' * Calculated outputs from init():
#'  - `phie`: unfished eggs per recruit at equilibrium,
#'  - `R_A`: stage-independent maximum survival (alpha parameter of Beverton-Holt recruitment),
#'  - `R_B`: stage-independent carrying capacity (beta parameter of Beverton-Holt recruitment),
#'  - `A_s`: stanza-specific maximum survival (alpha parameter of Beverton-Holt recruitment),
#'  - `B_s`: stanza-specifc carrying capacity (beta parameter of Beverton-Holt recruitment).
#'
#' * Objects created from PVA simulations
#'  - `Nt`: 3-dimensional abundance array (dimensions: time-steps, ages, simulations),
#'  - `Et`: matrix of eggs for each timestep for each year in each simulation,
#'  - `nest`: !!!!!!!!!!!! estimated numbers for each year in each simulation,
#'  - `Vfin`: vector of abundance in the final year across simulations,
#'  - `p_extinct`: vector of proportion of erdicated in each time step,
#'  - `p_extinct_50`: proportion of simulations where the population is eradicated by the 50th time step,
#'  - `p_extinct_100`: proportion of simulations where the population is eradicated by the 100th time step,
#'  - `p_extinct_200`: proportion of simulations where the population is eradicated by the 200th time step,
#'  - `t_extinct`: minimum number of timesteps needed for eradiction,
#'  - `yext_seq`:
#'  - `cost_1`: annual cost of sampling,
#'  - `cost_T`: total cost of sampling (up to the end of simulation or until 100% eradication),
#'  - `NPV`: net present value of sampling (taking into account intergenerational discounting),
#'  - `E_NPV`: expected mean present value (mean of NPV)
#'  - `NT`: abundance in the final time-step (reported as 5th percentile, mean, and 95 percentile of distributions),
#'  - `runtime`: time to execute the PVA,
#'  - `plot`: (Optional) if `create_plot=TRUE`, `plot` is returned as a ggplot object with multiple components. Created by the function `vwReg2.R`.
#'
#' @examples
#' # Run a simple PVA, no custom values or sensitivity testing.
#' pva(pva_params = inputParameterList)
#' @export

PVA <- function(params, custom_inits = NULL, sens_percent = NULL,
  sens_params = NULL, create_plot = FALSE, set_plot_y = NULL, quiet = FALSE){
  start <- Sys.time()
  if(!quiet){
    message("Calculating population projections ...\n")
  }
  if(!is.null(custom_inits)){
    inits <- custom_inits$initialized_params
  } else {
    inits <- init(params, quiet = quiet)$initialized_params
  }
  AR <- inits$AR
  A <- inits$A
  dt <- inits$dt
  R0 <- inits$R0_vec #!# Is this correct?
  age <- inits$age
  la <- inits$la
  wa <- inits$wa
  spn <- inits$spn
  fec <- inits$fec
  mat <- inits$mat
  Sa <- inits$Sa
  lx <- inits$lx
  sel <- inits$sel
  phie <- inits$phie
  R_A <- inits$R_A
  R_B <- inits$R_B
  A_s <- inits$A_s
  B_s <- inits$B_s
  can_a <- inits$can_a
  can_b <- inits$can_b
  M1 <- inits$M1
  Sa_M <- inits$Sa_M
  Ct <- inits$Ct
  Nt <- inits$Nt
  Et <- inits$Et
  nest <- inits$nest
  nT <- inits$nT                   # number of time-steps
  nS <- params$nS                   # number of pre-recruit stanzas
  n_sim <- params$n_sim             # number of simulations
  samp_A <- params$samp_A
  # samp_A <- vector()               # time-step within a year that each gear is fished
  r <- params$r                     # discounting rate = future generation discount factor
  G <- params$G                     # generation time (years)
  U_R <- params$U_R # vector()                  # proportion of pre-recruited animals removed per gear
  U_A <- params$U_A # vector()                  # proportion of recruited animals removed per gear
  q_R <- params$q_R # vector()                  # catchability of gear used to remove fish from each pre-recruit stanza
  q_A <- params$q_A # vector()                  # catchability of gear used to remove recruited animals
  n_gear <- params$n_gear           # number of capture gears applied to recruited animals
  t_start_R <- params$t_start_R # vector()            # time-step when sampling begins for each pre-recruit stanza
  t_start_A <- params$t_start_A # vector()            # time-step when sampling begins for each gear used on recruited animals
  E_R <- params$E_R # vector()                  # effort per time-step used to remove fish from each pre-recruit stanza
  E_A <- params$E_A # vector()                  # effort per time-step used to remove recruited animals
  v_a <- params$v_a # vector()                  # Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)
  v_b <- params$v_b # vector()                  # Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
  v_c <- params$v_c # vector()                  # Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)
  v_d <- params$v_d # vector()                  # Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
  C_f_R <- params$C_f_R # vector()
  C_f_A <- params$C_f_A # vector()
  C_E_R <- params$C_E_R # vector()
  C_E_A <- params$C_E_A # vector()
  reck <- params$reck             # recruitment compensation ratio
  p_can <- params$p_can           # proportion of recruit mortality at equilibrium due to cannibalism
  K <- params$K                   # von Bertalanffy metabolic parameter
  afec <- params$afec             # slope of fecundity-weight relationship
  Wmat <- params$Wmat             # weight at maturity
  t_spn <- params$t_spn           # range of time of year when spawning occurs
  V1 <- params$V1                 # initial vulnerable abundance (used to create initial population)
  cann_a <- params$cann_a         # age at which cannibalism on pre-recruits begins
  bet <- params$bet               # rate at which invasives disperse with abundance (between 0 (none) and greater)
  sd_S <- params$sd_S             # standard deviation of environmental effect on survival
  # Apply this when testing sensitivity to biological parameters
  if(!is.null(sens_params)){
    # Skip modification of parameters if these are already in inits
    if(!(sens_params %in% names(inits) | sens_params %in% c("Bs","Ms"))){
      stopifnot(!is.null(sens_percent))
      o.value <- get(sens_params)
      neW_val <- o.value*sens_percent
      assign(paste0(sens_params), neW_val)
    }
  }
  U_R <- pmin(U_R,0.999999999)
  U_A <- pmin(U_A,0.999999999)
  q_R <- -log(1-U_R)
  q_A <- -log(1-U_A)*V1^bet

  W_st <- rnorm(nT*n_sim,0,sd_S)
  anom <- matrix(data=exp(W_st),nrow=nT,ncol=n_sim)              # survival anomolies for each stanza
  go_R <- matrix(rep(0,nT*nS),nrow=nS)
  for(i in 1:nS){
    go_R[i,t_start_R[i]:nT] <- 1
  }
  go_A <- matrix(rep(0,nT*n_gear),nrow=n_gear)
  for(i in 1:n_gear){
    go_A[i,t_start_A[i]:nT] <- 1
  }

  # dynamics for subsequent years
  for(t in 2:nT){
    N_st <- matrix(nrow=nS+1,ncol=n_sim)    # numbers surviving through each stanza in a year
    Ct_st <- matrix(nrow=nS,ncol=n_sim)
    N_st[1,] <- Et[t-1,]
    V <- colSums(Nt[t-1,(cann_a/dt):A,])
    M0_t <- sweep(can_b,MARGIN=2,V,'*')
    M0_t <- sweep(M0_t,MARGIN=1,can_a,'+')  # density independent variable
    B_st <- M1/M0_t*(1-exp(-M0_t))              # density dependent variable
    for(st in 1:nS){
      Ct_st[st,] <- N_st[st,]*(1-exp(-E_R[st]*q_R[st]*go_R[st,t]))
      N_st[st+1,] <- rbinom(n_sim,round(pmax(0,N_st[st,]),0),
                            pmin(exp(-M0_t[st,]-E_R[st]*q_R[st]*go_R[st,t])*anom[t,]/
                                   (1+B_st[st,]*N_st[st,]),1))
    }
    Ct[t,1:nS] <- rowSums(Ct_st)
    Nt[t,AR,] <- N_st[nS+1,]                # numbers surviving through all stanzas become recruits to the population

    Ct_A <- array(dim=c(n_gear,A-AR+1,n_sim))
    Ft_A <- list()
    Vt <- colSums(Nt[t-1,AR:A,])
    for(i in 1:n_gear){
      q_N <- q_A[i]*pmax(1,Vt)^(-bet)
      if(t*dt-as.integer(t*dt)==samp_A[i]){
        Ft_A[[i]] <- go_A[i,t]*(E_A[i]*sel[i,] %o% q_N)
      } else {
        Ft_A[[i]] <- go_A[i,t]*matrix(rep(0,(A-AR+1)*n_sim),nrow=A-AR+1,ncol=n_sim)
      }
    }
    Ft <- Reduce('+',Ft_A)
    Zt <- Ft - log(Sa_M)
    for(i in 1:n_gear){
      Ct[t,nS+i] <- sum(Nt[t-1,AR:A,] * Ft_A[[i]]/Zt * ( 1 - exp( -Zt)))
    }
    # Get through the loop?
#      if(t<16) message(t," F\n",Ft[,1],"\nZt",Zt[,1],"\nSa",exp(-Zt[1:(A-AR),1]),"\n")
    Nt[t,(AR+1):A,] <- matrix(rbinom((A-AR)*n_sim,Nt[t-1,AR:(A-1),],exp(-Zt[1:(A-AR),])),ncol=n_sim)
    pairs <- matrix(as.integer(Nt[t,AR:A,]/2),ncol=n_sim)
    Et[t,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
    nest[t,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
  }
  # probabilities of extirpation
  p_extinct_50 <- NA
  p_extinct_100 <- NA
  p_extinct_200 <- NA
  stps <- as.integer(nT/c(4,2,1))
  p_extinct <- sapply(1:nT, function(ii)
    length(which(colSums(Nt[ii,,],na.rm=TRUE)==0))/n_sim) ## calculate p_extinct by time-step
  if( nT >= stps[1] ) p_extinct_50 <- p_extinct[stps[1]]     # probability of extinction in 50 time-steps
  if( nT >= stps[2] ) p_extinct_100 <- p_extinct[stps[2]]  # probability of extinction in 100 time-steps
  if( nT >= stps[3] ) p_extinct_200 <- p_extinct[stps[3]]  # probability of extinction in 200 time-steps
  t_extinct <- min(which(p_extinct==1),nT)
  nyrs <- nT*dt

  Nt[is.na(Nt)]<-0
  y_extinct <- sapply(seq(1/dt,nT,1/dt),function(ii)
    length(which(colSums(Nt[ii,,])==0)))
  y_extinct[2:nyrs] <- y_extinct[2:nyrs]-(y_extinct[1:(nyrs-1)])
  sum_extinct <- sum(y_extinct)
  y_extinct[nyrs] <- n_sim-sum_extinct
  if(is.na(sum_extinct)) y_extinct[nyrs]<-n_sim
  if(sum(y_extinct)<n_sim)y_extinct[nyrs]<-y_extinct[nyrs]+n_sim-sum(y_extinct)
  yext_seq <- rep(1:nyrs,y_extinct)
  y_extinct <- y_extinct / n_sim
  # numbers remaining
  NT <- quantile(colSums(Nt[nT,,]),probs=c(0.025,0.5,0.975))
  # costs of removal
  x_R <- rep(0,length(q_R))
  x_A <- rep(0,length(q_A))
  x_R[which(E_R>0)] <- 1
  x_A[which(E_A>0)] <- 1
  cost_1 <- sum(E_R*C_E_R) + sum(C_f_R[x_R]) +
    sum(E_A*C_E_A) + sum(C_f_A[x_A])
  y_ext <- t_extinct*dt
  cost_T <- mean(cost_1*(1:nyrs)*y_extinct)
  t_st <- (t_start_R[1]-1)*dt+1
  tseq <- (t_start_R[1]*dt):y_ext
  d <- 1/(1+r)
  df <- d
  NPV <- vector()
  for(i in 1:n_sim){
    NPV[i] <- cost_1*sum(d^(t_st:yext_seq[i])+(df^(t_st:yext_seq[i]))/G)
  }
  E_NPV <- mean(NPV)

  runtime <- Sys.time()-start
  # message(runtime)
  out <- list()
  out$initialized_params <- inits
  out$phie <- phie
  out$R_A <- R_A
  out$R_B <- R_B
  out$A_s <- A_s
  out$B_s <- B_s
  out$Nt <- Nt
  out$Et <- Et
  out$nest <- nest
  out$Vfin <- (Nt[nT,AR:A,])
  out$p_extinct <- p_extinct
  out$p_extinct_50 <- p_extinct_50
  out$p_extinct_100 <- p_extinct_100
  out$p_extinct_200 <- p_extinct_200
  out$t_extinct <- t_extinct
  out$yext_seq <- yext_seq
  out$cost_1 <- cost_1
  out$cost_T <- cost_T
  out$NPV <- NPV
  out$E_NPV <- E_NPV
  out$NT <- NT
  out$runtime <- runtime
  if(create_plot==T){
    Na <- apply(Nt,MARGIN=c(1,3),sum,na.rm=TRUE)
    Na <- Na[1:inits$nT,]
    data <- data.frame(Na)
    names(data) <- rep("y",n_sim)
    data$x <- (1:nT)*dt
    # Pulls in vwReg2
    out$plot <- vwReg2(data=data,input=params,set_ymax=set_plot_y)
  }
  return(out)
  gc()
}
