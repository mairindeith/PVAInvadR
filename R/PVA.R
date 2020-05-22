#' Using imported population and control parameters, run a population viability analysis (PVA) on the target species.
#'
#' @param params A list of initialized population and control parameters to inform the PVA. Parameters should be provided in the form of a named list. We suggest filling in a parameter template, which can be created and loaded using the \code{pva_template()} and \code{load_pva_parameters()} functions.
#' @param custom.inits (Optional, invoked by the `rankUncertainty` function) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pva.params}. Should be a named list. Can be be outputs of the \code{init()} function from \code{PVAInvas}.
#' @param sens.pcent (Optional, invoked by the `rankUncertainty` function) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param create.plot (Optional) Should a ggplot heatmap object also be provided? Default: false. Will return a list with the final entry being the created plot.
#' @return pva
#' * A named list of PVA outputs, including calculated parameters (from init())
#'   and outputs of the PVA.
#' * Calculated outputs from init():
#'  - `phie`: unfished eggs per recruit at equilibrium,
#'  - `R.A`: stage-independent maximum survival (alpha parameter of Beverton-Holt recruitment),
#'  - `R.B`: stage-independent carrying capacity (beta parameter of Beverton-Holt recruitment),
#'  - `A.s`: stanza-specific maximum survival (alpha parameter of Beverton-Holt recruitment),
#'  - `B.s`: stanza-specifc carrying capacity (beta parameter of Beverton-Holt recruitment).
#'
#' * Objects created from PVA simulations
#'  - `Nt`: 3-dimensional abundance array (dimensions: time-steps, ages, simulations),
#'  - `Et`: matrix of eggs for each timestep for each year in each simulation,
#'  - `nest`: !!!!!!!!!!!! estimated numbers for each year in each simulation,
#'  - `Vfin`: vector of abundance in the final year across simulations,
#'  - `p.extinct`: vector of proportion of erdicated in each time step,
#'  - `p.extinct.50`: proportion of simulations where the population is eradicated by the 50th time step,
#'  - `p.extinct.100`: proportion of simulations where the population is eradicated by the 100th time step,
#'  - `p.extinct.200`: proportion of simulations where the population is eradicated by the 200th time step,
#'  - `t.extinct`: minimum number of timesteps needed for eradiction,
#'  - `yext.seq`:
#'  - `cost.1`: annual cost of sampling,
#'  - `cost.T`: total cost of sampling (up to the end of simulation or until 100% eradication),
#'  - `NPV`: net present value of sampling (taking into account intergenerational discounting),
#'  - `E.NPV`: expected mean present value (mean of NPV)
#'  - `NT`: abundance in the final time-step (reported as 5th percentile, mean, and 95 percentile of distributions),
#'  - `runtime`: time to execute the PVA,
#'  - `plot`: (Optional) if `create_plot=TRUE`, `plot` is returned as a ggplot object with multiple components. Created by the function `vwReg2.R`.
#'
#' @examples
#' # Run a simple PVA, no custom values or sensitivity testing.
#' PVA(pva.params = inputParameterList)

PVA <- function(params, inits = NULL, custom.inits = NULL,
  sens.pcent = NULL, sens.params = NULL, create.plot = FALSE,
  set.plot.y = NULL, testing=TRUE){
  start <- Sys.time()
  cat("Calculating population projections ...\n")

  if(!is.null(custom.inits)){
    inits <- custom.inits
  } else {
    inits <- init(params)
  }
  AR <- inits$AR
  A <- inits$A
  dt <- inits$dt
  R0 <- inits$R0
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
  R.A <- inits$R.A
  R.B <- inits$R.B
  A.s <- inits$A.s
  B.s <- inits$B.s
  can.a <- inits$can.a
  can.b <- inits$can.b
  M1 <- inits$M1
  Sa.M <- inits$Sa.M
  Ct <- inits$Ct
  Nt <- inits$Nt
  Et <- inits$Et
  nest <- inits$nest

  nT <- inits$nT                   # number of time-steps
  nS <- params$nS                   # number of pre-recruit stanzas
  n.sim <- params$n.sim             # number of simulations
  samp.A <- params$sampA
  # samp.A <- vector()               # time-step within a year that each gear is fished
  r <- params$r                     # discounting rate = future generation discount factor
  G <- params$G                     # generation time (years)
  U.R <- params$U.R # vector()                  # proportion of pre-recruited animals removed per gear
  U.A <- params$U.A # vector()                  # proportion of recruited animals removed per gear
  q.R <- params$q.R # vector()                  # catchability of gear used to remove fish from each pre-recruit stanza
  q.A <- params$q.A # vector()                  # catchability of gear used to remove recruited animals
  n.gear <- params$n.gear           # number of capture gears applied to recruited animals
  t.start.R <- params$t.start.R # vector()            # time-step when sampling begins for each pre-recruit stanza
  t.start.A <- params$t.start.A # vector()            # time-step when sampling begins for each gear used on recruited animals
  E.R <- params$E.R # vector()                  # effort per time-step used to remove fish from each pre-recruit stanza
  E.A <- params$E.A # vector()                  # effort per time-step used to remove recruited animals
  v.a <- params$v.a # vector()                  # Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)
  v.b <- params$v.b # vector()                  # Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
  v.c <- params$v.c # vector()                  # Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)
  v.d <- params$v.d # vector()                  # Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
  C.f.R <- params$C.f.R # vector()
  C.f.A <- params$C.f.A # vector()
  C.E.R <- params$C.E.R # vector()
  C.E.A <- params$C.E.A # vector()
  reck <- params$reck             # recruitment compensation ratio
  p.can <- params$p.can           # proportion of recruit mortality at equilibrium due to cannibalism
  K <- params$K                   # von Bertalanffy metabolic parameter
  afec <- params$afec             # slope of fecundity-weight relationship
  Wmat <- params$Wmat             # weight at maturity
  t.spn <- params$t.spn           # range of time of year when spawning occurs
  V1 <- params$V1                 # initial vulnerable abundance (used to create initial population)
  cann.a <- params$cann.a         # age at which cannibalism on pre-recruits begins
  bet <- params$bet               # rate at which invasives disperse with abundance (between 0 (none) and greater)
  sd.S <- params$sd.S             # standard deviation of environmental effect on survival

  # Apply this when testing sensitivity to biological parameters
  if(!is.null(sens.params)){
    # Skip modification of parameters if these are already in inits
    if(sens.params %in% names(inits) | sens.params %in% c("Bs","Ms")){
      print("Pass")
    } else {
      stopifnot(!is.null(sens.pcent))
      o.value <- get(sens.params)
      new.val <- o.value*sens.pcent
      assign(paste0(sens.params), new.val)
    }
  }

  U.R <- pmin(U.R,0.999999999)
  U.A <- pmin(U.A,0.999999999)
  q.R <- -log(1-U.R)
  q.A <- -log(1-U.A)*V1^bet

  W.st <- rnorm(nT*n.sim,0,sd.S)
  anom <- matrix(data=exp(W.st),nrow=nT,ncol=n.sim)              # survival anomolies for each stanza
  go.R <- matrix(rep(0,nT*nS),nrow=nS)
  for(i in 1:nS)
    go.R[i,t.start.R[i]:nT] <- 1
  go.A <- matrix(rep(0,nT*n.gear),nrow=n.gear)
  for(i in 1:n.gear)
    go.A[i,t.start.A[i]:nT] <- 1

  # dynamics for subsequent years
  for(t in 2:nT){
    N.st <- matrix(nrow=nS+1,ncol=n.sim)    # numbers surviving through each stanza in a year
    Ct.st <- matrix(nrow=nS,ncol=n.sim)
    N.st[1,] <- Et[t-1,]
    V <- colSums(Nt[t-1,(cann.a/dt):A,])
    M0.t <- sweep(can.b,MARGIN=2,V,'*')
    M0.t <- sweep(M0.t,MARGIN=1,can.a,'+')  # density independent variable
    B.st <- M1/M0.t*(1-exp(-M0.t))              # density dependent variable
    for(st in 1:nS){
      Ct.st[st,] <- N.st[st,]*(1-exp(-E.R[st]*q.R[st]*go.R[st,t]))
      N.st[st+1,] <- rbinom(n.sim,round(pmax(0,N.st[st,]),0),
                            pmin(exp(-M0.t[st,]-E.R[st]*q.R[st]*go.R[st,t])*anom[t,]/
                                   (1+B.st[st,]*N.st[st,]),1))
    }
    Ct[t,1:nS] <- rowSums(Ct.st)
    Nt[t,AR,] <- N.st[nS+1,]                # numbers surviving through all stanzas become recruits to the population

    Ct.A <- array(dim=c(n.gear,A-AR+1,n.sim))
    Ft.A <- list()
    Vt <- colSums(Nt[t-1,AR:A,])
    for(i in 1:n.gear){
      q.N <- q.A[i]*pmax(1,Vt)^(-bet)
      if(t*dt-as.integer(t*dt)==samp.A[i]){
        Ft.A[[i]] <- go.A[i,t]*(E.A[i]*sel[i,] %o% q.N)
      } else {
        Ft.A[[i]] <- go.A[i,t]*matrix(rep(0,(A-AR+1)*n.sim),nrow=A-AR+1,ncol=n.sim)
      }
    }
    Ft <- Reduce('+',Ft.A)
    Zt <- Ft - log(Sa.M)
    for(i in 1:n.gear){
      Ct[t,nS+i] <- sum(Nt[t-1,AR:A,] * Ft.A[[i]]/Zt * ( 1 - exp( -Zt)))
    }
    # Get through the loop?
#      if(t<16) cat(t," F\n",Ft[,1],"\nZt",Zt[,1],"\nSa",exp(-Zt[1:(A-AR),1]),"\n")
    Nt[t,(AR+1):A,] <- matrix(rbinom((A-AR)*n.sim,Nt[t-1,AR:(A-1),],exp(-Zt[1:(A-AR),])),ncol=n.sim)
    pairs <- matrix(as.integer(Nt[t,AR:A,]/2),ncol=n.sim)
    Et[t,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
    nest[t,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
  }
  # probabilities of extirpation
  p.extinct.50 <- NA
  p.extinct.100 <- NA
  p.extinct.200 <- NA
  stps <- as.integer(nT/c(4,2,1))
  p.extinct <- sapply(1:nT, function(ii)
    length(which(colSums(Nt[ii,,],na.rm=TRUE)==0))/n.sim) ## calculate p.extinct by time-step
  if( nT >= stps[1] ) p.extinct.50 <- p.extinct[stps[1]]     # probability of extinction in 50 time-steps
  if( nT >= stps[2] ) p.extinct.100 <- p.extinct[stps[2]]  # probability of extinction in 100 time-steps
  if( nT >= stps[3] ) p.extinct.200 <- p.extinct[stps[3]]  # probability of extinction in 200 time-steps
  t.extinct <- min(which(p.extinct==1),nT)
  nyrs <- nT*dt
  Nt[is.na(Nt)]<-0
  y.extinct <- sapply(seq(1/dt,nT,1/dt),function(ii)
    length(which(colSums(Nt[ii,,])==0)))
  y.extinct[2:nyrs] <- y.extinct[2:nyrs]-(y.extinct[1:(nyrs-1)])
  sum.extinct <- sum(y.extinct)
  y.extinct[nyrs] <- n.sim-sum.extinct
  if(is.na(sum.extinct)){
    y.extinct[nyrs] <- n.sim
  }
  if(sum(y.extinct) < n.sim){
    y.extinct[nyrs] <- y.extinct[nyrs]+n.sim-sum(y.extinct)
  }
  yext.seq <- rep(1:nyrs,y.extinct)
  y.extinct <- y.extinct / n.sim

  # numbers remaining
  NT <- quantile(colSums(Nt[nT,,]),probs=c(0.025,0.5,0.975))

  # costs of removal
  x.R <- rep(0,length(q.R))
  x.A <- rep(0,length(q.A))
  x.R[which(E.R>0)] <- 1
  x.A[which(E.A>0)] <- 1
  cost.1 <- sum(E.R*C.E.R) + sum(C.f.R[x.R]) +
    sum(E.A*C.E.A) + sum(C.f.A[x.A])
  y.ext <- t.extinct*dt
  cost.T <- mean(cost.1*(1:nyrs)*y.extinct)
  t.st <- (t.start.R[1]-1)*dt+1
  tseq <- (t.start.R[1]*dt):y.ext
  d <- 1/(1+r)
  df <- d
  NPV <- vector()
  for(i in 1:n.sim){
    NPV[i] <- cost.1*sum(d^(t.st:yext.seq[i])+(df^(t.st:yext.seq[i]))/G)
  }
  E.NPV <- mean(NPV)

  runtime <- Sys.time()-start
  print(runtime)
  out <- list()
  out$phie <- phie
  out$R.A <- R.A
  out$R.B <- R.B
  out$A.s <- A.s
  out$B.s <- B.s
  out$Nt <- Nt
  out$Et <- Et
  out$nest <- nest
  out$Vfin <- (Nt[nT,AR:A,])
  out$p.extinct <- p.extinct
  out$p.extinct.50 <- p.extinct.50
  out$p.extinct.100 <- p.extinct.100
  out$p.extinct.200 <- p.extinct.200
  out$t.extinct <- t.extinct
  out$yext.seq <- yext.seq
  out$cost.1 <- cost.1
  out$cost.T <- cost.T
  out$NPV <- NPV
  out$E.NPV <- E.NPV
  out$NT <- NT
  out$runtime <- runtime
  if(testing){
    out$inits = inits
    return(out)
  } else {
    if(create.plot){
      Na <- apply(Nt,MARGIN=c(1,3),sum,na.rm=TRUE)
      Na <- Na[1:inits$nT,]
      data <- data.frame(Na)
      names(data) <- rep("y",n.sim)
      data$x <- (1:nT)*dt
      cat("Plotting PVA results\n...Probability of extirpation after:\n")
      cat("   ",nT/4*dt," years - ",out$p.extinct.50*100,"%\n")
      cat("   ",nT/2*dt," years - ",out$p.extinct.100*100,"%\n")
      cat("   ",nT*dt," years - ",out$p.extinct.200*100,"%\n")
      # Pulls in vwReg2
      out$plot <- vwReg2(data=data,input=inits,set.ymax=set.ymax)
    }
    return(out)
  }
  gc()
}
