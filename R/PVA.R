# @title PVA; population viability analysis
#
# @param inits        - Initialized population parameters previously created with the init() function
#                     - If left NULL (default value), init() is run within the function
# @param controls     - Control parameters describing the gear used to remove invasive individuals
# @param parameters   - Population parameters
#
# @return             -

"PVA" <- function(inits = NULL, controls,parameters){
  start <- Sys.time()
  cat("Calculating population projections ...\n")
  # controls is a list including number of time-steps (nT), number of ages (in time-steps) (A),
  #   number of juvenile stanzas (nS), age at recruitment (AR),
  #   and the number of simulations (n.sim)
  # parameters is a list of all parameters
  if(is.null(inits)){ inits = init() } #!# Modify this so that it all makes sense
  R0 <- inits$R0.vec
  age <- inits$age
  la <- inits$la
  wa <- inits$wa
  spn <- inits$spn
  fec <- inits$fec
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

  dt <- controls$dt                   # time-step in years
  nT <- controls$nT/dt                # number of time-steps
  nS <- controls$nS                   # number of pre-recruit stanzas
  AR <- controls$AR/dt                # age at recruitment in time-steps
  n.sim <- controls$n.sim             # number of simulations
  samp.A <- controls$samp.A           # time-step within a year that each gear is fished
  r <- controls$r                     # discounting rate = future generation discount factor
  G <- controls$G                     # generation time (years)
  q.R <- as.vector(controls$q.R)      # catchability of gear used to remove fish from each pre-recruit stanza
  q.A <- as.vector(controls$q.A)      # catchability of gear used to remove recruited fish
  n.gear <- controls$n.gear           # number of capture gears applied to recruited fish
  t.start.R <- controls$t.start.R     # time-step when sampling begins for each pre-recruit stanza
  t.start.A <- controls$t.start.A     # time-step when sampling begins for each gear used on recruited animals
  E.R <- controls$E.R                 # effort per time-step used to remove fish from each pre-recruit stanza
  E.A <- controls$E.A                 # effort per time-step used to remove recruited fish
  v.a <- as.vector(controls$v.a)      # alpha of Thompson 1991 dome-shaped vulnerability (roughly the ascending slope) of removal gear for recruited fish (as a proportion of Linf)
  v.b <- as.vector(controls$v.b)      # beta of Thompson 1991 dome-shaped vulnerability (roughly the length at full selectivity) of removal gear for recruited fish (as a proportion of Linf)
  v.c <- as.vector(controls$v.c)      # gamma of Thompson 1991 dome-shaped vulnerability (roughly the rate of descending decline) of removal gear for recruited fish (0 > v.c > 1)
  C.f.R <- controls$C.f.R
  C.f.A <- controls$C.f.A
  C.E.R <- controls$C.E.R
  C.E.A <- controls$C.E.A

  reck <- parameters$reck             # recruitment compensation ratio
  p.can <- parameters$p.can           # proportion of recruit mortality at equilibrium due to cannibalism
  A <- as.integer(parameters$A/dt)    # age at 1% survivorship
  K <- parameters$K                   # von Bertalanffy metabolic parameter
  afec <- parameters$afec             # slope of fecundity-weight relationship
  Wmat <- parameters$Wmat             # weight at maturity
  t.spn <- parameters$t.spn           # range of time of year when spawning occurs
  Ms <- as.vector(parameters$Ms)      # maximum survival by stanza
  Bs <- as.vector(parameters$Bs)      # stanza-specific density effect (can be interpretted as amount of available habitat)
  V1 <- parameters$V1                 # initial vulnerable abundance (used to create initial population)
  cann.a <- parameters$cann.a         # age at which cannibalism on pre-recruits begins
  bet <- parameters$bet               # rate at which invasives disperse with abundance (between 0 (none) and greater)
  sd.S <- parameters$sd.S             # standard deviation of environmental effect on survival

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
      ifelse(t*dt-as.integer(t*dt)==samp.A[i],
             Ft.A[[i]] <- go.A[i,t]*(E.A[i]*sel[i,] %o% q.N),
             Ft.A[[i]] <- go.A[i,t]*matrix(rep(0,(A-AR+1)*n.sim),nrow=A-AR+1,ncol=n.sim))
    }
    Ft <- Reduce('+',Ft.A)
    Zt <- Ft - log(Sa.M)
    for(i in 1:n.gear){
      Ct[t,nS+i] <- sum(Nt[t-1,AR:A,] * Ft.A[[i]]/Zt * ( 1 - exp( -Zt)))
    }
    Nt[t,(AR+1):A,] <- matrix(rbinom((A-AR)*n.sim,Nt[t-1,AR:(A-1),],exp(-Zt[1:(A-AR),])),ncol=n.sim)
    pairs <- matrix(as.integer(Nt[t,AR:A,]/2),ncol=n.sim)
    Et[t,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
  }
  # probabilities of extirpation
  p.extinct.50 <- NA
  p.extinct.100 <- NA
  p.extinct.200 <- NA
  p.extinct <- sapply(1:nT, function(ii)
    length(which(colSums(Nt[ii,,],na.rm=TRUE)==0))/n.sim) ## calculate p.extinct by year
  if( nT >= 50 ) p.extinct.50 <- p.extinct[50]
  if( nT >= 100 ) p.extinct.100 <- p.extinct[100]
  if( nT >= 200 ) p.extinct.200 <- p.extinct[200]
  t.extinct <- min(which(p.extinct==1),nT)
  nyrs <- nT*dt
  Nt[is.na(Nt)]<-0
  y.extinct <- sapply(seq(1/dt,nT,1/dt),function(ii)
    length(which(colSums(Nt[ii,,])==0)))
  y.extinct[2:nyrs] <- y.extinct[2:nyrs]-(y.extinct[1:(nyrs-1)])
  sum.extinct <- sum(y.extinct)
  y.extinct[nyrs] <- n.sim-sum.extinct
  if(is.na(sum.extinct)) y.extinct[nyrs]<-n.sim
  if(sum(y.extinct)<n.sim)y.extinct[nyrs]<-y.extinct[nyrs]+n.sim-sum(y.extinct)
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
  out$Vfin <- Nt[nT,AR:A,]
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

  return(out)
}
