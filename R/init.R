#' Initialize parameters for use in a \code{pva()} simulated run; called by the \code{pva()} function.
#'
#'

"init" <- function(input, input.params = NULL, p.cent.trans = NULL){
    #!#
    start <- Sys.time()
    cat("Initializing populations ...\n")
    if(!is.null(p.cent.trans)){
      stopifnot(!is.null(input.params))
      pFact <- p.cent.trans
    }
    if(is.null(p.cent.trans)) {
      pFact <- 1.0
    }
    # dt <- input$dt
    dt <- 1/input$dt               # time-step in years
    R0 <- input$R0                 # unfished equilibrium recruitment
    # TRY:
    V0 <- input$V0
    print(paste0("R0: ", R0))
    print(paste0("V0: ", V0))
    reck <- input$reck             # recruitment compensation ratio
    p.can <- input$p.can           # proportion of recruit mortality at equilibrium due to cannibalism
    A <- as.integer(input$A/dt)    # age at 1% survivorship
    K <- input$K                   # von Bertalanffy metabolic parameter
    afec <- input$afec             # slope of fecundity-weight relationship
    Wmat <- input$Wmat             # weight at maturity
    t.spn <- input$t.spn           # range of time of year when spawning occurs
    Ms <- vector()                 # maximum survival by stanza
    Bs <- vector()                 # stanza-specific density effect (can be interpretted as amount of available habitat)
    V1 <- input$V1                 # initial vulnerable abundance (used to create initial population)
    cann.a <- input$cann.a         # age at which cannibalism on pre-recruits begins
    bet <- input$bet               # rate at which invasives disperse with abundance (between 0 (none) and greater)
    sd.S <- input$sd.S             # standard deviation of environmental effect on survival
    init.Na <- vector()            # initial age-structure in the event those data exist

    nT <- input$nT/dt                # number of R0time-steps
    nS <- input$nS                   # number of pre-recruit stanzas
    AR <- input$AR/dt                # age at recruitment in time-steps
    n.sim <- input$n.sim             # number of simulations
    samp.A <- vector()        # time-step within a year that each gear is fished
    r <- input$r                     # discounting rate = future generation discount factor
    G <- input$G                     # generation time (years)
    U.R <- vector()                  # proportion of pre-recruited animals removed per gear
    U.A <- vector()                  # proportion of recruited animals removed per gear
    q.R <- vector()                  # catchability of gear used to remove fish from each pre-recruit stanza
    q.A <- vector()                  # catchability of gear used to remove recruited animals
    n.gear <- input$n.gear           # number of capture gears applied to recruited animals
    t.start.R <- vector()            # time-step when sampling begins for each pre-recruit stanza
    t.start.A <- vector()            # time-step when sampling begins for each gear used on recruited animals
    E.R <- vector()                  # effort per time-step used to remove fish from each pre-recruit stanza
    E.A <- vector()                  # effort per time-step used to remove recruited animals
    v.a <- vector()                  # Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)
    v.b <- vector()                  # Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
    v.c <- vector()                  # Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)
    v.d <- vector()                  # Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
    C.f.R <- vector()
    C.f.A <- vector()
    C.E.R <- vector()
    C.E.A <- vector()
    for(i in 1:nS){
      Ms[i] <- input$Ms[i]
      Bs[i] <- input$Bs[i]
      U.R[i] <- input$UR[i]
      t.start.R[i] <- input$t.start.R[i]
      E.R[i] <- input$E.R[i]
      C.f.R[i] <- input$C.f.R[i]
      C.E.R[i] <- input$C.E.R[i]
    }
    for(i in 1:n.gear){
      U.A[i] <- input$UA[i]
      t.start.A[i] <- input$t.start.A[i]
      samp.A[i] <- input$sampA[i]
      E.A[i] <- input$E.A[i]
      C.f.A[i] <- input$C.f.A[i]
      C.E.A[i] <- input$C.E.A[i]
      v.a[i] <- input$v.a[i]
      v.b[i] <- input$v.b[i]
      v.c[i] <- input$v.c[i]
      v.d[i] <- input$v.d[i]
    }

    for(i in input$AR:input$A){
      j <- i-input$AR+1
      init.try <- try(init.Na[j] <- input$init.NA[j])
      if (class(init.try) == "try-error") {
        cat("Caught an error reading init NAs, applying default `NA` value.\n")
        init.Na[j] <- "NA"
      }
    }

    # Modify population parameter values provided
    #   in a vector of strings of parameter names
    # This applies during sensitivity testing.
    if(!is.null(input.params)){
      for(params in input.params){
        if(params %in% c("Ms", "Bs")){
          for(i in 1:nS){
            o.value <- get(params)
            assign(params, o.value*pFact)
            Ms[i] <- input$Ms[i]
            Bs[i] <- input$Bs[i]
          }
        }
        o.value <- get(params)
        var.par <- c("reck","p.can","A",
                     "K","afec",
                     "Wmat","Ms","Bs","V1","bet","cann.a","sd.S")
        if(params == "A" || params == "cann.a"){
          newval <- round((o.value*pFact)) #/dt)*dt # round to the nearest time step for discrete measures
        } else {
          newval <- o.value * pFact
        }
        suppressWarnings(assign(params, newval))
      }
    }

    U.R <- pmin(U.R,0.999999999)
    U.A <- pmin(U.A,0.999999999)
    q.R <- -log(1-U.R)
    q.A <- -log(1-U.A)*V1^bet

    # define stanza-specific recruitment parameters
    age <- (AR:A)*dt                                # ages
    la <- 1-exp(-K*age)                             # unfished mean length-at-age
    wa <- la^3                                      # unfished mean weight-at-age
    l0 <- la[1]
    Minf <- log(0.01)*K/(log(l0)-log(l0+exp(K*(A-AR)*dt)-1))   # Lorenzen-based mortality rate for fish at Linf
    spn <- rep(0,A-AR+1)                            # time-steps when spawning occurs
    i <- rep(seq(1,1/dt,length=1/dt),A*dt)[AR:A]
#!#    if(t.spn[1] == t.spn[2]){
#!#      spn[which(i >= t.spn[1]/dt)] <- 1
#!#      spn[which(i >= t.spn[1])] <- 1
#!#    } else {
#!#      spn[which(i >= t.spn[1]/dt & i < t.spn[2]/dt)] <- 1
#!#      spn[which(i >= t.spn[1] & i < t.spn[2])] <- 1
#!#    }
    #!# ALTERNATE FORM THAT MAKES MORE SENSE TO ME
    spn[which(i >= t.spn[1] & i < t.spn[2])] <- 1
    fec <- (pmax(0,wa-Wmat)*afec)                   # unfished eggs at age
    mat <- rep(0,A-AR+1)
    mat[which(fec>0)] <- 1
    init.t <- rep(0,A-AR+1)
    init.t[which(age==as.integer(age))] <- 1        # identifies which season will be evaluated in initial state calculations
    Sa <- (la/(la+exp(K*dt)-1))^(Minf/K)            # length-based survival (from Lorenzen 2000)
    lx <- c(1,Sa[1:(A-AR)])                         # incomplete survivorship to age
    lx <- cumprod(lx)                               # survivorship
    sel <- matrix(nrow=n.gear,ncol=A-AR+1)
    V0 <- c(V0[1], V0[2])
    R0 <- V0/sum(lx*init.t)
    # print(paste0("R0: ", R0))
#!#    R0 <- runif(n.sim, min = 10^R0[1], max = 10^R0[2])        # Unfished recruitment (calculated from V0)
    R0 <- runif(n.sim, min = R0[1], max = R0[2])
    for(i in 1:n.gear){
      sel[i,] <- 1/(1+exp(-(la-v.b[i])/v.a[i]))-1/(1+exp(-(la-v.d[i])/v.c[i]))
      sel[i,] <- sel[i,]/max(sel[i,])
    }
      #sel[i,] <- 1/(1-v.c[i])*((1-v.c[i])/v.c[i])^v.c[i]*exp(v.a[i]*v.c[i]*(v.b[i]-la))/(1+exp(v.a[i]*(v.b[i]-la)))
    phie <- sum(lx*fec*spn/2)                       # unfished eggs per recruit
    R.A <- reck/phie                                # alpha of recruitment function
    R.B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function
    Rinit <- V1/sum(lx*init.t)                      # initial recruitment given initial vulnerable population
    amat <- round(-log(1-Wmat^(1/3))/K,0)           # age-at-maturity (used for differentiating subadults and adults)
    N1 <- V1                                        # initial abundance
    # stanza-specific recruitment parameters
    A.s <- exp(log(R.A)*Ms/sum(Ms))                 # maximum survival for each stanza
    A.s.temp <- c(1,A.s)
    den <- vector()
    for(i in 1:nS)
      den[i] <- Bs[i]*prod(A.s.temp[1:i])
    B.s <- Bs %o% (R.B/sum(den))                      # carrying capacity parameter for each stanza
    # recast BH model as N1=N0*exp(-M0)/(1+M1/M0*(1-exp(-M0))N)
    M0 <- -log(A.s)
    M1 <- sweep(B.s,MARGIN=1,M0/(1-A.s),'*')
    # state M0=a.can+b.can*V where V is sum of animals > Wmat
    V <- rowSums(R0 %o% (lx*init.t)[(cann.a/dt):A-AR+1])        # number of cannibals
    can.a <- M0*(1-p.can)                           # density independent parameter
    can.b <- (M0*p.can) %o% (1/V)                           # density dependent parameter
    sp.t <- rep(init.t[1:(length(age)-1)],100)
    # initialize population
    Nt <- array(rep(0,(nT+AR)*A*n.sim),dim=c((nT+AR),A,n.sim))   # numbers at time and age for each simulation
    init.t <- rep(0,A-AR+1)
    init.t[which(age==as.integer(age))] <- 1
    ifelse(is.na(init.Na[1]),{
      can.a.star <- -log(R.A)*(1-p.can)
      can.b.star <- -log(R.A)*p.can/V
      Ntmp <- matrix(data=c(rep(1e-15,n.sim),rep(0,(A-AR)*n.sim)),nrow=A-AR+1,ncol=n.sim,byrow=TRUE)
      pairs <- Ntmp/2
      Et <- colSums(sweep(pairs,MARGIN=1,fec*spn,'*'))
      M1.t <- R.B*log(R.A)/(R.A-1)
      Nprev <- Ntmp
      t<-2
      while.counter <- 0
      while(mean(colSums(Ntmp))<V1){ # remove , na.rm = T
        while.counter <- while.counter+1
        V <- colSums(Nprev[(cann.a/dt):A-AR+1,])*dt # approximately the number of cannibals per time-step
        M0.t <- can.a.star+can.b.star*V
        Ntmp[1,] <- Et*exp(-M0.t)/(1+M1.t/M0.t*(1-exp(-M0.t))*Et)
        Ntmp[2:(A-AR+1),] <- sweep(Nprev[1:(A-AR),],MARGIN=1,Sa[1:(A-AR)],'*')
        pairs <- Ntmp/2
        Et <- colSums(sweep(pairs,MARGIN=1,fec*sp.t[t],'*'))
        Nprev <- Ntmp
        t<-t+1
      }
      N1 <- V1                                        # initial abundance
      rec.dev <- matrix(exp(rnorm((A-AR+1)*n.sim,0,sd.S)),nrow=A-AR+1,ncol=n.sim)  # annual recruitment deviate
      ifelse(V1<(0.9*median(R0)*sum(lx*init.t)),{
        Nt[1,AR:A,] <- matrix(as.integer(Ntmp*rec.dev*init.t),ncol=n.sim)
      },{
        tmp<-matrix(nrow=A-AR+1,ncol=n.sim)
        for(i in 1:n.sim)
          tmp[,i] <- rmultinom(1,N1,lx*rec.dev[,i]*init.t)       # initial population
        Nt[1,AR:A,] <- tmp
      })
    },{
      N1tmp<-rep(0,A-AR+1)
      N1tmp[which(init.t==1)]<-init.Na
      Ntmp <- matrix(rep(N1tmp,n.sim),nrow=n.sim,
                     ncol=A-AR+1,byrow=TRUE)
      Nt[1,AR:A,] <- t(Ntmp)
    })
    Et <- matrix(data=rep(0,nT*n.sim),nrow=nT,ncol=n.sim)        # eggs for each year in each simulation
    nest <- matrix(data=rep(0,nT*n.sim),nrow=nT,ncol=n.sim)      # nests for each year in each simulation
    Ct <- matrix(nrow=nT,ncol=n.gear+nS)                         # annual catch from each gear
    pairs <- matrix(as.integer(Nt[1,AR:A,]/2),ncol=n.sim)
    ifelse(AR>1,
           {
             Et[AR-1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
             nest[AR-1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
           },  # eggs produced in first year
           {
             Et[1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
             nest[1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
           })
    Sa.M <- matrix(rep(Sa,n.sim),nrow=A-AR+1,ncol=n.sim)
    runtime <- Sys.time()-start
    print(runtime)

    out <- list()
    out$AR <- AR
    out$A <- A
    out$dt <- dt
    out$nT <- nT
    out$R0.vec <- R0
    out$age <- age
    out$la <- la
    out$wa <- wa
    out$spn <- spn
    out$fec <- fec
    out$mat <- mat
    out$Sa <- Sa
    out$Minf <- Minf
    out$lx <- lx
    out$sel <- sel
    out$phie <- phie
    out$R.A <- R.A
    out$R.B <- R.B
    out$A.s <- A.s
    out$B.s <- B.s
    out$can.a <- can.a
    out$can.b <- can.b
    out$M1 <- M1
    out$Sa.M <- Sa.M
    out$Ct <- Ct
    out$Nt <- Nt
    out$Et <- Et
    out$nest <- nest
    out$init.t <- init.t
    print("Done init")
    return(out)
  }
