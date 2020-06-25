#' Initialize parameters from input parameters. This converts parameters into a form to be used by PVA.
#'    Does not need to be called by the user, instead is called by several user-facing functions in the PVAInvadR suite

init <- function(input, input_params = NULL, pcent_trans = NULL, quiet = F){
    start <- Sys.time()
    if(quiet == F){
      message("Initializing populations ...\n")
    }
    if(!is.null(pcent_trans)){
      # stopifnot(!is.null(input_params))
      pFact <- pcent_trans
    }
    if(is.null(pcent_trans)) {
      pFact <- 1.0
    }
    # dt <- input$dt
    if(input$dt > 1){
      dt <- 1/input$dt
    } else {
      dt <- input$dt
    }               # time-step in years
    R0 <- input$R0                 # unfished equilibrium recruitment
    V0 <- input$V0
    reck <- input$reck             # recruitment compensation ratio
    p_can <- input$p_can           # proportion of recruit mortality at equilibrium due to cannibalism
    A <- as.integer(input$A/dt)    # age at 1% survivorship
    K <- input$K                   # von Bertalanffy metabolic parameter
    afec <- input$afec             # slope of fecundity-weight relationship
    Wmat <- input$Wmat             # weight at maturity
    t_spn <- input$t_spn           # range of time of year when spawning occurs
    Ms <- vector()                 # maximum survival by stanza
    Bs <- vector()                 # stanza-specific density effect (can be interpretted as amount of available habitat)
    V1 <- input$V1                 # initial vulnerable abundance (used to create initial population)
    cann_a <- input$cann_a         # age at which cannibalism on pre-recruits begins
    bet <- input$bet               # rate at which invasives disperse with abundance (between 0 (none) and greater)
    sd_S <- input$sd_S             # standard deviation of environmental effect on survival
    init_NA <- vector()            # initial age-structure in the event those data exist

    nT <- input$nT/dt                # number of R0 time-steps
    nS <- input$nS                   # number of pre-recruit stanzas
    AR <- input$AR/dt                # age at recruitment in time-steps
    n_sim <- input$n_sim             # number of simulations
    samp_A <- vector()        # time-step within a year that each gear is fished
    r <- input$r                     # discounting rate = future generation discount factor
    G <- input$G                     # generation time (years)
    U_R <- vector()                  # proportion of pre-recruited animals removed per gear
    U_A <- vector()                  # proportion of recruited animals removed per gear
    q_R <- vector()                  # catchability of gear used to remove fish from each pre-recruit stanza
    q_A <- vector()                  # catchability of gear used to remove recruited animals
    n_gear <- input$n_gear           # number of capture gears applied to recruited animals
    t_start_R <- vector()            # time-step when sampling begins for each pre-recruit stanza
    t_start_A <- vector()            # time-step when sampling begins for each gear used on recruited animals
    E_R <- vector()                  # effort per time-step used to remove fish from each pre-recruit stanza
    E_A <- vector()                  # effort per time-step used to remove recruited animals
    v_a <- vector()                  # Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)
    v_b <- vector()                  # Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
    v_c <- vector()                  # Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)
    v_d <- vector()                  # Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)
    C_f_R <- vector()
    C_f_A <- vector()
    C_E_R <- vector()
    C_E_A <- vector()
    for(i in 1:nS){
      # warning("...i in 1:nS - 61 - ...")
      Ms[i] <- input$Ms[i]
      Bs[i] <- input$Bs[i]
      U_R[i] <- input$U_R[i]
      t_start_R[i] <- input$t_start_R[i]
      E_R[i] <- input$E_R[i]
      C_f_R[i] <- input$C_f_R[i]
      C_E_R[i] <- input$C_E_R[i]
    }
    for(i in 1:n_gear){
      U_A[i] <- input$U_A[i]
      t_start_A[i] <- input$t_start_A[i]
      samp_A[i] <- input$samp_A[i]
      E_A[i] <- input$E_A[i]
      C_f_A[i] <- input$C_f_A[i]
      C_E_A[i] <- input$C_E_A[i]
      v_a[i] <- input$v_a[i]
      v_b[i] <- input$v_b[i]
      v_c[i] <- input$v_c[i]
      v_d[i] <- input$v_d[i]
    }
    warn_counter <- 0
    for(i in input$AR:input$A){
      j <- i-input$AR+1
      init_try <- try(init_NA[j] <- input$init_NA[j])
      if (class(init_try) == "try-error") {
        warn_counter <- warn_counter + 1
        if(warn_counter == 1){
          warning("Caught an error reading init NAs, applying default `NA` value.\n")
        }
        init_NA[j] <- "NA"
      }
    }

    # Modify population parameter values provided
    #   in a vector of strings of parameter names
    # This applies during sensitivity testing.
    if(!is.null(input_params)){
      for(params in input_params){
        if(params %in% c("Ms", "Bs")){
          for(i in 1:nS){
            o_value <- get(params)
            assign(params, o_value*pFact)
            Ms[i] <- input$Ms[i]
            Bs[i] <- input$Bs[i]
          }
        }
        o_value <- get(params)
        vaR_par <- c("reck","p_can","A",
                     "K","afec",
                     "Wmat","Ms","Bs","V1","bet","cann_a","sd_S")
        if(params == "A" || params == "cann_a"){
          newval <- round((o_value*pFact)) #/dt)*dt # round to the nearest time step for discrete measures
        } else {
          newval <- o_value * pFact
        }
        if(params %in% c("p_can", "Wmat")){
          if(newval > 1){
            newval <- 1
          } else if(newval < 0){
            newval <- 0
          }
        }
        assign(params, newval)
      }
    }
    U_R <- pmin(U_R,0.999999999)
    U_A <- pmin(U_A,0.999999999)
    q_R <- -log(1-U_R)
    q_A <- -log(1-U_A)*V1^bet

    # define stanza-specific recruitment parameters
    age <- (AR:A)*dt                                # ages
    la <- 1-exp(-K*age)                             # unfished mean length-at-age
    wa <- la^3                                      # unfished mean weight-at-age
    l0 <- la[1]
    Minf <- log(0.01)*K/(log(l0)-log(l0+exp(K*(A-AR)*dt)-1))   # Lorenzen-based mortality rate for fish at Linf
    spn <- rep(0,A-AR+1)                            # time-steps when spawning occurs
    i <- rep(seq(1,1/dt,length=1/dt),A*dt)[AR:A]
    spn[which(i >= t_spn[1] & i < t_spn[2])] <- 1
    fec <- (pmax(0,wa-Wmat)*afec)                   # unfished eggs at age
    mat <- rep(0,A-AR+1)
    mat[which(fec>0)] <- 1
    init_t <- rep(0,A-AR+1)
    init_t[which(age==as.integer(age))] <- 1        # identifies which season will be evaluated in initial state calculations
    Sa <- (la/(la+exp(K*dt)-1))^(Minf/K)            # length-based survival (from Lorenzen 2000)
    lx <- c(1,Sa[1:(A-AR)])                         # incomplete survivorship to age
    lx <- cumprod(lx)                               # survivorship
    sel <- matrix(nrow=n_gear,ncol=A-AR+1)
    V0 <- c(V0[1], V0[2])
    R0 <- V0/sum(lx*init_t)
    R0 <- runif(n_sim, min = R0[1], max = R0[2])
    for(i in 1:n_gear){
      sel[i,] <- 1/(1+exp(-(la-v_b[i])/v_a[i]))-1/(1+exp(-(la-v_d[i])/v_c[i]))
      sel[i,] <- sel[i,]/max(sel[i,])
    }
    phie <- sum(lx*fec*spn/2)                       # unfished eggs per recruit
    R_A <- reck/phie                                # alpha of recruitment function
    R_B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function
    Rinit <- V1/sum(lx*init_t)                      # initial recruitment given initial vulnerable population
    amat <- round(-log(1-Wmat^(1/3))/K,0)           # age-at-maturity (used for differentiating subadults and adults)
    N1 <- V1                                        # initial abundance
    # stanza-specific recruitment parameters
    A_s <- exp(log(R_A)*Ms/sum(Ms))                 # maximum survival for each stanza
    A_s_temp <- c(1,A_s)
    den <- vector()
    for(i in 1:nS)
      den[i] <- Bs[i]*prod(A_s_temp[1:i])
    B_s <- Bs %o% (R_B/sum(den))                      # carrying capacity parameter for each stanza
    # recast BH model as N1=N0*exp(-M0)/(1+M1/M0*(1-exp(-M0))N)
    M0 <- -log(A_s)
    M1 <- sweep(B_s,MARGIN=1,M0/(1-A_s),'*')
    # state M0=a.can+b.can*V where V is sum of animals > Wmat
    V <- rowSums(R0 %o% (lx*init_t)[(cann_a/dt):A-AR+1])        # number of cannibals
    can_a <- M0*(1-p_can)                           # density independent parameter
    can_b <- (M0*p_can) %o% (1/V)                           # density dependent parameter
    sp_t <- rep(init_t[1:(length(age)-1)],100)
    # initialize population
    Nt <- array(rep(0,(nT+AR)*A*n_sim),dim=c((nT+AR),A,n_sim))   # numbers at time and age for each simulation
    init_t <- rep(0,A-AR+1)
    init_t[which(age==as.integer(age))] <- 1
    ifelse(is.na(init_NA[1]),{
      can_a_star <- -log(R_A)*(1-p_can)
      can_b_star <- -log(R_A)*p_can/V
      Ntmp <- matrix(data=c(rep(1e-15,n_sim),rep(0,(A-AR)*n_sim)),nrow=A-AR+1,ncol=n_sim,byrow=TRUE)
      pairs <- Ntmp/2
      Et <- colSums(sweep(pairs,MARGIN=1,fec*spn,'*'))
      M1.t <- R_B*log(R_A)/(R_A-1)
      Nprev <- Ntmp
      t<-2
      while.counter <- 0
      while(mean(colSums(Ntmp))<V1){ # remove , na.rm = T
        while.counter <- while.counter+1
        V <- colSums(Nprev[(cann_a/dt):A-AR+1,])*dt # approximately the number of cannibals per time-step
        M0.t <- can_a_star+can_b_star*V
        Ntmp[1,] <- Et*exp(-M0.t)/(1+M1.t/M0.t*(1-exp(-M0.t))*Et)
        Ntmp[2:(A-AR+1),] <- sweep(Nprev[1:(A-AR),],MARGIN=1,Sa[1:(A-AR)],'*')
        pairs <- Ntmp/2
        Et <- colSums(sweep(pairs,MARGIN=1,fec*sp_t[t],'*'))
        Nprev <- Ntmp
        t<-t+1
      }
      N1 <- V1                                        # initial abundance
      rec_dev <- matrix(exp(rnorm((A-AR+1)*n_sim,0,sd_S)),nrow=A-AR+1,ncol=n_sim)  # annual recruitment deviate
      ifelse(V1<(0.9*median(R0)*sum(lx*init_t)),{
        Nt[1,AR:A,] <- matrix(as.integer(Ntmp*rec_dev*init_t),ncol=n_sim)
      },{
        tmp<-matrix(nrow=A-AR+1,ncol=n_sim)
        for(i in 1:n_sim)
          tmp[,i] <- rmultinom(1,N1,lx*rec_dev[,i]*init_t)       # initial population
        Nt[1,AR:A,] <- tmp
      })
    },{
      N1tmp<-rep(0,A-AR+1)
      N1tmp[which(init_t==1)]<-init_NA
      Ntmp <- matrix(rep(N1tmp,n_sim),nrow=n_sim,
                     ncol=A-AR+1,byrow=TRUE)
      Nt[1,AR:A,] <- t(Ntmp)
    })
    Et <- matrix(data=rep(0,nT*n_sim),nrow=nT,ncol=n_sim)        # eggs for each year in each simulation
    nest <- matrix(data=rep(0,nT*n_sim),nrow=nT,ncol=n_sim)      # nests for each year in each simulation
    Ct <- matrix(nrow=nT,ncol=n_gear+nS)                         # annual catch from each gear
    pairs <- matrix(as.integer(Nt[1,AR:A,]/2),ncol=n_sim)
    ifelse(AR>1,
           {
             Et[AR-1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
             nest[AR-1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
           },  # eggs produced in first year
           {
             Et[1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,fec*spn,'*')))
             nest[1,] <- as.integer(colSums(sweep(pairs,MARGIN=1,mat*spn,'*')))
           })
    Sa_M <- matrix(rep(Sa,n_sim),nrow=A-AR+1,ncol=n_sim)
    # Prepare outputs
    out <- list()
    out$input_params <- input

    out$initialized_params <- list()

    out$initialized_params$AR <- AR
    out$initialized_params$A <- A
    out$initialized_params$dt <- dt # dt is input
    out$initialized_params$nT <- nT
    out$initialized_params$R0_vec <- R0
    out$initialized_params$age <- age
    out$initialized_params$la <- la
    out$initialized_params$wa <- wa
    out$initialized_params$spn <- spn
    out$initialized_params$fec <- fec
    out$initialized_params$mat <- mat
    out$initialized_params$Sa <- Sa
    out$initialized_params$Minf <- Minf
    out$initialized_params$lx <- lx
    out$initialized_params$sel <- sel
    out$initialized_params$phie <- phie
    out$initialized_params$R_A <- R_A
    out$initialized_params$R_B <- R_B
    out$initialized_params$A_s <- A_s
    out$initialized_params$B_s <- B_s
    out$initialized_params$can_a <- can_a
    out$initialized_params$can_b <- can_b
    out$initialized_params$M1 <- M1
    out$initialized_params$Sa_M <- Sa_M
    out$initialized_params$Ct <- Ct
    out$initialized_params$Nt <- Nt
    out$initialized_params$Et <- Et
    out$initialized_params$nest <- nest
    out$initialized_params$init_t <- init_t
    # out$initialized_params$n_sim <- n_sim
    return(out)
  }
