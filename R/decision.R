#' Run multiple PVA simulations initialized with different control scenarios to compare their simulated costs and outcomes.
#'
#' @param custom.inits (Optional) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pva.params}. Should be a named list of po  Can be be outputs of the \code{init()} function from \code{PVAInvas}.
#' @param sens.pcent (Optional) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param create.plot (Optional) Should a ggplot heatmap object also be provided? Default: false. Will return a list with the final entry being the created plot.
#' @param direction (Optional)



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

"decision" <- function(custom.inits = NULL, sens.pcent = NULL, direction = NULL, sens.params = NULL, run_parallel = F){ #, save_pva = F){
  n_cores <- detectCores() - 1
  if(is.na(n_cores)){
    n_cores <- 1
  }
  registerDoParallel(n_cores)
  start <- Sys.time()
  if(!is.null(custom.inits)){
    inits <- custom.inits
  }
  # } else {
  #   inits <- init()
  #   #!# force(inits)
  # }
  dt <- inits$dt
  R0 <- inits$R0.vec
  nT <- inits$nT                   # number of time-steps
  nS <- input$nS                   # number of pre-recruit stanzas
  n.gear <- input$n.gear           # number of capture gears applied to recruited animals
  nR0s <- length(R0)
  scenNames <- values$DF$ScenarioName     # names of scenarios
  tmp.input <- input

  output <- data.frame(row.names=scenNames)
  cost.1 <- list()
  cost.T <- list()
  p.extirp <- list()
  t.extirp <- list()
  t.extirp.u <- list()
  t.extirp.l <- list()
  NT.lcl <- list()
  NT.med <- list()
  NT.ucl <- list()

  # Introduce parallel here
  # How many cores?
#    try(stopCluster(cl), silent = T)
  # New functionality with doParallel instead of parallel
  # Should no longer need OS-specific parallel calls
  if(run_parallel == T){
    n_cores <- detectCores() - 1
    if(is.na(n_cores)){
      n_cores <- 1
    }
    registerDoParallel(n_cores)
    # Don't know if this will work
    vals <- force(values)
    df <- foreach(sn = scenNames, .export = "vals") %dopar% { # .combine = "rbind"
      gc()
      ind <- which(scenNames == sn)
      newParams <- list()
      message("Scenario: ",sn)
      for(c in 2:ncol(vals$DF)){
        col <- colnames(vals$DF)[c]
        param.shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param.num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param.shortname]][[as.numeric(param.num)]] <- as.numeric(vals$DF[ind,c])
      }
      pva <- PVA(input.params = newParams, custom.inits = custom.inits, sens.pcent = sens.pcent, sens.params = sens.params)
      cost.1 <- format(pva$cost.1,big.mark=",",trim=TRUE)
      cost.T <- format(pva$E.NPV,big.mark=",",trim=TRUE)
      p.extirp <- round(pva$p.extinct[nT],2)
      t.extirp <- round(mean(pva$yext.seq),0)
      t.extirp.l <- round(min(pva$yext.seq),0)
      ifelse(max(pva$yext.seq)==nT*dt,
             t.extirp.u <- paste0(round(max(pva$yext.seq),0),"+"),
             t.extirp.u <- paste0(round(max(pva$yext.seq),0)))
      NT.lcl <- round(pva$NT[1],1)
      NT.med <- round(pva$NT[2],1)
      NT.ucl <- round(pva$NT[3],1)
      if(is.null(sens.params)){
        save(pva, file = file.path(as.character(dir_info), "EvaluateScenarios", paste0("Scenario",sn,"_PVAOutputs.RData")))
      } else {
        save(pva, file = file.path(as.character(dir_info), "RankUncertainty", paste0("Scenario",sn,"_",direction,"_",sens.params,"_PVAOutputs.RData")))
      }
      rm(pva)
      gc()
      #d <- cbind(sn,cost.1,cost.T,p.extirp,t.extirp.l,t.extirp,t.extirp.u,NT.lcl,NT.med,NT.ucl)
      #names(d) <- c("Scenario","Annual.cost","Max.total.cost","p.eradicate","min.teradicate",
      #              "t.eradicate","max.teradicate","5%.abundance","50%abundance","95%abundance")
      # d <- cbind(values$DF,d)
      # write.csv(as.data.frame(d),file=paste(input$species,"decision matrix.csv"),quote=FALSE,row.names=FALSE)
      decision_table <- cbind("Scenario Name" = sn,
                              "Annual\ncost ($)" = cost.1,
                              "Total expected\ncost ($)" = cost.T,
                              "Probability\nof eradication" = p.extirp,
                              "Expected time \nto eradication (95% quantiles)" = paste0(t.extirp," (", t.extirp.l,", ", t.extirp.u,")"),
                              "Median abundance\n(95% quantiles)" = paste0(NT.med," (",NT.lcl,", ",NT.ucl,")"))
      decision_table
    #  }
    }
  } else {
    df <- foreach(sn = scenNames, .export = "values") %do% { # .combine = "rbind"
      gc()
      ind <- which(scenNames == sn)
      newParams <- list()
      message("Scenario: ",sn)
      for(c in 2:ncol(values$DF)){
        col <- colnames(values$DF)[c]
        param.shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param.num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param.shortname]][[as.numeric(param.num)]] <- as.numeric(values$DF[ind,c])
      }
      # pva <- PVA(input.params = newParams, custom.inits = custom.inits, sens.pcent = sens.pcent, sens.params = sens.params)
      pva <- PVA(input.params = newParams, custom.inits = custom.inits, sens.pcent = sens.pcent, sens.params = sens.params)
      cost.1 <- format(pva$cost.1,big.mark=",",trim=TRUE)
      cost.T <- format(pva$E.NPV,big.mark=",",trim=TRUE)
      p.extirp <- round(pva$p.extinct[nT],2)
      t.extirp <- round(mean(pva$yext.seq),0)
      t.extirp.l <- round(min(pva$yext.seq),0)
      ifelse(max(pva$yext.seq)==nT*dt,
             t.extirp.u <- paste0(round(max(pva$yext.seq),0),"+"),
             t.extirp.u <- paste0(round(max(pva$yext.seq),0)))
      NT.lcl <- round(pva$NT[1],1)
      NT.med <- round(pva$NT[2],1)
      NT.ucl <- round(pva$NT[3],1)
      if(is.null(sens.params)){
        save(pva, file = file.path(as.character(dir_info), "EvaluateScenarios", paste0("Scenario",sn,"_PVAOutputs.RData")))
      } else {
        save(pva, file = file.path(as.character(dir_info), "RankUncertainty", paste0("Scenario",sn,"_",direction,"_",sens.params,"_PVAOutputs.RData")))
      }
      rm(pva)
      gc()
      decision_table <- cbind("Scenario Name" = sn,
                              "Annual\ncost ($)" = cost.1,
                              "Total expected\ncost ($)" = cost.T,
                              "Probability\nof eradication" = p.extirp,
                              "Expected time \nto eradication (95% quantiles)" = paste0(t.extirp," (", t.extirp.l,", ", t.extirp.u,")"),
                              "Median abundance\n(95% quantiles)" = paste0(NT.med," (",NT.lcl,", ",NT.ucl,")"))
      decision_table
    }
  }
  decision_df <- as.data.frame(df[[1]])
  for(i in 2:length(df)){
    decision_df <- rbind(decision_df, as.data.frame(df[[i]]))
  }
  rm(df)
  gc()
  message("Drawing decision table")
  message(paste0("Runtime: ",Sys.time()-start))
  return(decision_df)
  }
# }
