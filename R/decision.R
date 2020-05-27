#' Run multiple PVA simulations initialized with different control scenarios to compare their simulated costs and outcomes.
#' @import foreach
#'
#' @param custom.inits (Optional, invoked by the `rankUncertainty` function) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pva.params}.
#' @param sens.pcent (Optional, invoked by the `rankUncertainty` function) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param direction (Optional, invoked by the `rankUncertainty` function) Should biological parameters be increased or decreased by  `sens.percent`?

"decision" <- function(input, decision_csv = NULL, decision_list = NULL, sens.pcent = NULL, direction = NULL, sens.params = NULL, run_parallel = F){ #, save_pva = F){
  if(is.null(decision_csv) && is.null(decision_list)){
    stop("decision() requires one of decision_csv (path to the filled in decision_csv file) or decision_list (a named list with modified parameters)")
    return(NULL)
  }
  # Convert decision_ info into a readable format
  if(!is.null(decision_csv)){
    decision_setup = readr::read_csv(decision_csv)
  } else if(!is.null(decision_list)){
    scen_names = names(decision_list)
    param_names = names(decision_list[[1]])
    decision_setup = data.frame(matrix(nrow=length(scen_names), ncol=length(param_names)+1))
    for(s in 1:length(scen_names)){
      decision_setup[s,1] = scen_names[s]
      for(i in 1:length(param_names)){
        param_name = names(decision_list[[s]][i])
        param_val = as.numeric(decision_list[[s]][i])
        decision_setup[s,i+1] = param_val
        colnames(decision_setup) = c("scenario",param_names)
      }
    }
  }

  start <- Sys.time()

  dt <- 1/input$dt
  R0 <- input$R0.vec
  nT <- input$nT                   # number of time-steps
  nS <- input$nS                   # number of pre-recruit stanzas
  n.gear <- input$n.gear           # number of capture gears applied to recruited animals
  nR0s <- length(R0)
  scenNames <- decision_setup$scenario     # names of scenarios
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

  if(run_parallel == T){
    n_cores <- parallel::detectCores() - 1
    if(is.na(n_cores)){
      n_cores <- 1
    }
    doParallel::registerDoParallel(n_cores)
    # Don't know if this will work
    df <- foreach::foreach(sn = scenNames, .export = "decision_setup", .combine = "rbind") %dopar% {
      ind <- which(scenNames == sn)
      newParams <- list()
      message("Scenario: ",sn)
      for(c in 2:ncol(decision_setup)){
        col <- colnames(decision_setup)[c]
        param.shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param.num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param.shortname]][[as.numeric(param.num)]] <- as.numeric(decision_setup[ind,c])
      }
      cat(newParams)
      pva <- PVA(params = newParams, custom.inits = custom.inits, sens.pcent = sens.pcent, sens.params = sens.params)
      cat(data.frame(pva$cost.l))
      cat("\nPVA dopar DONE\n")
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
      # rm(pva)
      decision_table <- cbind("Scenario Name" = sn,
                              "Annual\ncost ($)" = cost.1,
                              "Total expected\ncost ($)" = cost.T,
                              "Probability\nof eradication" = p.extirp,
                              "Expected time \nto eradication (95% quantiles)" = paste0(t.extirp," (", t.extirp.l,", ", t.extirp.u,")"),
                              "Median abundance\n(95% quantiles)" = paste0(NT.med," (",NT.lcl,", ",NT.ucl,")"))
      data.frame(decision_table)
    }
  } else {
    df <- foreach(sn = scenNames, .combine = "rbind") %do% { # .combine = "rbind"
      gc()
      ind <- which(scenNames == sn)
      newParams <- list()
      message("Scenario: ",sn)
      for(c in 2:ncol(decision_setup)){
        col <- colnames(decision_setup)[c]
        param.shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param.num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param.shortname]][[as.numeric(param.num)]] <- as.numeric(decision_setup[ind,c])
      }
      cat(data.frame(newParams))
      pva <- PVA(params = newParams, custom.inits = custom.inits, sens.pcent = sens.pcent, sens.params = sens.params)
      cat(data.frame(pva$cost.l))
      cat("\nPVA dopar DONE\n")
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
      rm(pva)
      decision_table <- cbind("Scenario Name" = sn,
                              "Annual\ncost ($)" = cost.1,
                              "Total expected\ncost ($)" = cost.T,
                              "Probability\nof eradication" = p.extirp,
                              "Expected time \nto eradication (95% quantiles)" = paste0(t.extirp," (", t.extirp.l,", ", t.extirp.u,")"),
                              "Median abundance\n(95% quantiles)" = paste0(NT.med," (",NT.lcl,", ",NT.ucl,")"))
      data.frame(decision_table)
    }
  }
  return(df)
  }
# }
