#' Run multiple PVA simulations initialized with different control scenarios to compare their simulated costs and outcomes.
#' @import foreach
#'
#' @param input
#' @param decision_csv
#' @param decision_list
#' @param custom_inits (Optional, invoked by the `rankUncertainty` function) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pva.params}.
#' @param sens_percent (Optional, invoked by the `rankUncertainty` function) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param direction (Optional, invoked by the `rankUncertainty` function) Should biological parameters be increased or decreased by  `sens.percent`?

decision <- function(input, decision_csv = NULL, decision_list = NULL, custom_inits = NULL, sens_percent = NULL, direction = NULL, sens_params = NULL, run_parallel = F){ #, save_pva = F){
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
  n_gear<- input$n_gear          # number of capture gears applied to recruited animals
  nR0s <- length(R0)
  scenNames <- decision_setup$scenario     # names of scenarios
  tmp_input <- input

  output <- data.frame(row.names=scenNames)
  cost_1 <- list()
  cost_T <- list()
  p_extirp <- list()
  t_extirp <- list()
  t_extirp_u <- list()
  t_extirp_l <- list()
  Nt_lcl <- list()
  Nt_med <- list()
  Nt_ucl <- list()

  if(run_parallel == T){
    n_cores <- parallel::detectCores() - 1
    if(is.na(n_cores)){
      n_cores <- 1
    }
    doParallel::registerDoParallel(n_cores)
    df <- foreach::foreach(sn = scenNames, .export = "decision_setup", .combine = "rbind") %dopar% {
      ind <- which(scenNames == sn)
      newParams <- list()
      message("Scenario: ",sn)
      for(c in 2:ncol(decision_setup)){
        col <- colnames(decision_setup)[c]
        param_shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param_num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param_shortname]][[as.numeric(param_num)]] <- as.numeric(decision_setup[ind,c])
      }
      pva <- PVA(params = newParams, custom_inits = custom_inits, sens_percent = sens_percent, sens_params = sens_params)
      cost_1 <- format(pva$cost_1, big.mark=",", trim=TRUE)
      cost_T <- format(pva$E_NPV, big.mark=",", trim=TRUE)
      p_extirp <- round(pva$p_extinct[nT],2)
      t_extirp <- round(mean(pva$yext_seq),0)
      t_extirp_l <- round(min(pva$yext_seq),0)
      ifelse(max(pva$yext_seq)==nT*dt,
             t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"),
             t_extirp_u <- paste0(round(max(pva$yext_seq),0)))
      Nt_lcl <- round(pva$NT[1],1)
      Nt_med <- round(pva$NT[2],1)
      Nt_ucl <- round(pva$NT[3],1)
      # rm(pva)
      decision_table <- cbind("scenario_name" = sn,
                              "annual_cost" = cost_1,
                              "total_expected_cost" = cost_T,
                              "p_eradication" = p_extirp,
                              "time_to_eradication_95" = paste0(t_extirp," (", t_extirp_l,", ", t_extirp_u,")"),
                              "median_abundance_95" = paste0(Nt_med," (",Nt_lcl,", ",Nt_ucl,")"))
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
        param_shortname <- substr(col, start=1, stop=regexpr("\\.[0-9]", col, fixed=F)-1)
        param_num <- substr(col, start=regexpr("\\.[0-9]", col, fixed=F)+1, stop=nchar(col))
        newParams[[param_shortname]][[as.numeric(param_num)]] <- as.numeric(decision_setup[ind,c])
      }
      # cat(data.frame(newParams))
      pva <- PVA(params = newParams, custom_inits = custom_inits, sens_percent = sens_percent, sens_params = sens_params)
      # cat(data.frame(pva$cost_l))
      # cat("\nPVA dopar DONE\n")
      cost_1 <- format(pva$cost_1,big.mark=",",trim=TRUE)
      cost_T <- format(pva$E_NPV,big.mark=",",trim=TRUE)
      p_extirp <- round(pva$p_extinct[nT],2)
      t_extirp <- round(mean(pva$yext_seq),0)
      t_extirp_l <- round(min(pva$yext_seq),0)
      ifelse(max(pva$yext_seq)==nT*dt,
             t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"),
             t_extirp_u <- paste0(round(max(pva$yext_seq),0)))
      Nt_lcl <- round(pva$NT[1],1)
      Nt_med <- round(pva$NT[2],1)
      Nt_ucl <- round(pva$NT[3],1)
      rm(pva)
      decision_table <- cbind("scenario_name" = sn,
                              "annual_cost" = cost_1,
                              "total_expected_cost" = cost_T,
                              "p_eradication" = p_extirp,
                              "time_to_eradication_95" = paste0(t_extirp," (", t_extirp_l,", ", t_extirp_u,")"),
                              "median_abundance_95" = paste0(Nt_med," (",Nt_lcl,", ",Nt_ucl,")"))
      data.frame(decision_table)
    }
  }
  return(df)
  }
