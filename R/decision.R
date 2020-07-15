#' Run multiple PVA simulations initialized with different control scenarios to compare their simulated costs and outcomes.
#' @import foreach
#' @export
#' @param params Original population and control parameters for the target species and control gear. See `load_pva_parameters`.
#' @param decision_csv (One of decision_csv or decision_list must be provided) Path to a csv file containing scenario-specific control parameters to be used in decision making. May be created with `decision_setup`.
#' @param decision_list (One of decision_csv or decision_list must be provided) An R list containing named parameters and associated values for each scenario. May be created with `decision_setup`.
#' @param custom_inits (Optional, invoked by the `rankUncertainty` function) A vector containing the names of which parameters, if any, should differ from the values provided in \code{pva.params}.
#' @param sens_percent (Optional, invoked by the `rankUncertainty` function) For the sake of sensitivity analysis, how much should population parameters
#' @param sens_percent (Optional, invoked by the `rankUncertainty` function) For the sake of sensitivity analysis, how much should population parameters (i.e. \code{}, \code{}, \code{}, \code{}, \code{}, \code{}, \code{})
#' @param direction (Optional, invoked by the `rankUncertainty` function) Should biological parameters be increased or decreased by  `sens.percent`?
#' @param parallel (Optional), if TRUE decision simulations are run in parallel using  outputs are formatted with dollar signs and commas to "prettify")
#' @param pretty (Optional), if TRUE decision outputs are formatted as in the shiny app, with comma delimiters and dollar signs.

decision <- function(params, decision_csv = NULL, decision_list = NULL, custom_inits = NULL, sens_percent = NULL, direction = NULL, sens_params = NULL, parallel = F, pretty = F, quiet = F){ #, save_pva = F){
  if(is.null(decision_csv) && is.null(decision_list)){
    stop("decision() requires one of decision_csv (path to the filled in decision_csv file) or decision_list (a named list with modified parameters)")
    return(NULL)
  }
  # Convert decision_ info into a readable format
  if(!is.null(decision_csv)){
    decision_setup = readr::read_csv(decision_csv, col_types = readr::cols())
    
  } else if(!is.null(decision_list)){
    scen_names = names(decision_list)
    param_names = names(decision_list[[1]])
    decision_setup = as.data.frame(matrix(nrow=length(scen_names), ncol=length(param_names)+1))
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
  scenNames <- decision_setup$scenario
  if(parallel == T){
    n_cores <- parallel::detectCores() - 1
    if(is.na(n_cores)){
      n_cores <- 1
    }
    doParallel::registerDoParallel(n_cores)
    df <- foreach::foreach(sn = scenNames, .export = "decision_setup", .combine = "rbind") %dopar% {
      sn_idx <- which(scenNames == sn)
      scenParams <- paramss # for each scenario, copy the existing paramss and modify as needed
      if(quiet == F){
        message("Scenario ",sn_idx, ": ", sn)
      }
      for(col_idx in 2:ncol(decision_setup)){
        col <- colnames(decision_setup)[col_idx]
        param_shortname <- substr(col, start=1, stop=regexpr("\\_[0-9]", col, fixed=F)-1)
        param_num <- substr(col, start=regexpr("\\_[0-9]", col, fixed=F)+1, stop=nchar(col))
        scenParams[[param_shortname]][[as.numeric(param_num)]] <- as.numeric(decision_setup[sn_idx,col_idx])
        if(quiet == F){
          message("...param: ", param_shortname, "[", param_num,"] = ", decision_setup[sn_idx,col_idx])
        }
      }
      pva <- PVAInvadR::PVA(params = scenParams, custom_inits = custom_inits, sens_percent = sens_percent, sens_params = sens_params, quiet = quiet)
      if(pretty==T){
        cost_1 <- format(pva$cost_1, big.mark=",", trim=TRUE)
        cost_T <- format(pva$E_NPV, big.mark=",", trim=TRUE)
        p_extirp <- round(pva$p_extinct[params$nT],2)
        t_extirp <- round(mean(pva$yext_seq),0)
        t_extirp_l <- round(min(pva$yext_seq),0)
        ifelse(max(pva$yext_seq)==params$nT*params$dt,
          t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"), # if
          t_extirp_u <- paste0(round(max(pva$yext_seq),0))) # else
        Nt_lcl <- round(pva$NT[1],1)
        Nt_med <- round(pva$NT[2],1)
        Nt_ucl <- round(pva$NT[3],1)
        decision_table <- cbind(
          "Scenario Name" = sn,
          "Annual\ncost ($)" = paste0("$",cost_1),
          "Total expected\ncost ($)" = paste0("$",cost_T),
          "Probability\nof eradication" = p_extirp,
          "Expected time \nto eradication (95% quantiles)" = paste0(t_extirp," (", t_extirp_l,", ", t_extirp_u,")"),
          "Median abundance\n(95% quantiles)" = paste0(Nt_med," (",Nt_lcl,", ",Nt_ucl,")")
        )
      } else {
        cost_1 <- pva$cost_1
        cost_T <- pva$cost_T
        p_extirp <- pva$p_extinct[params$nT]
        t_extirp <- mean(pva$yext_seq)
        t_extirp_l <- min(pva$yext_seq)
        ifelse(max(pva$yext_seq)==params$nT*params$dt,
          t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"), # if
          t_extirp_u <- paste0(round(max(pva$yext_seq),0))) # else
        Nt_lcl <- pva$NT[1]
        Nt_med <- pva$NT[2]
        Nt_ucl <- pva$NT[3]
        decision_table <- data.frame(
          "scenario_name" = sn,
          "annual_cost" = as.numeric(cost_1),
          "total_expected_cost" = as.numeric(cost_T),
          "p_eradication" =  as.numeric(p_extirp),
          "time_to_eradication_mean" = as.numeric(t_extirp),
          "time_to_eradication_2.5" = as.numeric(t_extirp_l),
          "time_to_eradication_97.5" = as.numeric(t_extirp_u),
          "median_abundance" = as.numeric(Nt_med),
          "median_abundance_2.5" = as.numeric(Nt_lcl),
          "median_abundance_97.5" = as.numeric(Nt_ucl)
        )
      }
      decision_table
    }
  } else {
    df <- foreach(sn = scenNames, .combine = "rbind") %do% {
      sn_idx <- which(scenNames == sn)
      scenParams <- paramss
      if(quiet == F){
        message("Scenario: ",sn)
      }
      for(c in 2:ncol(decision_setup)){
        col <- colnames(decision_setup)[c]
        param_shortname <- substr(col, start=1, stop=regexpr("\\_[0-9]", col, fixed=F)-1)
        param_num <- substr(col, start=regexpr("\\_[0-9]", col, fixed=F)+1, stop=nchar(col))
        scenParams[[param_shortname]][[as.numeric(param_num)]] <- as.numeric(decision_setup[sn_idx,c])
      }
      if(pretty==T){
        cost_1 <- format(pva$cost_1, big.mark=",", trim=TRUE)
        cost_T <- format(pva$E_NPV, big.mark=",", trim=TRUE)
        p_extirp <- round(pva$p_extinct[params$nT],2)
        t_extirp <- round(mean(pva$yext_seq),0)
        t_extirp_l <- round(min(pva$yext_seq),0)
        ifelse(max(pva$yext_seq)==params$nT*params$dt,
          t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"), # if
          t_extirp_u <- paste0(round(max(pva$yext_seq),0))) # else
        Nt_lcl <- round(pva$NT[1],1)
        Nt_med <- round(pva$NT[2],1)
        Nt_ucl <- round(pva$NT[3],1)
        decision_table <- cbind(
          "Scenario Name" = sn,
          "Annual\ncost ($)" = paste0("$",cost_1),
          "Total expected\ncost ($)" = paste0("$",cost_T),
          "Probability\nof eradication" = p_extirp,
          "Expected time \nto eradication (95% quantiles)" = paste0(t_extirp," (", t_extirp_l,", ", t_extirp_u,")"),
          "Median abundance\n(95% quantiles)" = paste0(Nt_med," (",Nt_lcl,", ",Nt_ucl,")")
        )
      } else {
        cost_1 <- pva$cost_1
        cost_T <- pva$cost_T
        p_extirp <- pva$p_extinct[params$nT]
        t_extirp <- mean(pva$yext_seq)
        t_extirp_l <- min(pva$yext_seq)
        ifelse(max(pva$yext_seq)==params$nT*params$dt,
          t_extirp_u <- paste0(round(max(pva$yext_seq),0),"+"), # if
          t_extirp_u <- paste0(round(max(pva$yext_seq),0))) # else
        Nt_lcl <- pva$NT[1]
        Nt_med <- pva$NT[2]
        Nt_ucl <- pva$NT[3]
        decision_table <- data.frame(
          "scenario_name" = sn,
          "annual_cost" = as.numeric(cost_1),
          "total_expected_cost" = as.numeric(cost_T),
          "p_eradication" =  as.numeric(p_extirp),
          "time_to_eradication_mean" = as.numeric(t_extirp),
          "time_to_eradication_2.5" = as.numeric(t_extirp_l),
          "time_to_eradication_97.5" = as.numeric(t_extirp_u),
          "median_abundance" = as.numeric(Nt_med),
          "median_abundance_2.5" = as.numeric(Nt_lcl),
          "median_abundance_97.5" = as.numeric(Nt_ucl)
        )
      }
      decision_table
    }
  }
  row.names(df) <- NULL
  return(df)
}
