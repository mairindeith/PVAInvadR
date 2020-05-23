#' Using imported population and control parameters, run a population viability analysis (PVA) on the target species.
#'
#' @param input List of input parameters, usually created with the load_pva_parameters() function.
#' @param scen_names Character vector of scenario names.
#' @param list (Optional) Boolean; should the output be a named list? See below for list format. Either list or csv must be TRUE.
#' @param csv (Optional) Boolean; should the output be a CSV to be re-loaded into the session? If csv=TRUE, csv_path must be provided.
#' @param csv_path (Optional) If csv=TRUE, where should the decision() template be saved?
#' @param gui (Optional) If TRUE, parameters will be selected by the user via a pop-up menu. Otherwise, selected_params must be provided.
#' @param selected_params (Optional) If gui = FALSE, the user can provide valid named parameters in a character vector here. Alternatively, "all" will create a decision template including all control parameters. Must be provided if gui = FALSE.
#'
#' Outputs either a nested named list (if list = TRUE) or a csv (if csv = TRUE and csv_path provided).
#'  -  List output: a nested named list where the outermost list contains names of scenarios, and within each scenario is a named list of parameters to modify for decision-making.
#'  -  CSV output: a CSV file where the first column indicates scenario names, and subsequent columns indicate control parameters to be varied between scenarios.
#' When setting up input files for the decision, users should leave blank spaces where the parameters should be the same as the input parameters. Only those values that are filled in will be modified and compared with decision().
#'
#' @example
#' # Simple example to create a named list for scenarios "Scen1", "Scen2", and "Scen3".
#'
#' decision_setup(input = load_pva_parameters("parameters.csv", scen_names = c("Scen1", "Scen2"))

decision_setup = function(input, scen_names, list = T, csv = F, csv_path = NULL,
  gui = T, selected_params = NULL){
  # Identify all possible parameters that could be modified by scenarios
  all_params = c(
    paste0("t.start.R.", 1:input$nS),
    paste0("t.start.A.",1:input$n.gear),
    paste0("samp.A.",1:input$n.gear),
    paste0("E.R.",1:input$nS),
    paste0("U.R.",1:input$nS),
    paste0("E.A.",1:input$n.gear)
  )
  # Find original values for each parameter
  all_param_vals = vector(mode = "list", length = length(all_params))
  names(all_param_vals) = all_params

  for(p in names(all_param_vals)){
    param.shortname <- substr(p, start=1, stop=regexpr("\\.[0-9]", p, fixed=F)-1)
    param.num <- as.numeric(substr(p, start=regexpr("\\.[0-9]", p, fixed=F)+1, stop=nchar(p)))
    all_param_vals[[p]] = get(param.shortname)[param.num]
  }

  if(gui == F){
    # If the GUI does not pop up, selected parameters must be provided.
    if(is.null(selected_params)){
      stop("If gui = FALSE, which parameters to vary between scenarios must be provided as a vector of characters, selected_params.")
      return(NULL)
    } else if(selected_params == "all"){
      selected_params = all_params
    } else { # if(!is.null(selected_params))
      # If selected_parameters are not NULL, are they in the set?
      # Check to see if the parameters in selected_params are valid:
      diffs = selected_params[!(which(selected_params %in% all_params))]
      if(!identical(diffs, character(0))){
        stop("Parameter names in selected_params not valid. decision_setup() failed.")
        return(NULL)
      }
    } # if()
  } else { # if(gui == T)
    param_desc=c(
      paste0("Time-step to begin sampling pre-recruits, stanza ",1:input$nS,
        " (t.start.R.",1:input$nS,")"),
      paste0("Time-step to begin sampling adults, gear ",1:input$n.gear,
        " (t.start.A.",1:input$n.gear,")"),
      paste0("Time-step in the year that gear ",1:input$n.gear,
        " is fished (samp.A.",1:input$n.gear,")"),
      paste0("Effort expended by each pre-recruit sampling gear, stanza ",
        1:input$nS," (E.R.",1:input$nS,")"),
      paste0("Proportion of pre-recruits removed with 1 unit of effort, stanza ",
        1:input$nS," (U.R.",1:input$nS,")"),
      paste0("Effort expended by each adult sampling gear, gear ",
        1:input$n.gear," (E.A.",1:input$n.gear,")")
    )
    selected_params = select.list(all_params, multiple = T,
      title = "Select the parameters to vary between scenarios")
    if(length(selected_params) == 0){
      stop("No parameters selected, decision_setup() failed.")
      return(NULL)
    }
  }
  if(list==T){
    scen_list = vector("list", length(scen_names)) # Outermost list, named after scenario names
    names(scen_list) = scen_names
    for(n in scen_names){
      scen_list[[n]] = all_param_vals[names(all_param_vals) %in% selected_params]# Pre-populate with original values
    }
    return(scen_list)
  }
  if(csv==T){
    if(is.null(csv_file)){
      stop("No csv save path provided, decision_setup() failed.")
      return(NULL)
    }
    param_vals = unlist(all_param_vals[names(all_param_vals) %in% selected_params]) # Converts to a numeric here
    scen_df = data.frame(matrix(ncol=length(scen_names), nrow=length(selected_params)), stringsAsFactors=F)
    scen_df[,] = as.numeric(param_vals)
    scen_df = data.frame(cbind(scen_names, t(scen_df)), row.names=F)
    colnames(scen_df) = c("Scenario", selected_params)
    write_csv(scen_df, csv_file)
  }
}
