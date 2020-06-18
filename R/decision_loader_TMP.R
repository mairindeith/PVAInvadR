#' Using imported population and control parameters, run a population viability analysis (PVA) on the target species.
#'
#' @param input List of input parameters, usually created with the load_pva_parameters() function.
#'
#' decision_setup(input = load_pva_parameters("parameters.csv", scen_names = c("Scen1", "Scen2"))

decision_loader = function(csv_input = NULL, list_input = NULL){
  if(is.null(csv_input) && is.null(list_input)){
    stop("One of either csv_input or list_input must be provided.")
    return(NULL)
  }
}
