#' Load a .csv template of parameters. A blank template that matches the parameters needed by `PVAInvasR` can be created with the pva_template() function.
#' @export
#' @param filepath The path to the filled-in .csv template of population and control parameters.
#' @return pva_parameters A list of population and control parameters for use in other `PVAInvasR` functions.
#' @examples
#' load_pva_parameters('~/Documents/PVA_Examples/filledin_parameter_template.csv')

load_pva_parameters <- function(filepath) {
  if (is.null(filepath)){
    stop("Please specify the filepath for the input PVA parameters.",
      call. = FALSE)
  }
  # For old versions of parameters that use "." delimiters instead of "_"
  dot_params <- c(
      'species',
      'A',
      'AR',
      'nS',
      'nT',
      'dt',
      'n.sim',
      'n.gear',
      't.start.R', #
      't.start.A', #
      # Population parameters
      'V0',
      'reck',
      'p.can',
      'K',
      'afec',
      'Wmat',
      't.spn',
      'Ms',
      'Bs',
      'V1',
      'bet',
      'cann.a',
      'sd.S',
      'U.R', #
      'U.A', #
      'samp.A',
      'E.R', #
      'E.A', #
      'C.f.R', #
      'C.f.A', #
      'C.E.R', #
      'C.E.A', #
      'r',
      'G',
      'v.a', #
      'v.b', #
      'v.c', #
      'v.d', #
      'init.NA'
  )
    # For new versions of pva_template:
  names(dot_params) <- c(
      'species',
      'A',
      'AR',
      'nS',
      'nT',
      'dt',
      'n_sim',
      'n_gear',
      't_start_R', #
      't_start_A', #
      # Population parameters
      'V0',
      'reck',
      'p_can',
      'K',
      'afec',
      'Wmat',
      't_spn',
      'Ms',
      'Bs',
      'V1',
      'bet',
      'cann_a',
      'sd_S',
      'U_R', #
      'U_A', #
      'samp_A',
      'E_R', #
      'E_A', #
      'C_f_R', #
      'C_f_A', #
      'C_E_R', #
      'C_E_A', #
      'r',
      'G',
      'v_a', #
      'v_b', #
      'v_c', #
      'v_d', #
      'init_NA'
    )
  param_list <- list()
  read_dat <- read.csv(filepath, sep=',', stringsAsFactors=FALSE, check.names = F)
  # Replace "'"
  read_dat[] <- lapply(read_dat, function(x) gsub("'", "", x))
  colnames(read_dat) <- gsub("'", "", colnames(read_dat))
  params <- read_dat$Parameters
  for(p in 1:length(params)){
    param_name <- params[p]
    if(param_name %in% dot_params){
      rename_idx <- which(dot_params == param_name)
      param_name <- names(dot_params)[rename_idx]
    }
    if(param_name == "species"){
      param_value <- read_dat[p,2]
    } else {
      param_value <- suppressWarnings(as.numeric(unlist(strsplit(read_dat[p,2], split=";"))))
    }
    param_list[[param_name]] <- param_value
  }
  return(param_list)
}
