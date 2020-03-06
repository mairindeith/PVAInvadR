#' Load a .csv template of parameters. A blank template that matches the parameters needed by `PVAInvasR` can be created with the pva_template() function.
#'
#' @param filepath The path to the filled-in .csv template of population and control parameters.
#' @return pva_parameters A list of population and control parameters for use in other `PVAInvasR` functions.
#' @examples
#' load_pva_parameters('~/Documents/PVA_Examples/filledin_parameter_template.csv')

load_pva_parameters <- function(filepath) {
  if (is.null(filepath)){
    # RETURN WITH ERROR
  }
  param.list <- list()
  read.dat <- read_csv(inFile$datapath,
                       header = TRUE,
                       colClasses = "character",
                       sep = ",")

  rename.params <- list('U.R', 'U.A','samp.A','CfR','CfA','ER','EA','CEA','CER','va','vb','vc','vd','tstartA','tstartR') #,'t.start.A','E.R','E.A')
  names(rename.params) <- c('UR', 'UA','sampA','C.f.R','C.f.A','E.R','E.A','C.E.A','C.E.R','v.a','v.b','v.c','v.d','t.start.A','t.start.R') #, 'tstartA','ER','EA')
  params <- read.dat$Parameter
  for(p in 1:length(params)){
    param.name <- params[p]
    if(param.name %in% rename.params){
      rename.index <- which(rename.params == param.name)
      param.name <- names(rename.params)[rename.index]
    }
    param.value <- as.numeric(unlist(strsplit(read.dat[p,2], split=";")))
    param.list[[p]] <- param.value
    names(param.list[p]) <- param.name
    #if(length(param.value)[[1]]==1){
    #  if(is.na(param.value) || param.value == ""){
    #    next
    #  }
    # }
  }
  return(param.list)
}
