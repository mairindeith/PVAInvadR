#' Plot the outputs of a PVA simulation.
#' @export
#' @param pva List of outputs created by running the PVA() function.
#' @return List of
#' * A list of PVA outputs, including calculated parameters:
#'  - `phie`: unfished eggs per recruit,
#'  - `R.A` and `R.B`: alpha and beta paramters for a Beverton-Holt recruitment function,
#'  - `Rinit`: initial recruitment,
#'  - `amat`: age at maturity to differentiate adults and sub-adults,
#'  - `A.s`: maximum survival in each stanza, and
#'  - `B.s`: carrying capacity parameter for each stanza)
#'
#' and PVA-generated calculations for each iteration of the simulation
#'  - `N1`: initial population size per simulation,
#'
#' )

plot_pva <- function(pva, set.ymax=NULL){
# "heat.proj" <- function(pva){
  n.sim <- pva$n.sim
  dt <- pva$dt                   # time-step in years
  nT <- pva$nT
  Na <- apply(pva$Nt,MARGIN=c(1,3),sum,na.rm=TRUE)
  Na <- Na[1:nT,]
  data <- data.frame(Na)
  names(data) <- rep("y",n.sim)
  data$x <- (1:nT)*dt
  cat("Plotting PVA results\n...Probability of extirpation after:\n")
  cat("   ",nT/4*dt," years - ",pva$p.extinct.50*100,"%\n")
  cat("   ",nT/2*dt," years - ",pva$p.extinct.100*100,"%\n")
  cat("   ",nT*dt," years - ",pva$p.extinct.200*100,"%\n")
  out <- list()
  out$dat <- pva
  if(!is.null(set.ymax)){
    set.ymax<-input$set.ymax
  }
  out$plot <- vwReg2(data=data,input=pva,set.ymax=set.ymax)
  return(out)
}
