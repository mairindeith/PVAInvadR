# Helper function used by PVA
#' @import ggplot2
#' @export

vwReg2 <- function(data, input, quiet = FALSE, palette=colorRampPalette(c("purple4","blue","green","yellow","orange","red"), bias=2, space="rgb")(40), set_ymax=NULL){
  dt <- input$dt
  nT <- input$nT*dt
  n_sim <- input$n_sim
  # convert spaghettis to long format
  b2 <- reshape2::melt(as.matrix(data)[,1:n_sim])
  b2$x <- rep(1:nT,n_sim)
  colnames(b2) <- c("index", "rep", "value", "x")
  # Construct ggplot
  # All plot elements are constructed as a list, so they can be added to an existing ggplot
  p0 <- ggplot2::ggplot(data, ggplot2::aes_string(x=data$x, y=data$y)) + ggplot2::theme_bw()
  # p0 <- ggplot2::ggplot(data, x=x, y=y) + theme_bw()
  # initialize elements with NULL (if they are defined, they are overwritten with something meaningful)
  gg.tiles <- NULL
  ymax <- as.integer(log10(max(data)))-1
  ifelse(is.null(set_ymax),
         ylim <- c(0,as.integer(max(data)/ymax+1)*ymax),
         ylim <- c(0,set_ymax))
  if(!quiet){
    cat("Computing density estimates for each vertical cut ...\n")
    flush.console()
    cat("ymax")
    d2 <- plyr::ddply(b2[, c("x", "value")], "x", function(df) { #.(x),
      res <- data.frame(density(df$value, na.rm=TRUE, n=n_sim, bw=ylim[2]/100,from=ylim[1], to=ylim[2])[c("x", "y")])
      colnames(res) <- c("y", "dens")
      return(res)
    }, .progress="text")
  } else {
    d2 <- plyr::ddply(b2[, c("x", "value")], "x", function(df) { #.(x),
      res <- data.frame(density(df$value, na.rm=TRUE, n=n_sim, bw=ylim[2]/100,from=ylim[1], to=ylim[2])[c("x", "y")])
      colnames(res) <- c("y", "dens")
      return(res)
    })
  }
  maxdens <- max(d2$dens,na.rm=TRUE)
  mindens <- min(d2$dens,na.rm=TRUE)
  d2$Density <- (d2$dens - mindens)/maxdens
  d2$x<-data$x[d2$x]
  ## Tile approach
  shade.alpha = 0.1
  d2$alpha.factor <- d2$Density^shade.alpha
  gg.tiles <-  list(ggplot2::geom_tile(data=d2, ggplot2::aes(x=x, y=y, fill=Density,
                                           alpha=alpha.factor)),
                    ggplot2::scale_fill_gradientn("Density\n", colours=palette),
                    ggplot2::scale_alpha_continuous(range=c(0.001, 0.999),guide="none"))
  gg.elements <- list(gg.tiles)
  pOut = p0 + gg.elements + ggplot2::xlab("Year") + ggplot2::ylab("Recruited population numbers") +
           ggplot2::theme(text = ggplot2::element_text(size=18), legend.key.height=ggplot2::unit(2,"cm"))
  return(pOut)
}
