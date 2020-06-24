#' Rank proposed PVA decision scenarios when biological parameters are increased or decreased by some percentage.
#' @import ggplot2
#'
#' @param df
#' @param ylimit

ggplot_sensitivity <- function(df, ylimit = NaN){
  pars <- c(expression(kappa),
            expression(paste("p"["cann"])),
            expression("A"),
            expression("K"),
            expression(paste("a"["f"])),
            expression(paste("W"["m"])),
            expression(paste({"M"^{"*"}}["s"])),
            #            expression(paste({"M"^{"*"}}["s,2"])),
            expression(paste({"B"^{"*"}}["s"])),
            #            expression(paste({"B"^{"*"}}["s,2"])),
            expression(paste("V"[1])),
            expression(beta),
            expression(paste("cann"["a"])),
            expression(paste(sigma["R"])))
  if(is.na(ylimit)){
    rec_y <- max(abs(df$value))+(0.2*max(abs(df$value)))
    if(rec_y == 0){
      rec_y <- 0.05
    }
  } else {
    rec_y <- ylimit
  }
  tmpplot <- ggplot2::ggplot(df, ggplot2::aes(color = direction)) +
    ggplot2::ylim(-rec_y, rec_y) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0.5, xmax = 1.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 2.5, xmax = 3.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 4.5, xmax = 5.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 6.5, xmax = 7.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 8.5, xmax = 9.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 10.5, xmax = 11.5, ymin = -rec_y, ymax = rec_y), fill = "grey95", color = NaN) + #, alpha = 0.01, color = NaN) +
    ggplot2::geom_hline(yintercept = 0, alpha = 0.5, linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x=new_x, xend=new_x, y=0, yend=value)) +    #, color="orange") +
    #  geom_segment(ggplot2::aes(x=par_num, xend=par_num, y=0, yend=lower.rank.change)) + #, color="turquoise") +
    ggplot2::geom_point(ggplot2::aes(x=new_x, y=value), fill="white", size=4) +
    #  geom_point(ggplot2::aes(x=pars, y=lower.rank.change), color="turquoise", size=4) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      #        panel.grid.major.x = element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size=16),
      axis.title = ggplot2::element_text(size=18),
      strip.background = ggplot2::element_rect(fill = "#337AB7"),
      strip.text = ggplot2::element_text(size=14, colour = "white")
    ) +
    ggplot2::scale_x_discrete(name ="Parameter",
                     breaks = 1:length(pars),
                     limits = 1:length(pars),
                     labels = pars) +
    ggplot2::facet_wrap(~scen, ncol=1)
  return(tmpplot)
}
