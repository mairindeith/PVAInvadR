#' Rank proposed PVA decision scenarios when biological parameters are increased or decreased by some percentage.
#'
#' @import
#' @param percent The percentage by which biological parameters () should be modified. Defaults to 15%, such that biological parameters are multiplied by 115% and 85% for upper and lower estimates, respectively.
#' @param decision The output from a decision() function call. Provides the basis for comparing rankings of control scenarios.
#' @param parallel (Optional) Set parallel = TRUE to run analyses using multiple cores, implemented by  Can be useful on multiple-core personal computers.

rank_uncertainty <- function(percent, decision, parallel = F){
  scenNames <- decision$scenario.name
  upper.var <- 1+percent
  lower.var <- 1-percent

  sens.results <- list()
  cost.T.base <- as.numeric(gsub(",", "", decision[,2]))
  print("---      Base cost:      ---")
  print(cost.T.base)
  p.extirp.base <- as.numeric(decision[,3])
  print("---      Base p(Extirp):      ---")
  print(p.extirp.base)
  NT.med.base <- as.numeric(substr(decision[,5], start = 1, stop = regexpr(" ", decision[,5])-1))
  print("---      Base NT.med:    ---")
  print(NT.med.base)
  var.par <- c("reck","p.can","A",
               "K","afec",
               "Wmat","Ms","Bs","V1","bet","cann.a","sd.S")
  # Parameter labels
  pars <- c(expression(kappa),
              expression(paste("p"["cann"])),
              expression("A"),
              expression("K"),
              expression(paste("a"["f"])),
              expression(paste("W"["m"])),
              expression(paste({"M"^{"*"}}["s"])),
              #            expression(paste({"M"^{"*"}}["s,1"])),
              #            expression(paste({"M"^{"*"}}["s,2"])),
              expression(paste({"B"^{"*"}}["s"])),
              #            expression(paste({"B"^{"*"}}["s,1"])),
              #            expression(paste({"B"^{"*"}}["s,2"])),
              expression(paste("V"[1])),
              expression(beta),
              expression(paste("cann"["a"])),
              expression(paste(sigma["R"])))
    ### vary parameters
    # one row per scenario
    # Unique data frames for upper (u) and lower (l) uncertainty modification
    cost.T.u <- cost.T.l <- data.frame(
      base = cost.T.base)
    p.extirp.u <- p.extirp.l <- data.frame(
      base = p.extirp.base)
    NT.med.u <- NT.med.l <- data.frame(
      base = NT.med.base)

    for(p in var.par){
      message(paste("...Evaluating changes to ", p, "..."))
      ind <- which(var.par == p)
      # Apply upper transformation +X%
      tmp.inits.upper <- PVAInvasR::init(input.params = p, p.cent.trans = upper.var)
      tmp.decision <- as.data.frame(PVAInvasR::decision(custom.inits = tmp.inits.upper, direction = "upper", sens.pcent = upper.var, sens.params = p, run_parallel=parallel)) #[[1]])
      print("Upper:")
      print(tmp.decision)
      cost.T.u[,ind+1] <- as.numeric(gsub(",", "", tmp.decision[,2]))
      p.extirp.u[,ind+1] <- as.numeric(tmp.decision[,3])
      NT.med.u[,ind+1] <- as.numeric(substr(tmp.decision[,5], start = 1, stop = regexpr(" ", tmp.decision[,5])-1))
      colnames(cost.T.u)[ind+1] <- colnames(p.extirp.u)[ind+1] <- colnames(NT.med.u)[ind+1] <- p
      rm(tmp.decision)
      gc()

      # Apply lower transformation -X%
      tmp.inits.lower <- PVAInvasR::init(input.params = p, p.cent.trans = lower.var)
      tmp.decision.l <- as.data.frame(PVAInvasR::decision(custom.inits = tmp.inits.lower, direction = "lower", sens.pcent = lower.var, sens.params = p,run_parallel=parallel))
      print("Lower:")
      print(tmp.decision.l)
      cost.T.l[,ind+1] <- as.numeric(gsub(",", "", tmp.decision.l[,2]))
      p.extirp.l[,ind+1] <- as.numeric(tmp.decision.l[,3])
      NT.med.l[,ind+1] <- as.numeric(substr(tmp.decision.l[,5], start = 1, stop = regexpr(" ", tmp.decision.l[,5])-1))
      colnames(cost.T.l)[ind+1] <- colnames(p.extirp.l)[ind+1] <- colnames(NT.med.l)[ind+1] <- p
      rm(tmp.decision.l)
      gc()
    }

    rownames(cost.T.u) <- rownames(cost.T.l) <-
      rownames(p.extirp.u) <- rownames(p.extirp.l) <-
      rownames(NT.med.u) <- rownames(NT.med.l) <-
      scenNames
    print("---          Post-sensitivity DF         ---")
    print("--- Cost (upper) ---")
    print(cost.T.u)
    print("--- Cost (lower) ---")
    print(cost.T.l)

    # Create ranked and difference ranked data frames
    for(df_name in list("cost.T.u", "cost.T.l", "p.extirp.u", "p.extirp.l", "NT.med.u", "NT.med.l")){
      df <- get(df_name)
      rank_df <- data.frame(apply(df, 2, order), row.names = rownames(df))
      rank <- tidyr::gather(rank_df)
      print(paste0("rank DF ", df_name))
      print(rank)
      rankdiff <- tidyr::gather(data.frame(rank_df[,1] - rank_df[,2:ncol(rank_df)]))
      print(paste0("rankdiff DF ", df_name))
      print(rankdiff)

      diff <- tidyr::gather(data.frame(df[,1] - df[,2:ncol(df)], row.names = rownames(df)))
      print(paste0("diff DF ", df_name))
      print(diff)

      rank$scen <- paste0("Scenario: ",rownames(df))
      rankdiff$scen <- paste0("Scenario: ",rownames(df))
      diff$scen <- paste0("Scenario: ",rownames(df))
      if(length(grep("l", df_name) == 0)){
        rank$direction <- diff$direction <- rankdiff$direction <- paste0("-",percent)
      } else {
        rank$direction <- diff$direction <- rankdiff$direction <- paste0("+",percent)
      }
      #    save(list="rank", file = "../rankdf.RData")

      assign(paste0(df_name, ".rank"), rank)
      assign(paste0(df_name, ".rankdiff"), rankdiff)
      assign(paste0(df_name, ".diff"), diff)
    }
    save_list = list()
    ind <- 1
    for(cat in c("cost.T", "p.extirp", "NT.med")){
      for(type in c("rank", "rankdiff", "diff")){
        upper <- get(paste0(cat, ".u.", type))
        lower <- get(paste0(cat, ".l.", type))
        both_tmp <- rbind(upper,lower)
        both_tmp$new_x <- NaN
        both_tmp$xbase <- as.numeric(as.factor(both_tmp$key))
        for(r in 1:nrow(both_tmp)){
          if(both_tmp$direction[r] == paste0("+",percent)){
            both_tmp$new_x[r] <- both_tmp$xbase[r] + 0.2
          } else {
            both_tmp$new_x[r] <- both_tmp$xbase[r] - 0.2
          }
        }
        # Create numeric x-axis for each parameter
        assign(paste0(cat,".",type,".both"), both_tmp)
        save_list[[ind]] <- paste0(cat, ".u.", type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat, '.l.', type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat,".",type,".both")
        ind <- ind + 1
      }
    }
    # Modify cost so that higher cost = lower ranking
    cost.T.rankdiff.both.mod <- cost.T.rankdiff.both %>%
      mutate(value = value*-1)

    cost_rankdiff_plot <- ggplot_sensitivity(cost.T.rankdiff.both.mod)
    cost_rankdiff_plot <- cost_rankdiff_plot + ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    cost_diff_plot <- ggplot_sensitivity(cost.T.diff.both)
    cost_diff_plot <- cost_diff_plot + ylab("Absolute difference \n compared to the base case ($)") +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    # More animals = lower ranking
    NT.med.rankdiff.both.mod <- NT.med.rankdiff.both %>%
      mutate(value = value*-1)
    NT_rankdiff_plot <- ggplot_sensitivity(NT.med.rankdiff.both)
    NT_rankdiff_plot <- NT_rankdiff_plot + ylab("Rank difference \n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    NT_diff_plot <- ggplot_sensitivity(NT.med.diff.both)
    NT_diff_plot <- NT_diff_plot + ylab("Absolute difference \n compared to the base case (abundance)") +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    pextirp_rankdiff_plot <- ggplot_sensitivity(p.extirp.rankdiff.both)
    pextirp_rankdiff_plot <- pextirp_rankdiff_plot + ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    pextirp_diff_plot <- ggplot_sensitivity(p.extirp.diff.both)
    pextirp_diff_plot <- pextirp_diff_plot + ylab("Absolute difference\n compared to the base case (P(extirpation))") +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%"))) +
      theme(legend.title = element_text(size=16),
            legend.text = element_text(size=14),
            strip.text = element_text(size=16))

    universal_legend <- get_legend(pextirp_diff_plot + theme(panel.background = element_rect(fill = "white", colour = "white"),
                                                             plot.background = element_rect(fill = "white", colour = "white", size = 3))
    )
    # save("universal_legend", file = "../universal_legend.Rdata")
    plot_list <- list(
      "cost_rankdiff_plot" = cost_rankdiff_plot + theme(legend.position = "none"),
      "cost_absolutediff_plot" = cost_diff_plot + theme(legend.position = "none"),
      "nT_rankdiff_plot" = NT_rankdiff_plot + theme(legend.position = "none"),
      "nT_absolutediff_plot" = NT_diff_plot + theme(legend.position = "none"),
      "pExtirpation_rankdiff_plot" = pextirp_rankdiff_plot + theme(legend.position = "none"),
      "pExtirpation_absolutediff_plot" = pextirp_diff_plot + theme(legend.position = "none"),
      "plot_legend" = universal_legend
    )
    data_list <- list()
    data_list$cost_rankchange <- cost.T.rankdiff.both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, RankChange)

    data_list$cost_absolutechange <- cost.T.diff.both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, RankChange)

    data_list$nT_rankchange <- NT.med.rankdiff.both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, RankChange)

    data_list$nT_absolutechange <- NT.med.diff.both %>%
      dplyr::rename(AbsoluteChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, AbsoluteChange)

    data_list$pExtirpation_rankchange <- p.extirp.rankdiff.both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, RankChange)

    data_list$pExtirpation_absolutechange <- p.extirp.diff.both %>%
      dplyr::rename(AbsoluteChange = value, Scenario = scen, Parameter = key, Increase.or.Decrease = direction) %>%
      dplyr::select(Scenario, Increase.or.Decrease, Parameter, AbsoluteChange)

    out_list <- list()
    out_list$sensitivity_plots <- plot_list
    out_list$sensitivity_data <- data_list
    return(out_list)
  }
