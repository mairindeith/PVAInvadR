#' Rank proposed PVA decision scenarios when biological parameters are increased or decreased by some percentage.
#'
#' @param percent The percentage by which biological parameters () should be modified. Defaults to 15%, such that biological parameters are multiplied by 115% and 85% for upper and lower estimates, respectively.
#' @param decision The output from a decision() function call. Provides the basis for comparing rankings of control scenarios.
#' @param parallel (Optional) Set parallel = TRUE to run analyses using multiple cores, implemented by  Can be useful on multiple-core personal computers.

rank_uncertainty <- function(percent, base_decision, parallel = F){
  scenNames <- base_decision$scenario_name
  upper_var <- 1+percent
  lower_var <- 1-percent

  sens_results <- list()
  cost_T_base <- as.numeric(gsub(",", "", base_decision[,2]))
  print("---      Base cost:      ---")
  print(cost_T_base)
  p_extirp_base <- as.numeric(base_decision[,3])
  print("---      Base p(Extirp):      ---")
  print(p_extirp_base)
  Nt_med_base <- as.numeric(substr(base_decision[,5], start = 1, stop = regexpr(" ", base_decision[,5])-1))
  print("---      Base NT_med:    ---")
  print(Nt_med_base)
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
    cost_T_u <- cost_T_l <- data.frame(
      base = cost_T_base)
    p_extirp_u <- p_extirp_l <- data.frame(
      base = p_extirp_base)
    Nt_med_u <- Nt_med_l <- data.frame(
      base = Nt_med_base)

    for(p in var.par){
      message(paste("...Evaluating changes to ", p, "..."))
      ind <- which(var.par == p)
      # Apply upper transformation +X%
      tmp_inits_upperpper <- PVAInvasR::init(input.params = p, p.cent.trans = upper_var)
      tmp_decision <- as.data.frame(PVAInvasR::decision(custom_inits = tmp_inits_upperpper, direction = "upper", sens_pcent = upper_var, sens_params = p, run_parallel=parallel)) #[[1]])
      print("Upper:")
      print(tmp_decision)
      cost_T_u[,ind+1] <- as.numeric(gsub(",", "", tmp_decision[,2]))
      p_extirp_u[,ind+1] <- as.numeric(tmp_decision[,3])
      Nt_med_u[,ind+1] <- as.numeric(substr(tmp_decision[,5], start = 1, stop = regexpr(" ", tmp_decision[,5])-1))
      colnames(cost_T_u)[ind+1] <- colnames(p_extirp_u)[ind+1] <- colnames(Nt_med_u)[ind+1] <- p
      rm(tmp_decision)
      gc()

      # Apply lower transformation -X%
      tmp_inits_lower <- PVAInvasR::init(input.params = p, p.cent.trans = lower_var)
      tmp_decision_lower <- as.data.frame(PVAInvasR::decision(custom_inits = tmp_inits_lower, direction = "lower", sens_pcent = lower_var, sens_params = p,run_parallel=parallel))
      print("Lower:")
      print(tmp_decision_lower)
      cost_T_l[,ind+1] <- as.numeric(gsub(",", "", tmp_decision_lower[,2]))
      p_extirp_l[,ind+1] <- as.numeric(tmp_decision_lower[,3])
      Nt_med_l[,ind+1] <- as.numeric(substr(tmp_decision_lower[,5], start = 1, stop = regexpr(" ", tmp_decision_lower[,5])-1))
      colnames(cost_T_l)[ind+1] <- colnames(p_extirp_l)[ind+1] <- colnames(Nt_med_l)[ind+1] <- p
      rm(tmp_decision_lower)
      gc()
    }

    rownames(cost_T_u) <- rownames(cost_T_l) <-
      rownames(p_extirp_u) <- rownames(p_extirp_l) <-
      rownames(Nt_med_u) <- rownames(Nt_med_l) <-
      scenNames
    print("---          Post-sensitivity DF         ---")
    print("--- Cost (upper) ---")
    print(cost_T_u)
    print("--- Cost (lower) ---")
    print(cost_T_l)

    # Create ranked and difference ranked data frames
    for(df_name in list("cost_T_u", "cost_T_l", "p_extirp_u", "p_extirp_l", "Nt_med_u", "Nt_med_l")){
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

      assign(paste0(df_name, "_rank"), rank)
      assign(paste0(df_name, "_rankdiff"), rankdiff)
      assign(paste0(df_name, "_diff"), diff)
    }
    save_list = list()
    ind <- 1
    for(cat in c("cost_T", "p_extirp", "NT_med")){
      for(type in c("rank", "rankdiff", "diff")){
        upper <- get(paste0(cat, "_upper_", type))
        lower <- get(paste0(cat, "_lower_", type))
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
        assign(paste0(cat,".",type,"_both"), both_tmp)
        save_list[[ind]] <- paste0(cat, "_upper_", type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat, '_lower_', type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat,".",type,"_both")
        ind <- ind + 1
      }
    }
    # Modify cost so that higher cost = lower ranking
    cost_T_rankdiff_both.mod <- cost_T_rankdiff_both %>%
      mutate(value = value*-1)

    cost_rankdiff_plot <- ggplot_sensitivity(cost_T_rankdiff_both.mod)
    cost_rankdiff_plot <- cost_rankdiff_plot + ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    cost_diff_plot <- ggplot_sensitivity(cost_T_diff_both)
    cost_diff_plot <- cost_diff_plot + ylab("Absolute difference \n compared to the base case ($)") +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    # More animals = lower ranking
    Nt_med_rankdiff_both.mod <- Nt_med_rankdiff_both %>%
      mutate(value = value*-1)
    NT_rankdiff_plot <- ggplot_sensitivity(Nt_med_rankdiff_both)
    NT_rankdiff_plot <- NT_rankdiff_plot + ylab("Rank difference \n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    NT_diff_plot <- ggplot_sensitivity(Nt_med_diff_both)
    NT_diff_plot <- NT_diff_plot + ylab("Absolute difference \n compared to the base case (abundance)") +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    pextirp_rankdiff_plot <- ggplot_sensitivity(p_extirp_rankdiff_both)
    pextirp_rankdiff_plot <- pextirp_rankdiff_plot + ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      labs(color = "Direction of\nvariation") +
      scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    pextirp_diff_plot <- ggplot_sensitivity(p_extirp_diff_both)
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
    data_list$cost_rankchange <- cost_T_rankdiff_both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$cost_absolutechange <- cost_T_diff_both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$nT_rankchange <- Nt_med_rankdiff_both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$nT_absolutechange <- Nt_med_diff_both %>%
      dplyr::rename(AbsoluteChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, AbsoluteChange)

    data_list$pExtirpation_rankchange <- p_extirp_rankdiff_both %>%
      dplyr::rename(RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$pExtirpation_absolutechange <- p_extirp_diff_both %>%
      dplyr::rename(AbsoluteChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction) %>%
      dplyr::select(Scenario, Increase_or_Decrease, Parameter, AbsoluteChange)

    out_list <- list()
    out_list$sensitivity_plots <- plot_list
    out_list$sensitivity_data <- data_list
    return(out_list)
  }
