#' Rank proposed PVA decision scenarios when biological parameters are increased or decreased by some percentage.
#' @import ggplot2
#'
#' @param params The output from a decision() function call. Provides the basis for comparing rankings of control scenarios.
#' @param decision The output from a decision() function call. Provides the basis for comparing rankings of control scenarios.
#' @param percent The percentage by which biological parameters should be modified. Defaults to 15%, such that biological parameters are multiplied by 115% and 85% for upper and lower estimates, respectively.
#' @param parallel (Optional) Set parallel = TRUE to run analyses using multiple cores.

rank_uncertainty <- function(params, percent=0.15, decision_csv = NULL, decision_list = NULL, parallel = FALSE, quiet = FALSE){
  if(!quiet){
    message("Running base decision scenario (biological parameters at original values)")
  }
  base_decision <- PVAInvadR::decision(params, decision_csv = decision_csv, decision_list = decision_list, parallel = parallel, pretty = F, quiet = quiet)
  scenNames <- base_decision$scenario_name
  upper_var <- 1+percent
  lower_var <- 1-percent
  sens_results <- list()
  if(quiet == F){
    message(paste0(rep("-", 40)))
    message("Base costs: \n   ", paste0("Scenario ", scenNames, ": ", base_decision$annual_cost, "\n   "))
    message("Base probabilities of extirpation: \n   ", paste0("Scenario ", scenNames, ": ", base_decision$p_eradication, "\n   "))
    message("Base median abundances (Nt): \n   ", paste0("Scenario ", scenNames, ": ", base_decision$median_abundance, "\n   "))
  }
  cost_T_base <- base_decision$annual_cost
  p_extirp_base <- base_decision$p_eradication
  Nt_med_base <- base_decision$median_abundance

  # Jun 17 Debugging:
  ### Indicates functional parameters that do not cause errors
  # Indicates "no loop for break/next, jumping to top level"
  var_par <- c("reck",
               "p_can",
               "A",
               "K",
               "afec",
               "Wmat",
               "Ms",
               "Bs",
               "V1",
               "bet",
               # "cann_a",
               "sd_S"
               )
    ### vary parameters
    # one row per scenario
    # Unique data frames for upper (u) and lower (l) uncertainty modification
    cost_T_u <- cost_T_l <- data.frame(
      base = cost_T_base)
    p_extirp_u <- p_extirp_l <- data.frame(
      base = p_extirp_base)
    Nt_med_u <- Nt_med_l <- data.frame(
      base = Nt_med_base)
    for(p in var_par){
      if(!quiet){
        message("\n", paste0(rep("-", 40)))
        message(paste("--- EVALUATING CHANGES TO ", p, " --- "))
        message(paste0(rep("-", 40)), "\n")
      }
      ind <- which(var_par == p)
      # Apply upper transformation +X%
      tmp_inits_upper <- PVAInvadR::init(params, params_params = p, pcent_trans = upper_var, quiet = quiet)
      tmp_decision_upper <- PVAInvadR::decision(params, decision_csv = decision_csv, decision_list = decision_list, custom_inits = tmp_inits_upper, direction = "upper", sens_percent = upper_var, sens_params = p, parallel = parallel, quiet = quiet)

      cost_T_u[,ind+1] <- tmp_decision_upper$annual_cost
      p_extirp_u[,ind+1] <- tmp_decision_upper$p_eradication
      Nt_med_u[,ind+1] <- tmp_decision_upper$median_abundance
      colnames(cost_T_u)[ind+1] <- colnames(p_extirp_u)[ind+1] <- colnames(Nt_med_u)[ind+1] <- p
      if(quiet == F){
        message("--- Increase ", p, " by ", percent*100, "% --- ")
        message("New costs: \n      ", paste0("Scenario ", scenNames, ": ", tmp_decision_upper$annual_cost, "\n      "))
        message("New p(Extirpation)s: \n      ", paste0("Scenario ", scenNames, ": ", tmp_decision_upper$p_eradication, "\n      "))
        message("New median abundance (Nt): \n      ", paste0("Scenario ", scenNames, ": ", tmp_decision_upper$median_abundance, "\n      "))
      }
      rm(tmp_decision_upper)
      gc()
      # Apply lower transformation -X%
      tmp_inits_lower <- PVAInvadR::init(params, params_params = p, pcent_trans = lower_var, quiet = quiet)
      tmp_decision_lower <- PVAInvadR::decision(params, decision_csv = decision_csv, decision_list = decision_list, custom_inits = tmp_inits_lower, direction = "lower", sens_percent = lower_var, sens_params = p, parallel = parallel, quiet = quiet)
      cost_T_l[,ind+1] <- tmp_decision_lower$annual_cost
      p_extirp_l[,ind+1] <- tmp_decision_lower$p_eradication
      Nt_med_l[,ind+1] <- tmp_decision_lower$median_abundance
      colnames(cost_T_l)[ind+1] <- colnames(p_extirp_l)[ind+1] <- colnames(Nt_med_l)[ind+1] <- p
      if(quiet == F){
        message("--- Decrease ", p, " by ", percent*100, "% --- ")
        message("New costs:", paste0(" \n      Scenario ", scenNames, ": ", tmp_decision_lower$annual_cost))
        message("New p(Extirpation):", paste0(" \n      Scenario ", scenNames, ": ", tmp_decision_lower$p_eradication))
        message("New median abundance (Nt):", paste0("\n      Scenario ", scenNames, ": ", tmp_decision_lower$median_abundance))
      }
      rm(tmp_decision_lower)
      gc()
    }
    if(quiet == F){
      message(paste0(rep("-", 40)))
    }
    rownames(cost_T_u) <- rownames(cost_T_l) <-
      rownames(p_extirp_u) <- rownames(p_extirp_l) <-
      rownames(Nt_med_u) <- rownames(Nt_med_l) <-
      scenNames

    # Create ranked and difference ranked data frames
    for(df_name in list("cost_T_u", "cost_T_l", "p_extirp_u", "p_extirp_l", "Nt_med_u", "Nt_med_l")){
      df <- get(df_name)
      rank_df <- data.frame(apply(df, 2, order), row.names = rownames(df))
      rank <- data.frame(tidyr::gather(rank_df))
      rankdiff <- tidyr::gather(rank_df[,1] - rank_df[,2:ncol(rank_df)])
      diff <- tidyr::gather(df[,1] - df[,2:ncol(df)])
      rank$scen <- paste0("Scenario: ",rownames(df))
      rankdiff$scen <- paste0("Scenario: ",rownames(df))
      diff$scen <- paste0("Scenario: ",rownames(df))
      if(length(grep("l", df_name) == 0)){
        rank$direction <- diff$direction <- rankdiff$direction <- paste0("-",percent)
      } else {
        rank$direction <- diff$direction <- rankdiff$direction <- paste0("+",percent)
      }
      assign(paste0(df_name, "_rank"), rank)
      assign(paste0(df_name, "_rankdiff"), rankdiff)
      assign(paste0(df_name, "_diff"), diff)
    }
    save_list = list()
    ind <- 1
    for(cat in c("cost_T", "p_extirp", "Nt_med")){
      for(type in c("rank", "rankdiff", "diff")){
        upper <- get(paste0(cat, "_u_", type))
        lower <- get(paste0(cat, "_l_", type))
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
        assign(paste0(cat,"_",type,"_both"), both_tmp)
        save_list[[ind]] <- paste0(cat, "_u_", type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat, '_l_', type)
        ind <- ind + 1
        save_list[[ind]] <- paste0(cat,"_",type,"_both")
        ind <- ind + 1
      }
    }
    # Modify cost so that higher cost = lower ranking
    cost_T_rankdiff_both.mod <- dplyr::mutate(cost_T_rankdiff_both, value = value*-1)

    cost_rankdiff_plot <- ggplot_sensitivity(cost_T_rankdiff_both.mod)
    cost_rankdiff_plot <- cost_rankdiff_plot + ggplot2::ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    cost_diff_plot <- ggplot_sensitivity(cost_T_diff_both)
    cost_diff_plot <- cost_diff_plot + ggplot2::ylab("Absolute difference \n compared to the base case ($)") +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    # More animals = lower ranking
    Nt_med_rankdiff_both.mod <- dplyr::mutate(Nt_med_rankdiff_both, value = value*-1)
    NT_rankdiff_plot <- ggplot_sensitivity(Nt_med_rankdiff_both)
    NT_rankdiff_plot <- NT_rankdiff_plot + ggplot2::ylab("Rank difference \n compared to the base case") +
      #    ylim(c(-2,2)) +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    NT_diff_plot <- ggplot_sensitivity(Nt_med_diff_both)
    NT_diff_plot <- NT_diff_plot + ggplot2::ylab("Absolute difference \n compared to the base case (abundance)") +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))

    pextirp_rankdiff_plot <- ggplot_sensitivity(p_extirp_rankdiff_both)
    pextirp_rankdiff_plot <- pextirp_rankdiff_plot + ggplot2::ylab("Rank difference\n compared to the base case") +
      #    ylim(c(-2,2)) +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%")))
    pextirp_diff_plot <- ggplot_sensitivity(p_extirp_diff_both)
    pextirp_diff_plot <- pextirp_diff_plot + ggplot2::ylab("Absolute difference\n compared to the base case (P(extirpation))") +
      ggplot2::labs(color = "Direction of\nvariation") +
      ggplot2::scale_color_hue(labels = c(paste0("Decrease by ", percent*100, "%"),
                                 paste0("Increase by ", percent*100, "%"))) +
      ggplot2::theme(legend.title = ggplot2::element_text(size=16),
            legend.text = ggplot2::element_text(size=14),
            strip.text = ggplot2::element_text(size=16))

    universal_legend <- cowplot::get_legend(pextirp_diff_plot + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
                                                             plot.background = ggplot2::element_rect(fill = "white", colour = "white", size = 3))
    )
    # save("universal_legend", file = "../universal_legend.Rdata")
    plot_list <- list(
      "cost_rankdiff_plot" = cost_rankdiff_plot + ggplot2::theme(legend.position = "none"),
      "cost_absolutediff_plot" = cost_diff_plot + ggplot2::theme(legend.position = "none"),
      "nT_rankdiff_plot" = NT_rankdiff_plot + ggplot2::theme(legend.position = "none"),
      "nT_absolutediff_plot" = NT_diff_plot + ggplot2::theme(legend.position = "none"),
      "pExtirpation_rankdiff_plot" = pextirp_rankdiff_plot + ggplot2::theme(legend.position = "none"),
      "pExtirpation_absolutediff_plot" = pextirp_diff_plot + ggplot2::theme(legend.position = "none"),
      "plot_legend" = universal_legend
    )
    data_list <- list()
    data_list$cost_rankchange <- dplyr::rename(cost_T_rankdiff_both, RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$cost_rankchange <- dplyr::select(data_list$cost_rankchange, Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$cost_absolutechange <- dplyr::rename(cost_T_diff_both, RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$cost_absolutechange <-  dplyr::select(data_list$cost_absolutechange, Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$nT_rankchange <- dplyr::rename(Nt_med_rankdiff_both, RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$nT_rankchange <- dplyr::select(data_list$nT_rankchange, Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$nT_absolutechange <- dplyr::rename(Nt_med_diff_both, AbsoluteChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$nT_absolutechange <- dplyr::select(data_list$nT_absolutechange, Scenario, Increase_or_Decrease, Parameter, AbsoluteChange)

    data_list$pExtirpation_rankchange <- dplyr::rename(p_extirp_rankdiff_both, RankChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$pExtirpation_rankchange <- dplyr::select(data_list$pExtirpation_rankchange, Scenario, Increase_or_Decrease, Parameter, RankChange)

    data_list$pExtirpation_absolutechange <- dplyr::rename(p_extirp_diff_both, AbsoluteChange = value, Scenario = scen, Parameter = key, Increase_or_Decrease = direction)
    data_list$pExtirpation_absolutechange <- dplyr::select(data_list$pExtirpation_absolutechange, Scenario, Increase_or_Decrease, Parameter, AbsoluteChange)

    out_list <- list()
    out_list$sensitivity_plots <- plot_list
    out_list$sensitivity_data <- data_list
    return(out_list)
  }
