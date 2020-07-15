# Lazy temporary script for testing PVA functions individually.
# To be run in the root directory of the package.

library(devtools)
setwd("~/pCloudDrive/Shiny_vanPoorten/PVAInvasR")
load_all() # Load the PVAInavsR package
inputs = load_pva_parameters("../ShinyDropbox/TestInputs/FuseeLMB_population_parameters.csv")

# Jump right in with a test?
pvatest = PVA(inputs, create_plot = T)

# names(pvatest)
# [1] "initialized_params" "phie"               "R.A"
# [4] "R.B"                "A.s"                "B.s"
# [7] "Nt"                 "Et"                 "nest"
# [10] "Vfin"               "p.extinct"          "p.extinct.50"
# [13] "p.extinct.100"      "p.extinct.200"      "t.extinct"
# [16] "yext.seq"           "cost.1"             "cost.T"
# [19] "NPV"                "E.NPV"              "NT"
# [22] "runtime"

# scens = c("A","B")
csvpath = "decision_scen_template2.csv"

### Instead of using the GUI, provide a list:
# selected_params = c("t_start_R_1", "E_A_1")
### Can also provide a single parameter
# selected_params = "E.A.1"

### decision_setup(inputs, scens, csv = T, list = F, csv_path = csvpath, gui = F, selected_params = selected_params)

csvpath_filled = "decision_scen_template2_filled.csv"

dectest = decision(inputs, decision_csv = csvpath_filled, parallel = T, pretty = F)

ranktest = rank_uncertainty(inputs, percent=0.15, decision_csv = csvpath_filled, decision_list = NULL, parallel = T)

ranktest_quiet = rank_uncertainty(inputs, percent=0.15, decision_csv = csvpath_filled, decision_list = NULL, parallel = T, quiet = T)

# Next time, check out what's happening with the ranktest outputs:
### $sensitivity_data
### $sensitivity_data$cost_rankchange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                   Parameter RankChange
### 1 rank_df...1....rank_df...2.ncol.rank_df..          0
### 2 rank_df...1....rank_df...2.ncol.rank_df..          0
### 3 rank_df...1....rank_df...2.ncol.rank_df..          0
### 4 rank_df...1....rank_df...2.ncol.rank_df..          0
###
### $sensitivity_data$cost_absolutechange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                                           Parameter RankChange
### 1 as.numeric.as.character.df...1......as.numeric.as.character.df...          0
### 2 as.numeric.as.character.df...1......as.numeric.as.character.df...          0
### 3 as.numeric.as.character.df...1......as.numeric.as.character.df...          0
### 4 as.numeric.as.character.df...1......as.numeric.as.character.df...          0
###
### $sensitivity_data$nT_rankchange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                   Parameter RankChange
### 1 rank_df...1....rank_df...2.ncol.rank_df..          0
### 2 rank_df...1....rank_df...2.ncol.rank_df..          0
### 3 rank_df...1....rank_df...2.ncol.rank_df..          0
### 4 rank_df...1....rank_df...2.ncol.rank_df..          0
###
### $sensitivity_data$nT_absolutechange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                                           Parameter
### 1 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 2 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 3 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 4 as.numeric.as.character.df...1......as.numeric.as.character.df...
###   AbsoluteChange
### 1         -508.0
### 2         -428.5
### 3          649.5
### 4          674.5
###
### $sensitivity_data$pExtirpation_rankchange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                   Parameter RankChange
### 1 rank_df...1....rank_df...2.ncol.rank_df..          0
### 2 rank_df...1....rank_df...2.ncol.rank_df..          0
### 3 rank_df...1....rank_df...2.ncol.rank_df..          0
### 4 rank_df...1....rank_df...2.ncol.rank_df..          0
###
### $sensitivity_data$pExtirpation_absolutechange
###              Scenario Increase_or_Decrease
### 1 Scenario: EffortA13                +0.15
### 2 Scenario: EffortA12                +0.15
### 3 Scenario: EffortA13                -0.15
### 4 Scenario: EffortA12                -0.15
###                                                           Parameter
### 1 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 2 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 3 as.numeric.as.character.df...1......as.numeric.as.character.df...
### 4 as.numeric.as.character.df...1......as.numeric.as.character.df...
###   AbsoluteChange
### 1              0
### 2              0
### 3              0
### 4              0
