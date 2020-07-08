#' Save a .csv template of parameters. Once filled in, the csv template can be loaded using load_pva_parameters().
#'
#' @param filepath File path where the parameter template should be saved.
#' @return Nothing - a parameter template will be saved at the file location provided.
#' @examples
#' # Save the default template
#' pva_template('~/Documents/PVA_Examples/parameter_template.csv')

pva_template <- function(file){
  param_names <- c(
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
    notes <<- c(
      'Species abbreviation',
      'Maximum age (only 1% of fish survive when unfished)',
      'Age at recruitment (years)',
      'Number of pre-recruit stanzas',
      'Length of simulation (number of years)',
      'Time-step length as a proportion of one year (<1 means more than 1 time step per year)',
      'Number of simulations',
      'Number of gears applied to recruited animals',
      'Time-step when removal starts for each pre-recruit stanza (separate stanza-specific parameters using ;)',
      'Time-step when removal starts for each gear used on recruited animals (separate gear-specific parameters using ;)',
      'Range of unfished recruited abundance (separate lower and upper range using ;)',
      'Recruitment compensation ratio (difference in juvenile survival between unfished and near-zero density)',
      'Proportion of mortality at equilibrium that is due to cannibalism',
      'von Bertalanffy growth parameter',
      'Fecundity multiplier on weight',
      'Minimum weight at maturity (as a proportion of maximum weight)',
      'Time-steps of the year when spawning occurs (inclusive of the upper time step; i.e. a range of 1-2 indicates that spawning occurs in the first and second time steps)',
      'Maximum survival by stanza (separate stanza-specific parameters using ;)',
      'Stanza-specific maximum available habitat (separate stanza-specific parameters using ;)',
      'Initial number of vulnerable fish in the population',
      'Hyperstability parameter',
      'Age at which cannibalism begins',
      'Standard deviation in recruitment across years',
      'Proportion of pre-recruit animals removed by gear (proportion caught per unit effort; separate stanza-specific parameters using ;)',
      'Proportion of recruited animals removed by gear (proportion caught per unit effort; separate gear-specific parameters using ;)',
      'Time step(s) in the year that each gear is fished (enter "1" if there is only 1 time step per year; separate gear-specific parameters using ;)',
      'Effort expended by each pre-recruit capture gear (separate gear-specific parameters using ;)',
      'Effort expended by each recruit capture gear (separate gear-specific parameters using ;)',
      'Fixed cost associated with each pre-recruit capture gear (separate stanza-specific parameters using ;)',
      'Fixed cost associated with each recruit capture gear (separate gear-specific parameters using ;)',
      'Effort-based cost associated with each pre-recruit capture gear (separate stanza-specific parameters using ;)',
      'Effort-based cost associated with each recruit capture gear (separate gear-specific parameters using ;)',
      'Discount rate (optional)',
      'Generation time (optional)',
      #'Proportional to ascending slope of dome-shaped selectivity curve',
      #'Proportional to length at maximum selectivity',
      #'Rate of decline of dome-shaped selectivity (bounded between 0-1)',
      'Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)',
      'Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)',
      'Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)',
      'Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)',
      'Initial age structure in the population (default is NA; separate age-specific population sizes using ;)'
    )

    file_out <- data.frame(
      Parameters = param_names,
      Values = '',
      Description = notes
    )
    readr::write_csv(file_out, file)
}
