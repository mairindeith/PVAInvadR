param_names <<- c(
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
    'R0', #!# SHOULD THIS BE V0???
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
    'UR', #
    'UA', #
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
  notes <<- c(
    'Species abbreviation',
    'Maximum age (only 1% of fish survive when unfished)',
    'Age at recruitment (years)',
    'Number of pre-recruit stanzas',
    'Length of simulation (number of years)',
    'Time-step length as a proportion of one year (<1 means more than 1 time step per year)',
    'Number of simulations',
    'Number of gears applied to recruited animals',
    'Time-step when removal starts for each pre-recruit stanza (separate stanza-specific parameters with ";")',
    'Time-step when removal starts for each gear used on recruited animals (separate gear-specific parameters with ";")',
    'Range of recruitment at equilibrium',
    'Recruitment compensation ratio (difference in juvenile survival between unfished and near-zero density)',
    'Proportion of mortality at equilibrium that is due to cannibalism',
    'von Bertalanffy growth parameter',
    'Fecundity multiplier on weight',
    'Minimum weight at maturity (as a proportion of maximum weight)',
    'Time-steps of the year when spawning occurs (exclusive of the upper value if 2 are given; i.e. a range of 1-2 indicates that spawning occurs in the first time step only)',
    'Maximum survival by stanza (separate stanza-specific parameters with ";")',
    'Stanza-specific maximum available habitat (separate stanza-specific parameters with ";")',
    'Initial number of vulnerable fish in the population',
    'Hyperstability parameter',
    'Age at which cannibalism begins',
    'Standard deviation in recruitment across years',
    'Proportion of pre-recruit animals removed by gear (proportion caught per unit effort)',
    'Proportion of recruited animals removed by gear (proportion caught per unit effort)',
    'Time step in the year that each gear is fished (i.e. 1 if there is 1 time step per year; 1 or 2 if there are 2 time steps per year)',
    'Effort expended by each pre-recruit capture gear',
    'Effort expended by each recruit capture gear',
    'Fixed cost associated with each pre-recruit capture gear (separate stanza-specific parameters with ";")',
    'Fixed cost associated with each recruit capture gear (separate gear-specific parameters with ";")',
    'Effort-based cost associated with each pre-recruit capture gear (separate stanza-specific parameters with ";")',
    'Effort-based cost associated with each recruit capture gear (separate gear-specific parameters with ";")',
    'Discount rate',
    'Generation time',
    #'Proportional to ascending slope of dome-shaped selectivity curve',
    #'Proportional to length at maximum selectivity',
    #'Rate of decline of dome-shaped selectivity (bounded between 0-1)',
    'Logistic ascending slope of removal gear for recruited animals (as a proportion of Linf)',
    'Ascending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)',
    'Logistic descending slope of removal gear for recruited animals (as a proportion of Linf)',
    'Descending length at 50% selectivity of removal gear for recruited animals (as a proportion of Linf)',
    'Value to initialize (default is NA)'
  )

  file_out <<- data.frame(
    Parameters = param_names,
    values = "",
    Description = notes
  )
