# Lazy temporary script for testing PVA functions individually.
# To be run in the root directory of the package.

library(devtools)

load_all() # Load the PVAInavsR package
inputs = load_pva_parameters("../ShinyDropbox/TestInputs/DouglasYP_population_parameters_v0.csv")

# Jump right in with a test?

pvatest = PVA(inputs)
pvaplotted = plot_pva(pvatest)
