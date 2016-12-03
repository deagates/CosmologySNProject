# Runs "mcmc_nested_external.py"
# So multiple runs can be preformed changing mcmc run lengths to check convergence
#mcmc_nested_external.py cannot be run of its own accord


import mcmc_nested_external as SNe

single_run_length=3
single_run_length_mini=3
NUM_ITER = 1
ResultsName='DelilahResults/bez.txt'

SNe.main(NUM_ITER,single_run_length,single_run_length_mini,ResultsName)
