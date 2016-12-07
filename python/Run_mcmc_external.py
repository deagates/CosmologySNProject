# Runs "mcmc_nested_external.py"
# So multiple runs can be preformed changing mcmc run lengths to check convergence
#mcmc_nested_external.py cannot be run of its own accord


import mcmc as SNe

single_run_length=10
NUM_ITER = 3
ResultsName='results_test2.txt'

SNe.main(NUM_ITER,single_run_length,single_run_length_mini,ResultsName)
