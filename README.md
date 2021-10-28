# causalCSME
R package implementing methods that adjust for both confounding and additive measurement error based on conditional score functions

### Replicating Paper Simulations

To replicate the first simulation study, run the simulator() function in the sims_scen1.R file 2000 times on a computing cluster or locally nested within a for loop or apply function procedure. To match the results in the paper exactly, set the seed in R such that the seed is set to 1000s for simulation number s. To replicate the second and third simulation study and three Web appendix simulations, repeat the procedure using simulator() in each of the respective sims_scen*.R and sims_appendix*.R files.

### Replicating Paper Application Section Analysis

To replicate the application section results, run the file application2.R. If this R package is loaded, it should pull the data correctly on its own; if the package is not loaded the dataset can be found at https://atlas.scharp.org/cpas/project/HVTN\%20Public\%20Data/HVTN\%20505/begin.view?.
