# Simulation code for "Linear Mixed Models for Comparing Dynamic Treatment Regimens on a Longitudinal Outcome in Sequentially Randomized Trials"


Command line arguments for each simulation script (`sim1.R`, `sim2.R`, `sim3.R`): 

1. sample size
2. number of Monte Carlo replicates
3. seed for the random number generator
4. desired number of cores for parallelization
5. effect size (optional, either "small", "med", or "large")

## Simulation 1
 
To obtain raw simulation results in `.RData` format:  
```
Rscript --vanilla sim1.R 50 1000 304 3
Rscript --vanilla sim1.R 200 1000 304 3
Rscript --vanilla sim1.R 1000 1000 304 3
Rscript --vanilla sim1.R 5000 1000 304 3 small
Rscript --vanilla sim1.R 5000 1000 304 3 med
Rscript --vanilla sim1.R 5000 1000 304 3 large
```
    
To produce LaTeX tables:

```
Rscript --vanilla fmt-sim1.R
```

## Simulation 2

To obtain raw simulation results in `.RData` format:  
```
Rscript --vanilla sim2.R 50 1000 304 3
Rscript --vanilla sim2.R 200 1000 304 3
Rscript --vanilla sim2.R 1000 1000 304 3
Rscript --vanilla sim2.R 5000 1000 304 3 small
Rscript --vanilla sim2.R 5000 1000 304 3 med
Rscript --vanilla sim2.R 5000 1000 304 3 large
```
    
To produce LaTeX tables:

```
Rscript --vanilla fmt-sim2.R
```

## Simulation 3 (in supplementary materials)

To obtain raw simulation results in `.RData` format:  
```
Rscript --vanilla sim3.R 50 1000 304 3
Rscript --vanilla sim3.R 200 1000 304 3
Rscript --vanilla sim3.R 1000 1000 304 3
Rscript --vanilla sim3.R 5000 1000 304 3 small
Rscript --vanilla sim3.R 5000 1000 304 3 med
Rscript --vanilla sim3.R 5000 1000 304 3 large
```
    
To produce LaTeX tables:

```
Rscript --vanilla fmt-sim3.R
```

