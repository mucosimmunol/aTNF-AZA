The script `communityModelling.R` contains the R-code for microbial community metabolism simulations using Flux-Balance-Analysis (FBA). 

To run the script, the custom R-Package `MicrbiomeAGORA` needs to be installed; e.g. via:

```
R CMD INSTALL MicrobiomeAGORA_0.1.4.tar.gz
```

Please note, that this package has a couple of dependencies incl. the R-packages `data.table` (>= 1.9.6), `stringr` (>= 1.0.0), `sybil` (>= 2.0.4), and a working local installation of CPLEX Optimizer (Version 12.9, IBM). 