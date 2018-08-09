#!/bin/sh
export ARRAY=289
## for testing
qsub -A CPOX -n 1 -t 1:00:00 -q debug-cache-quad --mode script ~/Code/conformer/scripts/anl.optimizing_others.sh

## for actually running the job
##qsub -A CPOX -n $ARRAY -t 6:00:00 -q debug-cache-quad --mode script ~/Code/conformer/scripts/anl.optimizing_others.sh
##qsub -A CPOX -n 1 -t 3:00:00 --mode script ./nwchem.script.sh