#!/bin/sh
#BSUB -J calliope[1-1]
#BSUB -n 4
#BSUB -R "rusage[mem=50G]"
#BSUB -W 240
#BSUB -r
#BSUB -o log_%I.log

./eurospores.sh.array.sh ${LSB_JOBINDEX}
