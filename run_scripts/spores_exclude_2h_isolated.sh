#!/bin/sh
#BSUB -J spores_exclude_2h_isolated[1-29]
#BSUB -n 6
#BSUB -R "rusage[mem=100G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/spores_exclude_isolated/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/euro-spores"

mkdir -p  "${DIR}/logs/spores_exclude_isolated"

sh "${DIR}/run_scripts/spores_exclude_jobs.sh" ${LSB_JOBINDEX} "${DIR}/outputs/spores_2h/15-industry_fuel_isolated,spores_supply/" "${DIR}/run_scripts"
