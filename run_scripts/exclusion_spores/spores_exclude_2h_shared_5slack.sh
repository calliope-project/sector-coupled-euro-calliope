#!/bin/sh
#BSUB -J spores_exclude_2h_shared[1-31]
#BSUB -n 6
#BSUB -R "rusage[mem=120G, scratch=10000]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050/build/logs/spores_exclude_shared_5slack/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050"

mkdir -p  "${DIR}/build/logs/spores_exclude_shared_5slack"

cd $TMPDIR

sh "${DIR}/run_scripts/spores_exclude_jobs.sh" ${LSB_JOBINDEX} "${DIR}/build/ehighways/spores_2h_5slack/16-industry_fuel_shared,spores_supply/" "${DIR}/run_scripts" 0.05
