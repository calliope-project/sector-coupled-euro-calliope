#!/bin/sh
#BSUB -J spores_exclude_2h_shared_update_demand[1-31]
#BSUB -n 6
#BSUB -R "rusage[mem=120G, scratch=10000]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050/build/logs/spores_exclude_shared_demand_update/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050"

mkdir -p  "${DIR}/build/logs/spores_exclude_shared_demand_update"  "${DIR}/build/ehighways/spores_2h_demand_update/16-industry_fuel_shared,spores_supply"

cd $TMPDIR

sh "${DIR}/run_scripts/spores_exclude_demand_update_jobs.sh" ${LSB_JOBINDEX} "${DIR}/build/ehighways/spores_2h_demand_update/16-industry_fuel_shared,spores_supply/" "${DIR}/run_scripts" 0.1 "${DIR}/build/ehighways/spores_2h/16-industry_fuel_shared,spores_supply/spore_0.nc" "${DIR}/build/ehighways/demand_update.nc"
