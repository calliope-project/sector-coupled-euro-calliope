#!/bin/sh
#BSUB -J eurospores_spores[1-16]
#BSUB -n 6
#BSUB -R "rusage[mem=100G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050/build/logs/spores_2h/log_%I.log


DIR="/cluster/work/cpesm/brynmorp/sentinel-free-model-runs/2050"

declare -a array=("industry_fuel_isolated,spores_electricity" "industry_fuel_shared,spores_electricity"  "industry_fuel_isolated,spores_storage" "industry_fuel_shared,spores_storage" "industry_fuel_isolated,spores_transmission" "industry_fuel_shared,spores_transmission" "industry_fuel_isolated,spores_fuel" "industry_fuel_shared,spores_fuel" "industry_fuel_isolated,spores_transport" "industry_fuel_shared,spores_transport" "industry_fuel_isolated,spores_heat" "industry_fuel_shared,spores_heat" "industry_fuel_isolated,spores_all" "industry_fuel_shared,spores_all" "industry_fuel_isolated,spores_supply" "industry_fuel_shared,spores_supply")

mkdir -p "${DIR}/build/eurospores/spores_2h" "${DIR}/build/logs/spores_2h" "${DIR}/build/eurospores/spores_2h/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}"

calliope run --scenario "transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,link_cap_dynamic,res_2h,add-biofuel,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/build/eurospores/spores_2h/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"
