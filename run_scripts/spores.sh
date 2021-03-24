#!/bin/sh
#BSUB -J eurospores_spores[1-14]
#BSUB -n 6
#BSUB -R "rusage[mem=180G]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/spores_full/log_%I.log


DIR="/cluster/work/cpesm/brynmorp/euro-spores"

declare -a array=("industry_fuel_isolated,spores_electricity" "industry_fuel_shared,spores_electricity" "industry_fuel_isolated,spores_storage" "industry_fuel_shared,spores_storage" "industry_fuel_isolated,spores_transmission" "industry_fuel_shared,spores_transmission" "industry_fuel_isolated,spores_fuel" "industry_fuel_shared,spores_fuel" "industry_fuel_isolated,spores_transport" "industry_fuel_shared,spores_transport" "industry_fuel_isolated,spores_heat" "industry_fuel_shared,spores_heat" "industry_fuel_isolated,spores_all" "industry_fuel_shared,spores_all")

mkdir -p "${DIR}/outputs/spores_full" "${DIR}/logs/spores_full" "${DIR}/outputs/spores_full/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}"

calliope run --scenario "transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,link_cap_dynamic,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/spores_full/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"