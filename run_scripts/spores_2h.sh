#!/bin/sh
#BSUB -J eurospores_spores[1-16]
#BSUB -n 6
#BSUB -R "rusage[mem=100G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/spores_2h/log_%I.log


DIR="/cluster/work/cpesm/brynmorp/euro-spores"

declare -a array=("industry_fuel_isolated,spores_electricity" "industry_fuel_shared,spores_electricity"  "industry_fuel_isolated,spores_storage" "industry_fuel_shared,spores_storage" "industry_fuel_isolated,spores_transmission" "industry_fuel_shared,spores_transmission" "industry_fuel_isolated,spores_fuel" "industry_fuel_shared,spores_fuel" "industry_fuel_isolated,spores_transport" "industry_fuel_shared,spores_transport" "industry_fuel_isolated,spores_heat" "industry_fuel_shared,spores_heat" "industry_fuel_isolated,spores_all" "industry_fuel_shared,spores_all" "industry_fuel_isolated,spores_supply" "industry_fuel_shared,spores_supply")

mkdir -p "${DIR}/outputs/spores_2h" "${DIR}/logs/spores_2h" "${DIR}/outputs/spores_2h/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}"

calliope run --scenario "transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,link_cap_dynamic,res_2h,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/spores_2h/${LSB_JOBINDEX}-${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"
