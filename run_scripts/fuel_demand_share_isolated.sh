#!/bin/sh
#BSUB -J eurospores_fuel_demand_share_isolated[1-11]
#BSUB -n 6
#BSUB -R "rusage[mem=130G]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/eurospores_fuel_demand_share_isolated/log_%I.log

fuel="isolated"
fuel_share="$((10 * (${LSB_JOBINDEX}-1)))"

DIR="/cluster/work/cpesm/brynmorp/euro-spores"
mkdir -p "${DIR}/outputs/eurospores_fuel_demand_share_${fuel}" "${DIR}/logs/eurospores_fuel_demand_share_${fuel}"

calliope run --scenario "industry_fuel_${fuel},transport,heat,config_overrides,gas_storage,link_cap_dynamic,freeze-hydro-capacities,demand_share_fuel_${fuel_share}" --save_netcdf "${DIR}/outputs/eurospores_fuel_demand_share_${fuel}/${fuel_share}.nc" "${DIR}/build/model/eurospores/model.yaml"