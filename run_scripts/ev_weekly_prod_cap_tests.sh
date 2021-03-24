#!/bin/sh
#BSUB -J eurospores_weekly_ev_test[1-3]
#BSUB -n 6
#BSUB -R "rusage[mem=140G]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/weekly_ev_test/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/euro-spores"
declare -a array=("all_transport" "all_transport,weekly_transport_demand_range" "all_transport,weekly_transport_demand_equality")
mkdir -p "${DIR}/outputs/weekly_ev_test" "${DIR}/logs/weekly_ev_test"
calliope run --scenario "industry_fuel_isolated,annual_transport_distance,transport_demand_exists,transport_ev_exists,transport_ice_exists,heat,config_overrides,gas_storage,freeze-hydro-capacities,link_cap_dynamic,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/weekly_ev_test/${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"