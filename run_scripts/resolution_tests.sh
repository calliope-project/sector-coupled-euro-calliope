#!/bin/sh
#BSUB -J eurospores_resolution_test[1-5]
#BSUB -n 6
#BSUB -R "rusage[mem=140G]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/resolution_test/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/euro-spores"
declare -a array=("link_cap_dynamic" "link_cap_dynamic,res_2h" "link_cap_dynamic,res_3h" "link_cap_dynamic,res_6h" "link_cap_dynamic,res_12h")
mkdir -p "${DIR}/outputs/resolution_test" "${DIR}/logs/resolution_test"
calliope run --scenario "industry_fuel_isolated,transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/resolution_test/${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"