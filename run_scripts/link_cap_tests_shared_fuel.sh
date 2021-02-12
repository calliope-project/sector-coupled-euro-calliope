#!/bin/sh
#BSUB -J eurospores_link_cap_shared[1-11]
#BSUB -n 6
#BSUB -R "rusage[mem=80G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/link_cap_test_shared/log_%I.log

fuel="shared"
declare -a array=("link_cap_20" "link_cap_40" "link_cap_10_40" "link_cap_1x" "link_cap_2x" "link_cap_5x" "link_cap_10x" "link_cap_50x" "link_cap_100x" "link_cap_500x" "link_cap_dynamic")
DIR="/cluster/work/cpesm/brynmorp/euro-spores"
mkdir -p "${DIR}/outputs/link_cap_test_${fuel}" "${DIR}/logs/link_cap_test_${fuel}"
calliope run --scenario "industry_fuel_${fuel},transport,heat,config_overrides,res_2h,gas_storage,freeze-hydro-capacities,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/link_cap_test_${fuel}/${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"