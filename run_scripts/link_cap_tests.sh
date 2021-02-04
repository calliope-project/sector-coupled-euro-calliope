#!/bin/sh
#BSUB -J eurospores_link_cap[1-11]
#BSUB -n 6
#BSUB -R "rusage[mem=40G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/link_cap_test/log_%I.log

declare -a array=("link_cap_20" "link_cap_40" "link_cap_10_40" "link_cap_1x" "link_cap_2x" "link_cap_5x" "link_cap_10x" "link_cap_50x" "link_cap_100x" "link_cap_500x" "link_cap_dynamic")
DIR="/cluster/work/cpesm/brynmorp/euro-spores"
calliope run --scenario "industry_fuel,transport,heat,config_overrides,res_3h,gas_storage,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"