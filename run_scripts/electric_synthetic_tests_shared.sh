#!/bin/sh
#BSUB -J eurospores_synthetic_electric_shared[1-13]
#BSUB -n 6
#BSUB -R "rusage[mem=130G]"
#BSUB -W 7200
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/synthetic_electric_shared/log_%I.log

fuel="shared"

declare -a array=("transport,no_biofuel" "no_electric" "no_synthetic" "no_electric,no_biofuel" "no_synthetic,no_biofuel" "no_electric_transport" "no_synthetic_transport" "no_electric_transport,no_biofuel" "no_synthetic_transport,no_biofuel" "transport,no_electric_heat" "transport,no_synthetic_heat" "transport,no_electric_heat,no_biofuel" "transport,no_synthetic_heat,no_biofuel")

DIR="/cluster/work/cpesm/brynmorp/euro-spores"
mkdir -p "${DIR}/outputs/synthetic_electric_${fuel}" "${DIR}/logs/synthetic_electric_${fuel}"

calliope run --scenario "industry_fuel_${fuel},heat,config_overrides,gas_storage,link_cap_dynamic,freeze-hydro-capacities,${array[${LSB_JOBINDEX}-1]}" --save_netcdf "${DIR}/outputs/synthetic_electric_${fuel}/${array[${LSB_JOBINDEX}-1]}.nc" "${DIR}/build/model/eurospores/model.yaml"