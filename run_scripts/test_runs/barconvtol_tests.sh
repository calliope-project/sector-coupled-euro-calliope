#!/bin/sh
#BSUB -J eurospores_barconvtol[1-8]
#BSUB -n 6
#BSUB -R "rusage[mem=80G]"
#BSUB -W 1440
#BSUB -r
#BSUB -o /cluster/work/cpesm/brynmorp/euro-spores/logs/barconvtol/log_%I.log

DIR="/cluster/work/cpesm/brynmorp/euro-spores"
mkdir -p "${DIR}/outputs/barconvtol" "${DIR}/logs/barconvtol"
calliope run --scenario "industry_fuel_isolated,transport,heat,config_overrides,res_2h,gas_storage,freeze-hydro-capacities,link_cap_dynamic" --override_dict="{'run.solver_options.BarConvTol': 1e-${LSB_JOBINDEX}}" --save_netcdf "${DIR}/outputs/barconvtol/1e-${LSB_JOBINDEX}.nc" "${DIR}/build/model/eurospores/model.yaml"