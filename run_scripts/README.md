# Accompanying Calliope run scripts for the sector-coupled Euro-Calliope

This directory contains scripts to run batches of Calliope model runs on a high-performance computing cluster based on the LSF software.
These scripts are independent from the main workflow since they can take a very long time to run and sometimes take longer than a single cluster queue will allow (~5 days), requiring the run to be re-initiated part-way through.
These two elements make running the models unideal for use with snakemake.

The three subdirectories found here are used as follows:

1. `base_spores`: To queue the baseline SPORES runs for the Sector-Coupled Euro-Calliope, in which there are no explicit exclusions of technologies, but there are different groups of technologies to which "scoring" is applied.
You should first run `spores_2h.sh` (`bsub < spores_2h.sh`).
When these jobs have reached their time limit, you can then run each of the `spores_continue` runs, which will pick up where the previous runs left off.
In the case of `spores_continue_2h_5slack` and `spores_continue_2h_15slack` (which are cost relaxation sensitivity runs), the cost-optimal solution which is dropped out early on from running `spores_2h.sh` will be used as the starting point.

2. `exclusion_spores`: To queue the exclusion SPORES runs for the Sector-Coupled Euro-Calliope, in which spatial diversity of electricity supply technologies is sought after which explicitly trying to minimise the existence of specific technologies.
You should queue these runs only after running `spores_2h.sh` (and producing the relevant cost-optimal model file `./build/eurospores/spores_2h/16-industry_fuel_shared,spores_supply/spore_0.nc`).
The `demand_update` scripts are related to sensitivity runs with updated demands.
They expect the existence of the file `${DIR}/build/eurospores/spores_2h_demand_update/16-industry_fuel_shared,spores_supply/spore_0.nc`, which you will need to prepare separately, by running the workflow elsewhere with relevant configuration flags to use demand scaling and then building and running the cost-optimal model before feeding it back in.

3. `test_runs`: To run parameterisation tests for model temporal resolution, transmission capacity upper limits, and setting the value of the Gurobi solver paramter BarConvTol.
