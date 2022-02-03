
wildcard_constraints:
    co2_scenario = "neutral|current"

rule spores:
    message: "Running Calliope {wildcards.resolution} resolution SPORES model"
    input:
        model_yaml_path = "build/model/{resolution}/model.yaml"
    params:
        scenario = "industry_fuel_isolated,transport,heat,config_overrides,gas_storage,link_cap_dynamic,spores_electricity",
        threads = 6,
        mins = 1440,
        mem = "rusage[mem=90G]",
    output:
        model = "outputs/{resolution}/spores/electricity.nc",
        logs = "logs/{resolution}/spores/electricity.log"
    conda: "../envs/calliope.yaml"
    shell:
        "bsub -n {params.threads} -W {params.mins} -R {params.mem} -o {output.logs} 'calliope run --save_netcdf {output.model} --scenario {params.scenario} {input.model_yaml_path}'"

def get_scenario(base_scenario, projection_year, co2_scenario, add_link_cap=True):
    if int(projection_year) == 2030:
        scenario = base_scenario + f",heat_techs_2030,renewable_techs_2030,transformation_techs_2030,2030_{co2_scenario},coal_supply,fossil-fuel-supply"
        if add_link_cap:
            scenario += ",link_cap_1x"
    elif int(projection_year) == 2050:
        if add_link_cap:
            scenario = base_scenario + ",link_cap_dynamic"
        else:
            scenario = base_scenario
        if co2_scenario == "current":
            scenario += f",2050_{co2_scenario},fossil-fuel-supply"
    return scenario


rule build_eurocalliope:
    message: "Building Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year} and with co2 `{wildcards.co2_scenario}` limit"
    input:
        script = "src/run/create_input.py",
        model_yaml_path = "build/model/{resolution}/model-{year}.yaml"
    params:
        scenario = lambda wildcards: get_scenario(
            f"industry_fuel_isolated,transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,res_{wildcards.model_resolution}h,add-biofuel,synfuel_transmission",
            config["projection-year"], wildcards.co2_scenario, add_link_cap=True
        )
    output: "build/{resolution}/inputs/run_{year}_{model_resolution}H_{co2_scenario}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/create_input.py"

rule run_eurocalliope:
    message: "Running Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year} and with co2 `{wildcards.co2_scenario}` limit"
    input:
        script = "src/run/run.py",
        model = rules.build_eurocalliope.output[0]
    envmodules: "gurobi/9.1.1"
    output: "build/{resolution}/outputs/run_{year}_{model_resolution}H_{co2_scenario}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/run.py"
