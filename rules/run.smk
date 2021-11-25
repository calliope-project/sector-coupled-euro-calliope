
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

rule build_eurocalliope:
    message: "Building Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year}"
    input:
        script = "src/run/create_input.py",
        model_yaml_path = "build/model/{resolution}/model-{year}.yaml"
    params:
        scenario = lambda wildcards: f"industry_fuel_shared,transport,heat,config_overrides,gas_storage,link_cap_dynamic,freeze-hydro-capacities,res_{wildcards.model_resolution}h,add-biofuel",
    output: "build/{resolution}/inputs/run_{year}_{model_resolution}H.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/create_input.py"

rule run_eurocalliope:
    message: "Running Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year}"
    input:
        script = "src/run/run.py",
        model = rules.build_eurocalliope.output[0]
    envmodules: "gurobi/9.0.2"
    output: "build/{resolution}/outputs/run_{year}_{model_resolution}H.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/run.py"

rule build_storylines:
    message: "Building Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year} and {wildcards.storyline} social storyline"
    input:
        script = "src/run/create_input.py",
        model_yaml_path = "build/model/{resolution}/model-{year}.yaml"
    params:
        scenario = lambda wildcards: f"industry_fuel_shared,transport,heat,config_overrides,gas_storage,freeze-hydro-capacities,res_{wildcards.model_resolution}h,add-biofuel,{wildcards.storyline}-all",
    output: "build/{resolution}/inputs/run_{year}_{model_resolution}H_{storyline}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/create_input.py"

rule run_storylines:
    message: "Running Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year} and {wildcards.storyline} social storyline"
    input:
        script = "src/run/run_storylines.py",
        model = rules.build_storylines.output[0]
    params:
        max_transmission = lambda wildcards: config["storylines"]["max-transmission"][wildcards.storyline]
    envmodules: "gurobi/9.0.2"
    output: "build/{resolution}/outputs/run_{year}_{model_resolution}H_{storyline}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/run_storylines.py"
