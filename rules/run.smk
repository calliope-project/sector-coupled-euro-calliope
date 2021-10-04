
rule run_eurocalliope:
    message: "Running Calliope {wildcards.resolution} model with {wildcards.model_resolution} hourly temporal resolution for the model year {wildcards.year}"
    input:
        script = "src/run/run.py",
        model_yaml_path = "build/model/{resolution}/model-{year}.yaml"
    params:
        scenario = lambda wildcards: f"industry_fuel_shared,transport,heat,config_overrides,gas_storage,link_cap_dynamic,freeze-hydro-capacities,res_{wildcards.model_resolution}h,add-biofuel",
    envmodules: "gurobi/9.0.2"
    output: "build/{resolution}/outputs/run_{year}_{model_resolution}H.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/run/run.py"
