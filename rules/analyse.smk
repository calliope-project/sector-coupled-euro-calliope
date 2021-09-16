import glob

URL_NATIONAL_UNITS = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/NUTS_RG_20M_2016_4326_LEVL_0.geojson"

subworkflow landeligibility:
    workdir: "land-eligibility/"
    snakefile: "land-eligibility/Snakefile"
    configfile: "land-eligibility/config/default.yaml"
subworkflow eurocalliope:
    workdir: "euro-calliope/"
    snakefile: "euro-calliope/Snakefile"
    configfile: "config/default.yaml"

localrules: maps

rule maps:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = "src/analyse/maps.py",
        model = "build/{resolution}/model.nc",
        units = landeligibility("build/{resolution}/units.geojson")
    params:
        bounds = config["scope"]["bounds"]
    output: "build/figures/{resolution}/map.pdf"
    conda: "../envs/analyse.yaml"
    script: "../src/analyse/maps.py"


rule input_netcdf:
    message: "Creating Calliope {wildcards.resolution} resolution NetCDF model"
    input:
        src = "src/analyse/run.py",
        model_yaml_path = "build/model/{resolution}/model.yaml"
    params:
        scenario = "industry_fuel_isolated,transport,heat,config_overrides,gas_storage,link_cap_dynamic",
        run = False
    output: "build/{resolution}/model.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"


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


rule simplified_nuts_units:
    message: "Download low resolution NUTS0 units as geojson."
    output:
        "build/figures/national_units.geojson"
    params:
        url = URL_NATIONAL_UNITS
    shell:
        "curl -sLo {output} '{params.url}'"


rule figure_1:
    message: "Plotting first paper figure based on {wildcards.resolution} regions"
    input:
        units=landeligibility("build/{resolution}/units.geojson"),
        annual_demand="build/{resolution}/annual-demand.csv",
        national_units=rules.simplified_nuts_units.output[0],
        electricity_demand="build/model/{resolution}/electricity-demand.csv",
        current_electricity_demand=eurocalliope("build/model/{resolution}/electricity-demand.csv"),
        space_heat_demand="build/model/{resolution}/space-heat-demand.csv",
        water_heat_demand="build/model/{resolution}/water-heat-demand.csv",
        cooking_demand="build/model/{resolution}/cooking-demand.csv",
    params:
        energy_scaling_factor = config["scaling-factors"]["energy"],
        model_year = config["year"]
    conda: "../envs/analyse.yaml"
    output:
        out_path="build/figures/{resolution}/fig_1.pdf"
    script: "../src/analyse/paper_figs.py"


rule spores_result_metrics:
    message: "Get {wildcards.resolution} {wildcards.spore_dir} {wildcards.spore}.nc high-level metrics"
    input:
        src = "src/analyse/spores_metrics.py",
        src_util = "src/analyse/visualisation_util.py",
        cost_opt_model = "build/{resolution}/spores_2h/{spore_dir}/spore_0.nc",
        spore = "build/{resolution}/spores_2h/{spore_dir}/{spore}.nc",
        technical_potential_area=landeligibility("build/{resolution}/technical-potential/areas.csv"),
        technical_potential_protected_area=landeligibility("build/{resolution}/technical-potential-protected/areas.csv")
    params: config = config
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/spores_2h_metrics/{spore_dir}/{spore}")
    script: "../src/analyse/spores_metrics.py"


rule consolidate_spores_result_metrics:
    message: "Consolidate {wildcards.scenario} scenario {wildcards.resolution} SPORES into high-level metrics"
    input:
        src = "src/analyse/consolidate_spores.py",
        cost_optimal_model = lambda wildcards: "build/{0}/spores_2h_metrics/{1}-industry_fuel_{2},spores_electricity/spore_0".format(
            wildcards.resolution,
            1 if wildcards.scenario == "isolated" else 2,
            wildcards.scenario
        ),
        spores = lambda wildcards: expand(
            f"build/{wildcards.resolution}/spores_2h_metrics/{{spore}}",
            spore=[
                i.split("spores_2h/")[1].replace(".nc", "")
                for i in glob.glob(f"build/{wildcards.resolution}/spores_2h/*{wildcards.scenario}*/*.nc")
                if "spore_0.nc" not in i
            ]
        )
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/{scenario}_spores")
    script: "../src/analyse/consolidate_spores.py"
