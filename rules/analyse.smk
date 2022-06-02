import glob

wildcard_constraints:
    plot_suffix = "pdf|png|tif"

rule transmission_map:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = script_dir + "analyse/plot_transmission_map.py",
        model = "build/inputs/{resolution}/run_2018_2H_neutral.nc",
        units = "build/data/{resolution}/units.geojson"
    params:
        bounds = config["euro-calliope"]["scope"]["spatial"]["bounds"]
    output: "build/figures/{resolution}/map.{plot_suffix}"
    conda: "../envs/analyse.yaml"
    script: "../src/analyse/plot_transmission_map.py"


rule cost_opt_result_metrics:
    message: "Get {wildcards.resolution} {wildcards.model_resolution}H resolution cost optimal model high-level metrics for weather year {wildcards.year}"
    input:
        src = script_dir + "analyse/result_to_friendly.py",
        model_result = "build/outputs/{resolution}/run_{year}_{model_resolution}H_{co2_scenario}.nc",
        input_model = "build/inputs/{resolution}/run_{year}_{model_resolution}H_{co2_scenario}.nc",
        annual_demand = "build/data/{resolution}/annual-demand-2050.csv",
        industry_demand = "build/annual_industry_energy_demand_2050.csv",
        cfs = lambda wildcards: glob.glob(f"build/model/{wildcards.resolution}/capacityfactors-*.csv")
    params:
        initial_keywords = ["cost-optimal", "Subnational"],
        name = lambda wildcards: f"calliope-{wildcards.resolution}-cost-opt-{wildcards.year}-{wildcards.co2_scenario}",
        description = "Sector-coupled Euro-Calliope Cost optimal outputs",
        year = lambda wildcards: wildcards.year,
        scenario = lambda wildcards: wildcards.year,
        scenario_dim_name = "weather_year"
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/cost_opt/{year}_{model_resolution}H_{co2_scenario}")
    script: "../src/analyse/result_to_friendly.py"


rule consolidate_cost_opt_metrics:
    message: "Consolidate all years cost optimal {wildcards.resolution} {wildcards.model_resolution} results into high-level metrics"
    input:
        src = script_dir + "analyse/consolidate_cost_opt.py",
        all_friendly_files = expand(
            "build/{{resolution}}/cost_opt/{year}_{{model_resolution}}H_{{co2_scenario}}",
            year=range(config["euro-calliope"]["scope"]["temporal"]["first-year"],
            config["euro-calliope"]["scope"]["temporal"]["final-year"] + 1)
        ),
    params:
        initial_keywords = ["cost-optimal", "Subnational"],
        name = lambda wildcards: f"calliope-{wildcards.resolution}-cost-opt-2050-{wildcards.co2_scenario}",
        description = "Sector-coupled Euro-Calliope cost optimal outputs",
    conda: "../envs/analyse.yaml"
    output:
        friendly_data = directory("build/{resolution}/cost_opt_metrics_{model_resolution}H_{co2_scenario}"),
    script: "../src/analyse/consolidate_cost_opt.py"
