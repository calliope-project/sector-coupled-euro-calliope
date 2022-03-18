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

wildcard_constraints:
    plot_suffix = "pdf|png"

rule maps:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = "src/analyse/plot_transmission_map.py",
        model = "build/{resolution}/inputs/run_2018_2H.nc",
        units = landeligibility("build/{resolution}/units.geojson")
    params:
        bounds = config["scope"]["spatial"]["bounds"]
    output: "build/figures/{resolution}/map.{plot_suffix}"
    conda: "../envs/analyse.yaml"
    script: "../src/analyse/plot_transmission_map.py"


rule simplified_nuts_units:
    message: "Download low resolution NUTS0 units as geojson."
    output:
        "build/figures/national_units.geojson"
    params:
        url = URL_NATIONAL_UNITS
    shell:
        "curl -sLo {output} '{params.url}'"


rule figure_1:
    message: "Plotting first paper figure based on eurospores regions"
    input:
        src = "src/analyse/plot_annual_energy_demand.py",
        units = landeligibility("build/eurospores/units.geojson"),
        annual_demand = "build/eurospores/annual-demand.csv",
        electricity_demand = "build/model/eurospores/electricity-demand.csv",
        current_electricity_demand = eurocalliope("build/model/eurospores/electricity-demand.csv"),
        space_heat_demand = "build/model/eurospores/space-heat-demand.csv",
        water_heat_demand = "build/model/eurospores/water-heat-demand.csv",
        cooking_demand = "build/model/eurospores/cooking-demand.csv",
    params:
        transport_efficiency = config["parameters"]["transport"]["efficiency"],
        energy_scaling_factor = config["scaling-factors"]["energy"],
        model_year = config["plot-year"]
    conda: "../envs/analyse.yaml"
    output: "build/figures/fig_1.{plot_suffix}"
    script: "../src/analyse/plot_annual_energy_demand.py"


rule figure_2:
    message: "Plotting second paper figure based on eurospores regions, shared scenario SPORES results"
    input:
        src = "src/analyse/plot_gross_available_energy.py",
        energy_balances = "build/annual_energy_balances.csv",
        spores = "build/eurospores/shared_spores_metrics",
    params:
        countries = config["scope"]["spatial"]["countries"],
        model_year = config["plot-year"]
    conda: "../envs/analyse.yaml"
    output: "build/figures/fig_2.{plot_suffix}"
    script: "../src/analyse/plot_gross_available_energy.py"


rule spores_result_metrics:
    message: "Get {wildcards.resolution} {wildcards.spore_dir} slack{wildcards.slack} {wildcards.spore}.nc high-level metrics"
    input:
        src = "src/analyse/result_to_friendly.py",
        input_model = "build/{resolution}/spores_2h{slack}/{spore_dir}/spore_0.nc",
        model_result = "build/{resolution}/spores_2h{slack}/{spore_dir}/{spore}.nc",
        annual_demand = "build/{resolution}/annual-demand.csv",
        industry_demand = "build/annual_industry_energy_demand.csv.old",
        cfs = lambda wildcards: glob.glob(f"build/model/{wildcards.resolution}/capacityfactors-*.csv")
    params:
        initial_keywords = ["SPORES", "Subnational"],
        name = "calliope-subnational-spores-2050-carbon-neutral",
        description = "Sector-coupled Euro-Calliope SPORES outputs",
        year = 2018,
        scenario = 1,
        scenario_dim_name = "spore"
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/spores_2h{slack,.*}_metrics/{spore_dir}/{spore}")
    script: "../src/analyse/result_to_friendly.py"


rule cost_opt_result_metrics:
    message: "Get {wildcards.resolution} {wildcards.model_resolution}H resolution cost optimal model high-level metrics for weather year {wildcards.year}"
    input:
        src = "src/analyse/result_to_friendly.py",
        model_result = "build/{resolution}/outputs/run_{year}_{model_resolution}H.nc",
        input_model = "build/{resolution}/inputs/run_{year}_{model_resolution}H.nc",
        annual_demand = "build/{resolution}/annual-demand.csv",
        industry_demand = "build/annual_industry_energy_demand.csv.old",
        cfs = lambda wildcards: glob.glob(f"build/model/{wildcards.resolution}/capacityfactors-*.csv")
    params:
        initial_keywords = ["cost-optimal", "Subnational"],
        name = "calliope-subnational-cost-optimal-2050-carbon-neutral",
        description = "Sector-coupled Euro-Calliope Cost optimal outputs",
        year = lambda wildcards: wildcards.year,
        scenario = lambda wildcards: wildcards.year,
        scenario_dim_name = "weather_year"
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/cost_opt/{year}_{model_resolution}H")
    script: "../src/analyse/result_to_friendly.py"


rule consolidate_spores_result_metrics:
    message: "Consolidate {wildcards.scenario} scenario, slack{wildcards.slack} {wildcards.resolution} SPORES into high-level metrics"
    input:
        src = "src/analyse/consolidate_spores.py",
        all_friendly_files = lambda wildcards: expand(
            f"build/{wildcards.resolution}/spores_2h{wildcards.slack}_metrics/{{spore}}",
            spore=[
                i.split(f"spores_2h{wildcards.slack}/")[1].replace(".nc", "")
                for i in glob.glob(f"build/{wildcards.resolution}/spores_2h{wildcards.slack}/*{wildcards.scenario}*/*.nc")
                if "spore_0.nc" not in i
            ]
        )
    params:
        initial_keywords = ["SPORES", "Subnational"],
        name = "calliope-subnational-spores-2050-carbon-neutral",
        description = "Sector-coupled Euro-Calliope SPORES outputs",
    conda: "../envs/analyse.yaml"
    output:
        friendly_data = directory("build/{resolution}/{scenario}_spores{slack,.*}_metrics"),
        processed_spores = "build/{resolution}/{scenario}_spores{slack,.*}_list.csv"
    script: "../src/analyse/consolidate_spores.py"


rule consolidate_cost_opt_metrics:
    message: "Consolidate all years cost optimal {wildcards.resolution} {wildcards.model_resolution} results into high-level metrics"
    input:
        src = "src/analyse/consolidate_cost_opt.py",
        all_friendly_files = expand(
            "build/{{resolution}}/cost_opt/{year}_{{model_resolution}}H",
            year=range(config["scope"]["temporal"]["first-year"],
            config["scope"]["temporal"]["final-year"] + 1)
        ),
    params:
        initial_keywords = ["cost-optmal", "Subnational"],
        name = "calliope-subnational-cost-opt-2050-carbon-neutral",
        description = "Sector-coupled Euro-Calliope cost optimal outputs",
    conda: "../envs/analyse.yaml"
    output:
        friendly_data = directory("build/{resolution}/cost_opt_metrics_{model_resolution}H"),
    script: "../src/analyse/consolidate_cost_opt.py"


rule plot_map_metrics:
    message: "Plot {wildcards.spore} map metrics for interactive visualisation"
    input:
        src = "src/analyse/plot_spores_map_metrics.py",
        friendly_data = "build/eurospores/shared_spores",
        units = landeligibility("build/eurospores/units.geojson"),
        unit_groups = "data/plotting_unit_groups.csv",
    conda: "../envs/analyse.yaml"
    output: "build/figures/map_fig/{spore}.jpg"
    script: "../src/analyse/plot_spores_map_metrics.py"

rule plot_all_map_metrics:
    message: "Plot all SPORES map metrics for interactive visualisation"
    input: expand("build/figures/map_fig/{spore}.jpg", spore=range(1, 441))
