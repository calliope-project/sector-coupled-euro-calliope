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
        model = "build/{resolution}/inputs/run_2018_1H.nc",
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
        spores = "build/eurospores/shared_spores",

    params:
        countries = config["scope"]["spatial"]["countries"],
        model_year = config["plot-year"]
    conda: "../envs/analyse.yaml"
    output: "build/figures/fig_2.{plot_suffix}"
    script: "../src/analyse/plot_gross_available_energy.py"


rule spores_result_metrics:
    message: "Get {wildcards.resolution} {wildcards.spore_dir} slack{wildcards.slack} {wildcards.spore}.nc high-level metrics"
    input:
        src = "src/analyse/spores_metrics.py",
        src_util = "src/analyse/visualisation_util.py",
        cost_opt_model = "build/{resolution}/spores_2h{slack}/{spore_dir}/spore_0.nc",
        spore = "build/{resolution}/spores_2h{slack}/{spore_dir}/{spore}.nc",
        technical_potential_area=landeligibility("build/{resolution}/technical-potential/areas.csv"),
        technical_potential_protected_area=landeligibility("build/{resolution}/technical-potential-protected/areas.csv")
    params: config = config
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/spores_2h{slack,.*}_metrics/{spore_dir}/{spore}")
    script: "../src/analyse/spores_metrics.py"


rule consolidate_spores_result_metrics:
    message: "Consolidate {wildcards.scenario} scenario, slack{wildcards.slack} {wildcards.resolution} SPORES into high-level metrics"
    input:
        src = "src/analyse/consolidate_spores.py",
        cost_optimal_model = lambda wildcards: "build/{0}/spores_2h_metrics/{1}-industry_fuel_{2},spores_electricity/spore_0".format(
            wildcards.resolution,
            1 if wildcards.scenario == "isolated" else 2,
            wildcards.scenario
        ),
        spores = lambda wildcards: expand(
            f"build/{wildcards.resolution}/spores_2h{wildcards.slack}_metrics/{{spore}}",
            spore=[
                i.split(f"spores_2h{wildcards.slack}/")[1].replace(".nc", "")
                for i in glob.glob(f"build/{wildcards.resolution}/spores_2h{wildcards.slack}/*{wildcards.scenario}*/*.nc")
                if "spore_0.nc" not in i
            ]
        )
    conda: "../envs/analyse.yaml"
    output: directory("build/{resolution}/{scenario}_spores{slack,.*}")
    script: "../src/analyse/consolidate_spores.py"
