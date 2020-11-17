URL_NATIONAL_UNITS = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/NUTS_RG_20M_2016_4326_LEVL_0.geojson"

subworkflow landeligibility:
    workdir: "land-eligibility/"
    snakefile: "land-eligibility/Snakefile"
    configfile: "land-eligibility/config/default.yaml"

localrules: maps

rule maps:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = "src/analyse/maps.py",
        model = "inputs/{resolution}/model.nc",
        units = eurocalliope("build/data/eurospores/units.geojson")
    params:
        bounds = config["scope"]["bounds"]
    output:
        "inputs/{resolution}/map.pdf"
    conda: "../envs/plots.yaml"
    script: "../src/analyse/maps.py"


rule input_netcdf:
    message: "Creating Calliope {wildcards.resolution} resolution NetCDF model"
    input:
        src = "src/analyse/run.py",
        model_yaml_path = "build/model/{resolution}/model.yaml"
    params:
        scenario = "directional-rooftop-pv,industry_fuel,transport,heat,config_overrides",
        run = False
    output:
        output_model_path = "inputs/{resolution}/model.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"


rule simplified_nuts_units:
    message: "Download low resolution NUTS0 units as geojson."
    output:
        "build/figures/national_units.geojson"
    params:
        url = URL_NATIONAL_UNITS
    shell:
        "curl -sLo {output} '{params.url}'"


rule figure_1:
    message: "Plotting first paper figure"
    input:
        units=landeligibility("build/eurospores/units.geojson"),
        annual_demand="build/model/eurospores/annual-demand.csv",
        national_units=rules.simplified_nuts_units.output[0],
        electricity_demand="build/model/eurospores/electricity-demand.csv",
        space_heat_demand="build/model/eurospores/space-heat-demand.csv",
        water_heat_demand="build/model/eurospores/water-heat-demand.csv",
        cooking_demand="build/model/eurospores/cooking-demand.csv",
    params:
        energy_scaling_factor = config["scaling-factors"]["energy"],
        model_year = config["year"]
    conda: "../envs/plots.yaml"
    output:
        out_path="build/figures/fig_1.pdf"
    script: "../src/analyse/paper_figs.py"
