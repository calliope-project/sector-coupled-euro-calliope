URL_WHEN2HEAT = "https://github.com/oruhnau/when2heat"  # should this point to a commit hash?

subworkflow landeligibility:
    workdir: "./land-eligibility"
    snakefile: "./land-eligibility/Snakefile"
    configfile: "./land-eligibility/config/default.yaml"

# Process is:
# energy consumption for heating -> annual national heat demand (technology efficiencies)
# annual national heat demand -> hourly national heat demand (when2heat method)
# hourly national heat demand -> arbritary shape heat demand (based on population)

#rule cop:  # maybe this isn't necessary, if direct electric consumption = 100% efficiency in producing heat
#    message:"Calculate the annual average COP of heat pumps in each country, and the relative contribution of air source and ground source heat pumps."
#    input:
#        population = "build/population_map.csv",
#        soil_temp = landeligibility("build/capacityfactors/tsoil5.nc"),
#        air_temp = landeligibility("build/capacityfactors/temperature.nc"),
#        annual_consumption = "build/annual_heat_consumption.csv"
#    params:
#        heat_tech_params = config['parameters']['heat-techs']
#    conda: "../envs/default.yaml"
#    output:
#        cop="build/cop.csv",
#        annual_consumption="build/annual_heat_consumption_inc_hp.csv"
#    script: "src/construct/cop.py"


rule annual_heat_demand:
    message: "Calculate national heat demand, normalised based on population"
    input:
        hh_end_use = "data/automatic/hh_end_use.tsv.gz",
        ch_end_use = "data/automatic/ch_hh_end_use.xlsx",
        energy_balance = "build/annual_energy_balances.csv",
        jrc_industry_end_use = "data/industry/JRC_IDEES_industry_end_use_consumption.csv",
        jrc_commercial_end_use = "data/commercial/JRC_IDEES_commercial_end_use_consumption.csv"
        carrier_names = "data/energy_balance_carrier_names.csv"
        population = landeligibility("build/national/population.csv")
    params:
        countries = config["scope"]["countries"],
        heat_tech_params = config["parameters"]["heat-end-use"]
    conda: "../envs/default.yaml"
    output:
        demand=temp("build/annual_heat_demand.csv"),
        electricity=temp("build/annual_heat_electricity_consumption.csv"),
    script: "../src/construct/annual_heat_demand.py"


rule when2heat:
    message: "Clone when2heat github repo"
    output: "data/automatic/when2heat/"
    shell:
        "git clone {URL_WHEN2HEAT} {output}"


rule hourly_heat_demand:  # TODO: how to handle cooking demand? TODO: scale demand from TWh to whatever 'scaling factor' expects
    message: "Calculate {wildcard.resolution} heat demand at an hourly resolution."
    input:
        population = landeligibility("build/population-europe.tif"),
        units = landeligibility("build/{resolution}/units.geojson"),
        air_temp = "data/weather/temperature.nc",
        wind_10m = "data/weather/wind10m.nc",
        annual_demand = rules.annual_heat_demand.output.demand,
        when2heat = "data/automatic/when2heat/",
        dwellings = "data/automatic/dwellings.tsv.gz"
    params:
        scaling_factor = config["scaling-factors"]["heat"],  # should this be here?
        model_year = config["year"]
    conda: "../envs/default.yaml"
    output:
        space_heat="build/model/{resolution}/space-heat-demand.csv",
        water_heat="build/model/{resolution}/water-heat-demand.csv"
    script: "../src/construct/hourly_heat_demand.py"