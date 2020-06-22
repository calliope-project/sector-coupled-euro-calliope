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


rule annual_heat_demand_national:
    message: "Calculate the contribution of each technology to heat demand at a national resolution."
    input:
        population = "build/population_map.csv",
        soil_temp = landeligibility("build/capacityfactors/tsoil5.nc"),
        air_temp = landeligibility("build/capacityfactors/temperature.nc"),
        annual_consumption = "build/annual_heat_consumption.csv"
    params:
        heat_tech_params = config['parameters']['heat-techs']
    conda: "../envs/default.yaml"
    output:
        demand="build/annual_heat_demand.csv",
    script: "src/construct/annual_heat_demand.py"


rule hourly_heat_demand_national:
    message: "Calculate national heat demand at an hourly resolution."
    input:
        population = "build/population_map.csv",
        soil_temp = landeligibility("build/capacityfactors/tsoil5.nc"),
        air_temp = landeligibility("build/capacityfactors/temperature.nc"),
        wind_10m = landeligibility("build/capacityfactors/wind_10m.nc"),
        annual_demand = rules.annual_heat_demand_national.output,
        dwellings = "data/automatic/dwellings.tsv.gz"  # update when2heat hardcoded SFH:MFH ratio
    params: scaling_factor = config["scaling-factors"]["heat"]  # should this be here?
    conda: "../envs/default.yaml"
    output: "build/national-heat-demand.csv"
    script: "src/construct/heat_demand.py"


rule heat_demand:
    message: "Calculate {wildcard.resolution} heat demand at an hourly resolution."
    input:
        heat_demand = rules.hourly_heat_demand_national.output,
        population = landeligibility('build/{resolution}/population.csv')
    params: scaling_factor = config["scaling-factors"]["heat"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/heat-demand.csv"
    script: "src/construct/heat_demand.py"