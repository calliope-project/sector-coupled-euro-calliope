localrules: copy_euro_calliope, copy_resolution_specific_euro_calliope, model, eurostat_data_tsv, ch_data_xlsx, when2heat, download_weather_data, download_ramp_data

ruleorder: model > links > copy_from_template > copy_euro_calliope > annual_national_demand > annual_subnational_demand > heat_demand_profiles > cooking_heat_demand > scaled_heat_demand_profiles > scaled_public_transport_demand_profiles > update_electricity_with_other_sectors > heat_pump_characteristics > ev_energy_cap > annual_fuel_demand_constraints > annual_vehicle_constraints > annual_heat_constraints > gas_storage > copy_fuel_supply_techs > copy_biofuel_techs > copy_resolution_specific_euro_calliope
wildcard_constraints:
    definition_file = "[^\/]*", # must not traverse into directories
    year = "|".join([str(i) for i in range(2010, 2019)])

rule copy_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.{suffix}"),
    output: "build/model/{definition_file}.{suffix}"
    shell: "ln {input} {output}"


rule copy_resolution_specific_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input:
        eurocalliope("build/model/{resolution}/{definition_file}.{suffix}"),
    output: "build/model/{resolution}/{definition_file}.{suffix}"
    wildcard_constraints:
        definition_file = "((locations)|(directional-rooftop)|{})".format("|".join([
                "(capacityfactors-open-field-pv)", "(capacityfactors-rooftop-pv-n)",
                "(capacityfactors-rooftop-pv-e-w)", "(capacityfactors-rooftop-pv-s-flat)",
                "(capacityfactors-wind-offshore)", "(capacityfactors-wind-onshore)",
                "(capacityfactors-hydro-ror)", "(capacityfactors-hydro-reservoir-inflow)",
                "(capacityfactors-rooftop-pv)"
            ]))
    shell: "ln {input} {output}"


rule links:
    message: "Create links for {wildcards.resolution} resolution."
    input:
        src = "src/construct/template_links.py",
        gtc = "data/transmission.csv"
    params:
        scaling_factors = config["scaling-factors"],
        costs = config["parameters"]["transmission-costs"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/links.yaml"
    script: "../src/construct/template_links.py"


rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"][f"eurostat-{wildcards.dataset}"]
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.url}"


rule download_weather_data:
    message: "Download data/automatic/weather/{wildcards.filename}."
    params: url = lambda wildcards: config["data-sources"]["weather"].format(filename=wildcards.filename)
    output: protected("data/automatic/weather/{filename}")
    shell: "curl -sLo {output} '{params.url}'"


rule download_ramp_data:
    message: "Download data/automatic/ramp/{wildcards.ramp_profile}.csv.gz."
    params: url = lambda wildcards: config["data-sources"]["ramp-data"].format(ramp_profile=wildcards.ramp_profile)
    output: protected("data/automatic/ramp/{ramp_profile}.csv.gz")
    shell: "curl -sLo {output} '{params.url}'"


rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"][f"ch-{wildcards.dataset}"]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        src = "src/construct/annual_energy_balance.py",
        energy_balance = "data/automatic/eurostat-energy-balance.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/default.yaml"
    script: "../src/construct/annual_energy_balance.py"


rule annual_industry_demand:
    message: "Calculate {wildcards.projection_year} industry energy demand, following electrification and replacement of fossil feedstocks"
    input:
        src = "src/construct/annual_industry_demand_{projection_year}.py",
        energy_balances = rules.annual_energy_balances.output[0],
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv",
        jrc_industry_end_use = "data/industry/jrc_idees_processed_energy.csv.gz",
        jrc_industry_production = "data/industry/jrc_idees_processed_production.csv.gz",
        demand_scales = config["data-sources"]["industry-demand-scale"] if config["scale-demands"] else []
    params:
        scale_demand = config["scale-demands"]
    conda: "../envs/default.yaml"
    output:
        new_demand="build/annual_industry_energy_demand_{projection_year}.csv",
        bau_electricity="build/annual_industry_bau_electricity_{projection_year}.csv"
    wildcard_constraints:
        projection_year = "2030|2050"
    script: "../src/construct/annual_industry_demand_{wildcards.projection_year}.py"


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        src = "src/construct/annual_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_road_energy = "data/transport/jrc_idees_processed_road_energy.csv",
        jrc_road_distance = "data/transport/jrc_idees_processed_road_distance.csv",
        jrc_road_vehicles = "data/transport/jrc_idees_processed_road_vehicles.csv",
        jrc_rail_energy = "data/transport/jrc_idees_processed_rail_energy.csv",
        jrc_rail_distance = "data/transport/jrc_idees_processed_rail_distance.csv",
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/annual_road_transport_distance_demand.csv"),
        vehicles=temp("build/annual_road_transport_vehicles.csv"),
        efficiency=temp("build/annual_road_transport_efficiency.csv"),
        rail_energy=temp("build/annual_rail_transport_energy_demand.csv"),
        air_energy=temp("build/annual_air_transport_energy_demand.csv"),
        marine_energy=temp("build/annual_marine_transport_energy_demand.csv"),
        road_bau_electricity=temp("build/annual_road_transport_bau_electricity.csv"),
        rail_bau_electricity=temp("build/annual_rail_transport_bau_electricity.csv"),
    script: "../src/construct/annual_transport_demand.py"


rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        src = "src/construct/annual_heat_demand.py",
        hh_end_use = "data/automatic/eurostat-hh-end-use.tsv.gz",
        ch_end_use = "data/automatic/ch-end-use.xlsx",
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "data/commercial/jrc_idees_processed_energy.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    params:
        countries = config["scope"]["spatial"]["countries"],
        heat_tech_params = config["parameters"]["heat-end-use"]
    conda: "../envs/default.yaml"
    output:
        demand=temp("build/annual_heat_demand.csv"),
        electricity=temp("build/annual_heat_electricity_consumption.csv"),
    script: "../src/construct/annual_heat_demand.py"


rule annual_waste_supply:
    message: "Calculate the energy content of waste available in each {wildcards.resolution} region for incineration and energy recovery"
    input:
        energy_balance = rules.annual_energy_balances.output[0],
        population = eurocalliope("build/data/{resolution}/population.csv"),
        units = eurocalliope("build/data/{resolution}/units.geojson"),
    conda: "../envs/geodata.yaml"
    output:
        renewable = "build/{resolution}/annual_renewable_waste_supply.csv",
        total = "build/{resolution}/annual_waste_supply.csv"
    script: "../src/construct/annual_waste_supply.py"


rule annual_national_demand:
    message: "Scale {wildcards.projection_year} national demand to national resolution"
    input:
        src = "src/construct/annual_national_demand.py",
        annual_heat_demand = rules.annual_heat_demand.output.demand,
        annual_heat_electricity_consumption = rules.annual_heat_demand.output.electricity,
        industry_demand = rules.annual_industry_demand.output.new_demand,
        road_distance = rules.annual_transport_demand.output.distance,
        road_vehicles = rules.annual_transport_demand.output.vehicles,
        rail_demand = rules.annual_transport_demand.output.rail_energy,
        air_demand = rules.annual_transport_demand.output.air_energy,
        marine_demand = rules.annual_transport_demand.output.marine_energy,
        road_bau_electricity=rules.annual_transport_demand.output.road_bau_electricity,
        rail_bau_electricity=rules.annual_transport_demand.output.rail_bau_electricity,
        industry_bau_electricity=rules.annual_industry_demand.output.bau_electricity,
        demand_scales = config["data-sources"]["annual-demand-scale"] if config["scale-demands"] else []
    conda: "../envs/default.yaml"
    params:
        scaling_factors = config["scaling-factors"],
        industry_config = config["parameters"]["industry"],
        countries = config["scope"]["spatial"]["countries"],
        scale_demand = config["scale-demands"],
        demand_scale_scenario = config["demand-scale-scenario"]
    output: "build/national/annual-demand-{projection_year}.csv",
    script: "../src/construct/annual_national_demand.py"


rule annual_subnational_demand:
    message: "Scale {wildcards.projection_year} national demand to {wildcards.resolution} resolution"
    input:
        src = "src/construct/annual_subnational_demand.py",
        population = eurocalliope("build/data/{resolution}/population.csv"),
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        annual_heat_demand = rules.annual_heat_demand.output.demand,
        annual_heat_electricity_consumption = rules.annual_heat_demand.output.electricity,
        industry_demand = rules.annual_industry_demand.output.new_demand,
        road_distance = rules.annual_transport_demand.output.distance,
        road_vehicles = rules.annual_transport_demand.output.vehicles,
        rail_demand = rules.annual_transport_demand.output.rail_energy,
        air_demand = rules.annual_transport_demand.output.air_energy,
        marine_demand = rules.annual_transport_demand.output.marine_energy,
        road_bau_electricity=rules.annual_transport_demand.output.road_bau_electricity,
        rail_bau_electricity=rules.annual_transport_demand.output.rail_bau_electricity,
        industry_bau_electricity=rules.annual_industry_demand.output.bau_electricity,
        emissions = "data/industry/all_industrial_ets_eprtr_sites.geojson",
        freight = "data/automatic/eurostat-freight.tsv.gz",
        employees = "data/automatic/eurostat-employees.tsv.gz",
        gva = "data/automatic/eurostat-gva.tsv.gz",
        ch_gva = "data/automatic/ch-gva.xlsx",
        nuts_to_regions = "data/statistical_units_to_ehighways_regions.csv",
        industry_activity_codes = "data/industry/industry_activity_codes.csv",
        demand_scales = config["data-sources"]["annual-demand-scale"] if config["scale-demands"] else []
    conda: "../envs/geodata.yaml"
    params:
        scaling_factors = config["scaling-factors"],
        industry_config = config["parameters"]["industry"],
        scale_demand = config["scale-demands"],
        demand_scale_scenario = config["demand-scale-scenario"]
    output: "build/{resolution}/annual-demand-{projection_year}.csv",
    script: "../src/construct/annual_subnational_demand.py"


rule regions:
    message: "Get country codes corresponding to each {wildcards.resolution} region"
    input:
        src = "src/construct/regions.py",
        units = eurocalliope("build/data/{resolution}/units.geojson")
    conda: "../envs/geodata.yaml"
    output: "build/{resolution}/regions.csv"
    script: "../src/construct/regions.py"


rule raw_population_zipped:
    message: "Download population data."
    output:
        protected("data/automatic/raw-population-data.zip")
    params: url = config["data-sources"]["population"]
    conda: "../euro-calliope/envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule raw_population:
    message: "Extract population data TIF."
    input: rules.raw_population_zipped.output
    output: "build/population.tif"
    conda: "../euro-calliope/envs/shell.yaml"
    shell: "unzip {input} '*.tif' -d ./build/"


rule weather_and_population:
    message: "Determine weather conditions and population within each {wildcards.resolution} region, at a MERRA-2 gridcell resolution"
    input:
        src = "src/construct/weather.py",
        population = rules.raw_population.output[0],
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        air_temp = "data/automatic/weather/temperature.nc",
        wind_10m = "data/automatic/weather/wind10m.nc",
        soil_temp = "data/automatic/weather/tsoil5.nc",
    conda: "../envs/geodata.yaml"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"]
    output:
        weather_pop = "build/{resolution}/weather_pop.csv.gz",
    script: "../src/construct/weather.py"


rule when2heat:
    message: "Clone when2heat github repo"
    output: directory("data/automatic/when2heat")
    params: url = config["data-sources"]["when2heat"]
    shell:
        "git clone --depth 1 --branch 2019-08-06 {params.url} {output}"


rule heat_demand_profiles:
    message: "Generate unscaled hourly profiles at {wildcards.resolution} resolution."
    input:
        src = "src/construct/hourly_heat_profiles.py",
        dep_src = "src/construct/weather.py",
        weather_pop = rules.weather_and_population.output.weather_pop,
        when2heat = rules.when2heat.output[0],
    conda: "../envs/geodata.yaml"
    output:
        space_heat="build/{resolution}/space-heat-profile.csv",
        water_heat="build/{resolution}/water-heat-profile.csv",
        heat="build/{resolution}/heat-profile.csv",
    script: "../src/construct/hourly_heat_profiles.py"


rule regional_dwelling_ratio:
    message: "Get ratio of single family vs multi family homes in each {wildcards.resolution} region"
    input:
        src = "src/construct/dwellings.py",
        regions = rules.regions.output[0],
        dwellings = "data/automatic/eurostat-dwellings.tsv.gz",
        nuts_to_regions = "data/statistical_units_to_ehighways_regions.csv",
    conda: "../envs/geodata.yaml"
    output: "build/{resolution}/dwellings.csv"
    script: "../src/construct/dwellings.py"


rule cooking_heat_demand:
    message: "Clean RAMP-Cooking cooking{wildcards.demand_key}-demand {wildcards.resolution} profiles to match structure of other heat profiles"
    input:
        cooking_profiles = "data/automatic/ramp/cooking.csv.gz",
        regions = rules.regions.output[0],
        annual_demand = "build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"]),
    conda: "../envs/default.yaml"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        demand_key = "cooking{demand_key}"
    output: "build/model/{resolution}/cooking{demand_key,.*}-demand.csv"
    script: "../src/construct/scale_hourly_cooking_profiles.py"


rule scaled_heat_demand_profiles:
    message: "Scale {wildcards.end_use}heat{wildcards.demand_key} profiles at {wildcards.resolution} resolution according to annual demand."
    input:
        src = "src/construct/scale_hourly_heat_profiles.py",
        annual_demand = "build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"]),
        dwelling_ratio = rules.regional_dwelling_ratio.output[0],
        profile = "build/{resolution}/{end_use}heat-profile.csv",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        key = "{end_use,.*}heat{demand_key,.*}"  # ,.* allows the wildcard to be empty
    conda: "../envs/geodata.yaml"
    output:
        "build/model/{resolution}/{end_use,.*}heat{demand_key,.*}-demand.csv",
    script: "../src/construct/scale_hourly_heat_profiles.py"


rule scaled_public_transport_demand_profiles:
    message: "Scale hourly transport profiles at {wildcards.resolution} resolution according to annual demand."
    input:
        src = "src/construct/scale_hourly_transport_profiles.py",
        regions = rules.regions.output[0],
        annual_demand = "build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"]),
        rail_profiles = "data/transport/rail_daily_profiles_destinee.csv",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"]
    conda: "../envs/default.yaml"
    output: "build/{resolution}/public-transport-demand.csv"
    script: "../src/construct/scale_hourly_transport_profiles.py"


rule update_electricity_with_other_sectors:
    message: "Creating new electricity {wildcards.resolution} timeseries without heat demand"
    input:
        src = "src/construct/electricity_with_other_sectors.py",
        space_heat = "build/model/{resolution}/space-heat-bau-electricity-demand.csv",
        water_heat = "build/model/{resolution}/water-heat-bau-electricity-demand.csv",
        cooking = "build/model/{resolution}/cooking-bau-electricity-demand.csv",
        public_transport = "build/{resolution}/public-transport-demand.csv",
        annual_demand = "build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"]),
        hourly_electricity = eurocalliope("build/model/{resolution}/electricity-demand.csv"),
        regions = rules.regions.output[0],
        demand_scales = config["data-sources"]["building-electricity-demand-scale"] if config["scale-demands"] else []
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        scale_demand = config["scale-demands"],
        demand_scale_scenario = config["demand-scale-scenario"],
        projection_year = config["projection_year"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/electricity-demand.csv"
    script: "../src/construct/electricity_with_other_sectors.py"


rule heat_pump_characteristics:
    message: "Generate {wildcards.resolution} timeseries of {wildcards.tech}hp {wildcards.characteristic} for {wildcards.sink}heat carrier"
    input:
        src = "src/construct/heat_pump_characteristics.py",
        dep_src = "src/construct/hourly_heat_profiles.py",
        weather_pop = rules.weather_and_population.output.weather_pop,
        hp_characteristics = "data/heat_pump_characteristics.csv",
        annual_demand = "build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"])
    params:
        heat_tech_params = config["parameters"]["heat-end-use"],
        characteristic = "{characteristic}",
        tech = "{tech}hp",
        sink = "{sink}heat",
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/{characteristic}-{tech,.*}hp-{sink,.*}heat.csv"
    script: "../src/construct/heat_pump_characteristics.py"


rule ev_energy_cap:
    message: "Restructing RAMP-mobility EV {wildcards.dataset_name} profiles for use in {wildcards.resolution} Calliope model"
    input:
        src = "src/construct/hourly_ev_profiles.py",
        regions = rules.regions.output[0],
        ev_profiles = lambda wildcards: "data/automatic/ramp/ev-consumption.csv.gz" if "demand" in wildcards.dataset_name else f"data/automatic/ramp/ev-{wildcards.dataset_name}.csv.gz",
    params:
        demand_range = config["parameters"]["transport"]["monthly_demand"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/{dataset_name}-ev.csv"
    script: "../src/construct/hourly_ev_profiles.py"


rule annual_fuel_demand_constraints:
    message: "create all {wildcards.resolution} group constraints associated with annual fuel demand in the year {wildcards.year}"
    input:
        src = "src/construct/template_fuel_demand.py",
        annual_demand="build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"])
    params:
        scaling_factors = config["scaling-factors"],
        industry_carriers = config["parameters"]["industry"]["carriers"],
        model_time = config["calliope-parameters"]["model.subset_time"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/fuel_group_constraints_{year}.yaml"
    script: "../src/construct/template_fuel_demand.py"


rule annual_vehicle_constraints:
    message: "create all {wildcards.resolution} constraints associated with annual road vehicle demand in the year {wildcards.year}"
    input:
        src = "src/construct/template_vehicle_demand.py",
        annual_demand="build/{{resolution}}/annual-demand-{}.csv".format(config["projection_year"])
    params:
        transport = config["parameters"]["transport"],
        scaling_factors = config["scaling-factors"],
        model_time = config["calliope-parameters"]["model.subset_time"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/vehicle_group_constraints_{year}.yaml"
    script: "../src/construct/template_vehicle_demand.py"


rule annual_heat_constraints:
    message: "create all {wildcards.resolution} constraints associated with annual heat demand in the year {wildcards.year}"
    input:
        src = "src/construct/template_heat_demand.py",
        space_heat_demand = "build/model/{resolution}/space-heat-demand.csv",
        water_heat_demand = "build/model/{resolution}/water-heat-demand.csv",
        heat_demand = "build/model/{resolution}/heat-demand.csv",
        waste_supply = rules.annual_waste_supply.output.total
    params:
        storage_period = 48,  # there can only be as much storage as is reasonable for 48hrs of demand
        scaling_factors = config["scaling-factors"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/heat_group_constraints_{year}.yaml"
    script: "../src/construct/template_heat_demand.py"


rule gas_storage_xlsx:
    message: "Get data on underground gas storage available in each country"
    params: url = config["data-sources"]["gas-storage"]
    output: "data/automatic/gas_storage.xlsx"
    shell: "curl -sLo {output[0]} '{params.url}'"


rule gas_storage:
    message: "Assign gas storage facilities to {wildcards.resolution} regions"
    input:
        src = "src/construct/gas_storage.py",
        gas_storage_data = rules.gas_storage_xlsx.output[0],
        units = eurocalliope("build/data/{resolution}/units.geojson"),
    output:
        table = "build/{resolution}/gas_storage.csv",
        yaml = "build/model/{resolution}/gas_storage.yaml"
    params:
        scaling_factors = config["scaling-factors"]
    conda: "../envs/geodata.yaml"
    script: "../src/construct/gas_storage.py"


rule copy_from_template:
    message: "copy {wildcards.template} template"
    input:
        src = "src/construct/template_scenarios.py",
        template = "src/template/{template}"
    output: "build/model/{template}"
    params:
        shares = [i / 10 for i in range(11)],
        subset_time = config["calliope-parameters"]["model.subset_time"]
    wildcard_constraints:
        template = "((spores.yaml)|(config_overrides.yaml))"
    conda: "../envs/default.yaml"
    script: "../src/construct/template_scenarios.py"


rule fuel_cost_xlsx:
    message: "Download Heatroadmap fuel cost dataset"
    params: url = config["data-sources"]["fuel-costs"]
    output: "data/automatic/fuel_costs.xlsx"
    shell: "curl -sLo {output} {params.url}"


rule copy_fuel_supply_techs:
    message: "Build {wildcards.resolution} fuel supply YAML"
    input:
        src = "src/construct/template_fossil_fuel_supply.py",
        regions = rules.regions.output[0],
        fuel_costs = rules.fuel_cost_xlsx.output[0]
    output: "build/model/{resolution}/fossil-fuel-supply.yaml"
    params:
        scaling_factors = config["scaling-factors"],
        fuel_cost_source = config["parameters"]["fossil-fuel-cost"]["source"],
        fuel_cost_year = config["projection_year"]
    conda: "../envs/default.yaml"
    script: "../src/construct/template_fossil_fuel_supply.py"


rule copy_fuel_distribution_techs:
    message: "Build {wildcards.resolution} fuel distribution YAML"
    input:
        src = "src/construct/template_fuel_distribution_network.py",
        regions = rules.regions.output[0],
    output: "build/model/{resolution}/fuel-distribution.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/template_fuel_distribution_network.py"


rule copy_biofuel_techs:
    message: "Build {wildcards.resolution} biofuel supply YAML, subtracting {wildcards.year} renewable waste consumption to power CHPs"
    input:
        src = "src/construct/template_biofuel_supply.py",
        biofuel_potential = eurocalliope(
            "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv"
            .format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
        ),
        biofuel_costs = eurocalliope(
            "build/data/{{resolution}}/biofuel/{scenario}/costs-eur-per-mwh.csv"
            .format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
        ),
        renewable_waste_consumption_for_chp = rules.annual_waste_supply.output.renewable
    output: "build/model/{resolution}/biofuel-supply-{year}.yaml"
    params:
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    script: "../src/construct/template_biofuel_supply.py"


rule copy_2030_overrides:
    message: "Copy {wildcards.template} template describing 2030 overrides"
    input:
        src = "src/construct/template_2030.py",
        template = "src/template/2030/{template}"
    params:
        scaling_factors = config["scaling-factors"],
        heat = config["parameters"]["heat-end-use"],
    wildcard_constraints:
        template = "((heat-techs.yaml)|(renewable-techs.yaml)|(storage-techs.yaml)|(transformation-techs.yaml))"
    conda: "../envs/default.yaml"
    output: "build/model/overrides-2030/{template}"
    script: "../src/construct/template_2030.py"


rule emissions_scenario_yaml:
    message: "Generate Calliope {wildcards.resolution} emission target scenario YAML."
    input:
        src = "src/construct/template_emissions.py",
        emissions_targets = config["data-sources"]["emissions-targets"],
        regions = rules.regions.output[0]
    params:
        scaling_factors = config["scaling-factors"],
        projection_year = config["projection_year"],
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/emissions-scenarios.yaml"
    script: "../src/construct/template_emissions.py"


rule coal_supply_yaml:
    message: "Generate Calliope {wildcards.resolution} coal supply technology, with limit to today's capacities."
    input:
        src = "src/construct/template_coal_supply.py",
        power_plants = eurocalliope("data/automatic/JRC_OPEN_UNITS.csv"),
        units = eurocalliope("build/data/{resolution}/units.geojson"),
    params:
        scaling_factors = config["scaling-factors"],
    conda: "../envs/geodata.yaml"
    output: "build/model/{resolution}/coal_supply.yaml"
    script: "../src/construct/template_coal_supply.py"


rule model:
    message: "Build entire model on resolution {wildcards.resolution} for model year {wildcards.year}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/heat-techs.yaml",
        "build/model/demand-techs.yaml",
        "build/model/transformation-techs.yaml",
        "build/model/transport-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/legacy-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/directional-rooftop.yaml",
        "build/model/{resolution}/gas_storage.yaml",
        rules.copy_fuel_supply_techs.output,
        rules.copy_fuel_distribution_techs.output,
        rules.copy_biofuel_techs.output,
        rules.emissions_scenario_yaml.output,
        rules.coal_supply_yaml.output,
        rules.annual_fuel_demand_constraints.output,
        rules.annual_vehicle_constraints.output,
        rules.annual_heat_constraints.output,
        rules.links.output,
        "build/model/{resolution}/demand-equals-light-ev.csv",
        "build/model/{resolution}/demand-min-light-ev.csv",
        "build/model/{resolution}/demand-max-light-ev.csv",
        "build/model/{resolution}/demand-equals-heavy-ev.csv",
        "build/model/{resolution}/demand-min-heavy-ev.csv",
        "build/model/{resolution}/demand-max-heavy-ev.csv",
        "build/model/{resolution}/plugin-ev.csv",
        expand(
            "build/model/{{resolution}}/{characteristic}-{tech}-{sink}.csv",
            characteristic=["energy-cap", "cop"], tech=["ashp", "gshp", "hp"],
            sink=[config["parameters"]["heat-end-use"]["carriers"]] if isinstance(config["parameters"]["heat-end-use"]["carriers"], str) else config["parameters"]["heat-end-use"]["carriers"]
        ),
        expand(
            "build/model/{{resolution}}/{end_use}-demand.csv",
            end_use=["electricity", "cooking"] + ([config["parameters"]["heat-end-use"]["carriers"]] if isinstance(config["parameters"]["heat-end-use"]["carriers"], str) else config["parameters"]["heat-end-use"]["carriers"])
        ),
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=[
                "open-field-pv", "rooftop-pv-n", "rooftop-pv-e-w",
                "rooftop-pv-s-flat", "wind-offshore", "wind-onshore",
                "hydro-ror", "hydro-reservoir-inflow", "rooftop-pv"
            ],
        ),
        "build/model/config_overrides.yaml",
        "build/model/spores.yaml",
        "build/model/overrides-2030/heat-techs.yaml",
        "build/model/overrides-2030/renewable-techs.yaml",
        "build/model/overrides-2030/storage-techs.yaml",
        "build/model/overrides-2030/transformation-techs.yaml",
        template = "src/template/model.yaml",
        src = "src/construct/template_model_in_year.py"
    params:
        subset_time = config["calliope-parameters"]["model.subset_time"]
    output: "build/model/{resolution}/model-{year}.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/template_model_in_year.py"

rule model_all_years:
    message: "Build all years of model files"
    input: expand("build/model/{{resolution}}/model-{year}.yaml", year=range(config["scope"]["temporal"]["first-year"], config["scope"]["temporal"]["final-year"] + 1))
