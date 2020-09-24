
URL_ENERGY_BALANCE = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/nrg_bal_c.tsv.gz"
URL_HH_END_USE = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/nrg_d_hhq.tsv.gz"
URL_CH_END_USE = "https://www.bfe.admin.ch/bfe/en/home/versorgung/statistik-und-geodaten/energiestatistiken/energieverbrauch-nach-verwendungszweck.exturl.html/aHR0cHM6Ly9wdWJkYi5iZmUuYWRtaW4uY2gvZGUvcHVibGljYX/Rpb24vZG93bmxvYWQvOTg1NA==.html"
URL_CH_ENERGY_BALANCE = "https://www.bfe.admin.ch/bfe/en/home/versorgung/statistik-und-geodaten/energiestatistiken/gesamtenergiestatistik.exturl.html/aHR0cHM6Ly9wdWJkYi5iZmUuYWRtaW4uY2gvZGUvcHVibGljYX/Rpb24vZG93bmxvYWQvNzUxOQ==.html"
URL_CH_INDUSTRY_ENERGY_BALANCE = "https://www.bfe.admin.ch/bfe/en/home/versorgung/statistik-und-geodaten/energiestatistiken/teilstatistiken.exturl.html/aHR0cHM6Ly9wdWJkYi5iZmUuYWRtaW4uY2gvZGUvcHVibGljYX/Rpb24vZG93bmxvYWQvODc4OA==.html"
URL_DWELLINGS = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/cens_11dwob_r3.tsv.gz"
URL_WHEN2HEAT = "https://github.com/oruhnau/when2heat"  # should this point to a commit hash?
URL_FREIGHT = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/road_go_na_rl3g.tsv.gz"
URL_EMPLOYEES = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/sbs_r_nuts06_r2.tsv.gz"
URL_GVA = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/nama_10r_3gva.tsv.gz"
URL_CH_GVA = "https://www.bfs.admin.ch/bfsstatic/dam/assets/10647597/master"

subworkflow eurocalliope:
    workdir: "euro-calliope"
    snakefile: "euro-calliope/Snakefile"
    configfile: "config/default.yaml"

subworkflow landeligibility:
    workdir: "land-eligibility/"
    snakefile: "land-eligibility/Snakefile"
    configfile: "land-eligibility/config/default.yaml"

localrules: copy_euro_calliope, copy_resolution_specific_euro_calliope, model, links, outer_countries, eurostat_data_tsv, ch_data_xlsx, when2heat
ruleorder: model > links > outer_countries > copy_euro_calliope > annual_subnational_demand > heat_demand_profiles > scaled_heat_demand_profiles > scaled_bau_electricity_heat_demand_profiles > scaled_public_transport_demand_profiles > update_electricity_with_other_sectors > copy_resolution_specific_euro_calliope
wildcard_constraints:
    definition_file = "[^\/]*" # must not travers into directories


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
    shell: "ln {input} {output}"


rule links:
    message: "Create links for {wildcards.resolution} resolution."
    input:
        src = "../src/construct/links.py",
        gtc = "data/eurospores.xlsx"
    params: scaling_factor = config["scaling-factors"]["power"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/links.yaml"
    script: "../src/construct/links.py"


rule outer_countries:
    message: "Create neighbouring country links {wildcards.resolution} scenario."
    input:
        src = "../src/construct/outer_countries.py",
        gtc = "data/eurospores.xlsx"
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/outer-countries.yaml"
    script: "../src/construct/outer_countries.py"


rule eurostat_data_tsv:
    message: "Get various datasets from Eurostat"
    output:
        energy_balance = protected("data/automatic/annual_energy_balances.tsv.gz"),
        hh_end_use = protected("data/automatic/hh_end_use.tsv.gz"),
        freight = protected("data/automatic/freight.tsv.gz"),
        employees = protected("data/automatic/employees.tsv.gz"),
        gva = protected("data/automatic/gva.tsv.gz"),
        dwellings = protected("data/automatic/dwellings.tsv.gz"),
    shell:
        """
        curl -sLo {output.energy_balance} '{URL_ENERGY_BALANCE}'
        curl -sLo {output.hh_end_use} '{URL_HH_END_USE}'
        curl -sLo {output.freight} '{URL_FREIGHT}'
        curl -sLo {output.employees} '{URL_EMPLOYEES}'
        curl -sLo {output.gva} '{URL_GVA}'
        curl -sLo {output.dwellings} '{URL_DWELLINGS}'
        """


rule ch_data_xlsx:
    message: "Get Swiss annual energy balances and household end uses"
    output:
        energy_balance = protected("data/automatic/ch_annual_energy_balances.xlsx"),
        industry_energy_balance = protected("data/automatic/ch_annual_industry_energy_balances.xlsx"),
        end_use = protected("data/automatic/ch_hh_end_use.xlsx"),
        gva = protected("data/automatic/ch_gva.xlsx")
    shell:
        """
        curl -sLo {output.energy_balance} '{URL_CH_ENERGY_BALANCE}'
        curl -sLo {output.industry_energy_balance} '{URL_CH_INDUSTRY_ENERGY_BALANCE}'
        curl -sLo {output.end_use} '{URL_CH_END_USE}'
        curl -sLo {output.gva} '{URL_CH_GVA}'
        """


rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        energy_balance = rules.eurostat_data_tsv.output.energy_balance,
        ch_energy_balance = rules.ch_data_xlsx.output.energy_balance,
        ch_industry_energy_balance = rules.ch_data_xlsx.output.industry_energy_balance,
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    shadow: "minimal"
    script: "../src/construct/annual_energy_balance.py"


rule annual_industry_demand:
    message: "Calculate future industry energy demand, following electrification and replacement of fossil feedstocks"
    input:
        energy_balances = rules.annual_energy_balances.output[0],
        cat_names = "data/energy_balance_category_names.csv",
        jrc_industry_end_use = "data/industry/jrc_idees_processed_energy.csv.gz",
        jrc_industry_production = "data/industry/jrc_idees_processed_production.csv.gz",
    conda: "../envs/default.yaml"
    output:
        new_demand=temp("build/annual_industry_energy_demand.csv"),
        bau_electricity=temp("build/annual_industry_bau_electricity.csv")
    script: "../src/construct/annual_industry_demand.py"


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_road_energy = "data/transport/jrc_idees_processed_road_energy.csv",
        jrc_road_distance = "data/transport/jrc_idees_processed_road_distance.csv",
        jrc_rail_energy = "data/transport/jrc_idees_processed_rail_energy.csv",
        jrc_rail_distance = "data/transport/jrc_idees_processed_rail_distance.csv",
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/annual_road_transport_distance_demand.csv"),
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
        hh_end_use = rules.eurostat_data_tsv.output.hh_end_use,
        ch_end_use = rules.ch_data_xlsx.output.end_use,
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "data/commercial/jrc_idees_processed_energy.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    params:
        countries = config["scope"]["countries"],
        heat_tech_params = config["parameters"]["heat-end-use"]
    conda: "../envs/default.yaml"
    output:
        demand=temp("build/annual_heat_demand.csv"),
        electricity=temp("build/annual_heat_electricity_consumption.csv"),
    script: "../src/construct/annual_heat_demand.py"


rule annual_subnational_demand:
    message: "Scale national demand to {wildcards.resolution} resolution"
    input:
        population = landeligibility("build/{resolution}/population.csv"),
        units = landeligibility("build/{resolution}/units.geojson"),
        annual_demand = rules.annual_heat_demand.output.demand,
        annual_heat_electricity_consumption = rules.annual_heat_demand.output.electricity,
        industry_demand = rules.annual_industry_demand.output[0],
        road_distance = rules.annual_transport_demand.output.distance,
        rail_demand = rules.annual_transport_demand.output.rail_energy,
        air_demand = rules.annual_transport_demand.output.air_energy,
        marine_demand = rules.annual_transport_demand.output.marine_energy,
        road_bau_electricity=rules.annual_transport_demand.output.road_bau_electricity,
        rail_bau_electricity=rules.annual_transport_demand.output.rail_bau_electricity,
        industry_bau_electricity=rules.annual_industry_demand.output.bau_electricity,
        emissions = "data/industry/all_industrial_ets_eprtr_sites.geojson",
        freight = rules.eurostat_data_tsv.output.freight,
        employees = rules.eurostat_data_tsv.output.employees,
        gva = rules.eurostat_data_tsv.output.gva,
        ch_gva = rules.ch_data_xlsx.output.gva,
        nuts_to_regions = "data/eurospores.xlsx",
        industry_activity_codes = "data/industry/industry_activity_codes.csv"
    conda: "../envs/geodata.yaml"
    params:
        scaling_factors = config["scaling-factors"],
    output:
        all_annual = "build/model/{resolution}/annual-demand.csv",
    script: "../src/construct/annual_subnational_demand.py"


rule when2heat:
    message: "Clone when2heat github repo"
    output: directory("data/automatic/when2heat")
    shell:
        "git clone {URL_WHEN2HEAT} {output}"


rule heat_demand_profiles:
    message: "Generate unscaled hourly profiles at {wildcards.resolution} resolution."
    input:
        population = landeligibility("build/population-europe.tif"),
        units = landeligibility("build/{resolution}/units.geojson"),
        air_temp = "data/weather/temperature.nc",
        wind_10m = "data/weather/wind10m.nc",
        when2heat = rules.when2heat.output[0],
    params:
        model_year = config["year"]
    conda: "../envs/geodata.yaml"
    output:
        space_heating=temp("build/model/{resolution}/space-heat-profile.csv"),
        water_heating=temp("build/model/{resolution}/water-heat-profile.csv")
    script: "../src/construct/hourly_heat_profiles.py"


rule scaled_heat_demand_profiles:
    message: "Scale hourly profiles at {wildcards.resolution} resolution according to annual demand."
    input:
        units = landeligibility("build/{resolution}/units.geojson"),
        annual_demand = rules.annual_subnational_demand.output.all_annual,
        dwellings = rules.eurostat_data_tsv.output.dwellings,
        nuts_to_regions = "data/eurospores.xlsx",
        space_profiles = rules.heat_demand_profiles.output.space_heating,
        water_profiles = rules.heat_demand_profiles.output.water_heating,
        cooking_profiles = 'data/cooking_profiles.csv.gz'
    params:
        model_year = config["year"],
        space_heating_key = 'space_heating',
        water_heating_key = 'water_heating',
        cooking_key = 'cooking',
    conda: "../envs/geodata.yaml"
    output:
        space_heating="build/model/{resolution}/space-heat-demand.csv",
        water_heating="build/model/{resolution}/water-heat-demand.csv",
        cooking="build/model/{resolution}/cooking-demand.csv",
    script: "../src/construct/scale_hourly_heat_profiles.py"


rule scaled_bau_electricity_heat_demand_profiles:
    message: "Scale hourly profiles at {wildcards.resolution} resolution according to annual demand."
    input:
        units = landeligibility("build/{resolution}/units.geojson"),
        annual_demand = rules.annual_subnational_demand.output.all_annual,
        dwellings = rules.eurostat_data_tsv.output.dwellings,
        nuts_to_regions = "data/eurospores.xlsx",
        space_profiles = rules.heat_demand_profiles.output.space_heating,
        water_profiles = rules.heat_demand_profiles.output.water_heating,
        cooking_profiles = 'data/cooking_profiles.csv.gz'
    params:
        model_year = config["year"],
        space_heating_key = 'space_heating_bau_electricity',
        water_heating_key = 'water_heating_bau_electricity',
        cooking_key = 'cooking_bau_electricity',
    conda: "../envs/geodata.yaml"
    output:
        space_heating=temp("build/model/{resolution}/space-heat-bau-electricity-demand.csv"),
        water_heating=temp("build/model/{resolution}/water-heat-bau-electricity-demand.csv"),
        cooking=temp("build/model/{resolution}/cooking-bau-electricity-demand.csv")
    script: "../src/construct/scale_hourly_heat_profiles.py"


rule scaled_public_transport_demand_profiles:
    message: "Scale hourly profiles at {wildcards.resolution} resolution according to annual demand."
    input:
        annual_demand = rules.annual_subnational_demand.output.all_annual,
        rail_profiles = "data/transport/rail_daily_profiles_destinee.csv",
    params:
        model_year = config["year"]
    conda: "../envs/default.yaml"
    output: temp("build/model/{resolution}/public-transport-demand.csv")
    script: "../src/construct/scale_hourly_transport_profiles.py"


rule update_electricity_with_other_sectors:
    message: "Creating new electricity {wildcards.resolution} timeseries without heat demand"
    input:
        space_heating = rules.scaled_bau_electricity_heat_demand_profiles.output.space_heating,
        water_heating = rules.scaled_bau_electricity_heat_demand_profiles.output.water_heating,
        cooking = rules.scaled_bau_electricity_heat_demand_profiles.output.cooking,
        public_transport_profiles = rules.scaled_public_transport_demand_profiles.output[0],
        annual_demand = rules.annual_subnational_demand.output.all_annual,
        hourly_electricity = eurocalliope("build/model/{resolution}/electricity-demand.csv")
    params:
        model_year = config["year"]
    conda: "../envs/geodata.yaml"
    output: "build/model/{resolution}/electricity-demand.csv"
    script: "../src/construct/electricity_with_other_sectors.py"


rule heat_pump_characteristics:
    message: "Generate timeseries of heat pump coefficients of performance (COPs)"
    input:
        units = landeligibility("build/{resolution}/units.geojson"),
        path_to_population_tif = landeligibility("build/population-europe.tif"),
        soil_temp = "data/weather/tsoil5.nc",
        air_temp = "data/weather/temperature.nc",
        hp_characteristics = "data/heat_pump_characteristics.csv",
    params:
        heat_tech_params = config["parameters"]["heat-end-use"],
        model_year = config["year"]
        characteristic = "{characteristic}"
    conda: "../envs/geodata.yaml"
    output:
        gshp_water_heat = "build/model/{resolution}/{characteristic}-gshp-water-heat.csv",
        ashp_water_heat = "build/model/{resolution}/{characteristic}-ashp-water-heat.csv",
        gshp_space_heat = "build/model/{resolution}/{characteristic}-gshp-space-heat.csv",
        ashp_space_heat = "build/model/{resolution}/{characteristic}-ashp-space-heat.csv"
    script: "../src/construct/heat_pump_characteristics.py"


rule ev_energy_cap:
    message: "Restructing RAMP-mobility EV plug-in profiles for use in Calliope"
    input:
        units = landeligibility("build/{resolution}/units.geojson"),
        ev_profiles = "data/transport/ev_profiles_ramp.csv",
    params:
        model_year
    conda: "../envs/geodata.yaml"
    output:
        ev_energy_cap = "build/{resolution}/energy-cap-ev.csv"
    script: "../src/construct/hourly_ev_profiles.py"


rule model:
    message: "Build entire model on resolution {wildcards.resolution}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/directional-rooftop.yaml",
        rules.links.output,
        rules.outer_countries.output,
        "build/{resolution}/energy-cap-ev.csv",
        expand(
            "build/model/{{resolution}}/{characteristic}-{technology}-{end-use}.csv",
            characteristic=["energy-cap", "cop"], technology=["ahp", "gshp"],
            end_use=["space-heat", "water-heat"]
        )
        expand(
            "build/model/{{resolution}}/{end_use}-demand.csv",
            end_use=["electricity", "space-heat", "water-heat", "cooking", "annual"]
        )
        "build/model/{resolution}/cop.csv",
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=[
                "open-field-pv", "rooftop-pv-n", "rooftop-pv-e-w",
                "rooftop-pv-s-flat", "wind-offshore", "wind-onshore",
                "hydro-ror", "hydro-reservoir-inflow", "rooftop-pv"
            ],
        ),
        definition = "src/template/model.yaml"
    output:
        model = "build/model/{resolution}/model.yaml"
    shell:
        "cp {input.definition} {output}"
