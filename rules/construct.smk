
URL_ENERGY_BALANCE = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/nrg_bal_c.tsv.gz"
URL_HH_END_USE = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/nrg_bal_c.tsv.gz"
URL_CH_HH_END_USE = "https://www.bfe.admin.ch/bfe/en/home/versorgung/statistik-und-geodaten/energiestatistiken/energieverbrauch-nach-verwendungszweck.exturl.html/aHR0cHM6Ly9wdWJkYi5iZmUuYWRtaW4uY2gvZGUvcHVibGljYX/Rpb24vZG93bmxvYWQvOTg1NA==.html"
URL_CH_ENERGY_BALANCE = "https://www.bfe.admin.ch/bfe/en/home/versorgung/statistik-und-geodaten/energiestatistiken/gesamtenergiestatistik.exturl.html/aHR0cHM6Ly9wdWJkYi5iZmUuYWRtaW4uY2gvZGUvcHVibGljYX/Rpb24vZG93bmxvYWQvNzUxOQ==.html"
URL_DWELLINGS = "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/cens_11dwob_r3.tsv.gz"

subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

subworkflow landeligibility:
    workdir: "./land-eligibility"
    snakefile: "./land-eligibility/Snakefile"
    configfile: "./land-eligibility/config/default.yaml"

localrules: copy_euro_calliope, copy_resolution_specific_euro_calliope, model, links, outer_countries, heat_demand, eurostat_data_tsv, ch_data_xlsx, dwelling_type_tsv
ruleorder: model > links > outer_countries > heat_demand > copy_euro_calliope > copy_resolution_specific_euro_calliope
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
    message: "Get annual energy balances and household end uses from Eurostat"
    output:
        energy_balance = protected("data/automatic/annual_energy_balances.tsv.gz"),
        hh_end_use = protected("data/automatic/hh_end_use.tsv.gz")
    shell:
        """
        curl -sLo {output.energy_balance} '{URL_ENERGY_BALANCE}'
        curl -sLo {output.hh_end_use} '{URL_HH_END_USE}'
        """


rule ch_data_xlsx:
    message: "Get Swiss annual energy balances and household end uses"
    output:
        energy_balance = protected("data/automatic/ch_annual_energy_balances.xlsx"),
        hh_end_use = protected("data/automatic/ch_hh_end_use.xlsx")
    shell:
        """
        curl -sLo {output.energy_balance} '{URL_CH_ENERGY_BALANCE}'
        curl -sLo {output.hh_end_use} '{URL_CH_HH_END_USE}'
        """


rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        energy_balance = rules.eurostat_data_tsv.output.energy_balance,
        ch_energy_balance = rules.ch_data_xlsx.output.energy_balance,
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output:
        "build/annual_energy_balances.nc"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    shadow: "minimal"
    script: "../src/construct/annual_energy_balance.py"


rule annual_heat_consumption:
    message: "Compile total annual consumption for different end uses and energy carriers at national level"
    input:
        hh_end_use = rules.eurostat_data_tsv.output.hh_end_use,
        ch_end_use = rules.ch_data_xlsx.output.hh_end_use
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    output:
        "build/annual_heat_consumption.csv"
    script: "../src/construct/annual_heat_consumption.py"


rule pop_hdd_merra2:
    message: "Map population and heating degree days according to MERRA-2 gridcells"
    input:
        population = landeligibility("build/population-europe.tif"),
        air_temp = landeligibility("data/capacityfactors/temperature.nc"),
        units = landeligibility("build/national/units.geojson"),
    params:
        temperature_threshold = config["parameters"]["temperature-threshold"]
    conda: "../envs/geodata.yaml"
    output:
        "build/pop_hdd_merra2.csv"
    script: "../src/construct/map_to_merra2.py"


rule dwelling_type_tsv:
    message: "Get data on single-family vs multi-family dwellings"
    output:
        protected("data/automatic/dwellings.tsv.gz"),
    shell:
        """
        curl -sLo {output} '{URL_DWELLINGS}'
        """


rule model:
    message: "Build entire model on resolution {wildcards.resolution}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/electricity-demand.csv",
        "build/model/{resolution}/directional-rooftop.yaml",
        rules.links.output,
        rules.outer_countries.output,
        "build/model/{resolution}/heat-demand.csv",
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=[
                "open-field-pv", "rooftop-pv-n", "rooftop-pv-e-w",
                "rooftop-pv-s-flat", "wind-offshore", "wind-onshore",
                "hydro-ror", "hydro-reservoir-inflow", "rooftop-pv"
            ],
        ),
        definition = "../src/template/model.yaml"
    output:
        model = "build/model/{resolution}/model.yaml"
    shell:
        "cp {input.definition} {output}"
