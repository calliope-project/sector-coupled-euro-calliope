# Industry datasets
This directory provides a number of datasets relevant to building end use energy demand data for industry. Many of the datasets are ultimately not used in Euro-Calliope, but are available for comparison, since no dataset is representative of the 'ground truth'. This README gives more detail on data sources and the relevant Jupyter notebooks in which processing / data comparison is undertaken.

# Datasets
1. FORECAST-industry_2012_end_use_consumption.csv
   - countries: EU28
   - subsectors: Eurostat compatible
   - energy carriers: Partially aggregated relative to Eurostat
   - years: 2012
   - end uses: space heating, process heating (inc. temperature levels), cooling. No end use electricity consumption

    This dataset is a preprocessed version of the data published by [Rehfeldt, 2018] and available in its original Excel format as an output of the project "Mapping and analyses of the current and future (2020 - 2030) heating/cooling fuel deployment (fossil/renewables)" [here](https://ec.europa.eu/energy/studies/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment_en). This was preprocessed by hand since the Excel has deleted data that was accessed from a pivot table that had luckily not been refreshed.

2. FORECAST-industry_2015_end_use_consumption.csv
   - countries: EU28
   - subsectors: None
   - energy carriers: Partially aggregated relative to Eurostat
   - years: 2015
   - end uses: space heating, process heating (inc. temperature levels), cooling. No end use electricity consumption

    This dataset is available in raw form as an output of the HeatRoadmap Europe project [here](https://heatroadmap.eu/heating-and-cooling-energy-demand-profiles/). It is structured similarly to (1), but without the need to preprocess by hand. Data is instead processed in "process_FORECAST_industry_2015.ipynb". Note that unlike the 2012 data, there is no disaggregation of data into industry subsectors.

1. jrc_idees_processed_energy.csv.gz
   - countries: EU28
   - subsectors: Eurostat compatible
   - energy carriers: Partially aggregated relative to Eurostat
   - years: 2000-2015
   - end uses: All, including many processes

    This dataset is only available in its raw form following application to JRC for access to the database [here](https://ec.europa.eu/jrc/en/potencia/jrc-idees). Once downloaded it includes detailed information on end use energy consumption in the commercial, household, and industry sectors. The data on end uses in each industry subsector is highly detailed, including individual process (e.g. grinding, annealing, rolling). 'space heating' is not directly referred to, but there is a category for 'low enthalpy heat'. Processing of the data is undertaken in "process_jrc_idees.ipynb", and requires some manual allocation of processes to end use energy demand. The manual allocation is given in "jrc_idees_industry_process_end_uses.csv". To extract data from the industry spreadsheets, the colour and indent of the text is extracted using the [StyleFrame](https://styleframe.readthedocs.io/en/latest/) package, with a selection of text colours and indentations referring to particular end use categories and carriers.

2. published_industry_end_use_consumption.csv
   - countries: **AUT**, **GRB**, **DEU**, **CHE**
   - subsectors: **AUT**: Eurostat compatible; **GRB**: partially aggregated relative to Eurostat; **DEU**: quite different to Eurostat; **CHE**: None
   - energy carriers: **AUT**: partially aggregated relative to Eurostat; **GRB**: aggregated to main carriers (e.g. 'Gas', 'Electricity'); **DEU**: partially aggregated relative to Eurostat; **CHE**: only 'fuel' and 'electricity'
   - years: **AUT**: 2005-2017; **GRB**: 2016-2018; **DEU**: 2013-2018; **CHE**: 2000-2018
   - end uses: **AUT**: 12, inc. 'space heating', 'water heating' and process heating processes; **GRB**: 'space heating' and different temperature level 'process heating'; **DEU**: 'space heating', 'process heating', 'water heating', and various mechanical/electrical end uses; **CHE**: 'space heating', 'process heating', and various mechanical/electrical end uses

    Datasets published by national statistics offices are collated in this file, having been processed in "process_published_industry_end_use_datasets.ipynb". Links to datasets are given in the Jupyter notebook, except for Germany, for which the Excel spreadsheet was provided by Clemens Rohde at ISI Fraunhofer, to save having to scrape the data from a PDF. For further details on datasets, see the following documents for [AUT](https://www.statistik.at/web_en/statistics/EnergyEnvironmentInnovationMobility/energy_environment/energy/useful_energy_analysis/index.html), [GRB](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/820843/Energy_Consumption_in_the_UK__ECUK__MASTER_COPY.pdf), [DEU](https://ag-energiebilanzen.de/index.php?article_id=29&fileName=isi_industrie_ghd_18.pdf), and [CHE](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiIo8_foYbrAhUD_aQKHYtfCkcQFjABegQIAxAB&url=https%3A%2F%2Fpubdb.bfe.admin.ch%2Fde%2Fpublication%2Fdownload%2F9853&usg=AOvVaw2qPDCnrxZcjv_I1r1zZmy-). There is no commonality in any of chosen energy carrier, temporal scope, end use categorisation, or industry subsectors. When comparing datasets, the data is aggregated to the lowest common denominator, or gaps in data are explicit.

# Comparisons between datasets
Published datasets from individual countries are compared in "process_published_industry_end_use_datasets.ipynb". Published datasets are compared to modelling efforts in "compare_industry_end_use_datasets.ipynb". Datasets which could be used to spatially disaggregate industry data are compared to Austrian subnational data in "compare_AT_subregional_consumption.ipynb".
