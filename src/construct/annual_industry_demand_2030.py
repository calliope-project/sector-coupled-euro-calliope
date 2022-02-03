import pandas as pd
import numpy as np

import util

JRC_TO_EUROSTAT_CARRIER_MAPPING = {
    'Diesel oil (incl. biofuels)': "oil",
    'Electricity': "electricity",
    'Natural gas (incl. biogas)': "natural_gas",
    'Solar and geothermal': "renewable_heat"
}
EUROSTAT_CARRIER_SIMPLIFICATION = {
    'biofuel': "biofuel",
    'electricity': "electricity",
    'heat': "electricity",
    'manufactured_gas': "methane",
    'natural_gas': "methane",
    'oil': "diesel",
    # 'renewable_heat': None, We just discount this demand entirely
    'solid_fossil': "coal"
}

YEAR_RANGE = slice(2000, 2018)

final_projected_demand_level_order = ["subsector", "country_code", "unit", "carrier", "year"]
final_bau_demand_level_order = ["subsector", "country_code", "unit", "year"]


def get_industry_demand(
    path_to_energy_balances, path_to_cat_names, path_to_carrier_names,
    path_to_jrc_industry_end_use,
    path_to_new_output, path_to_bau_output
):
    energy_df = pd.read_csv(path_to_jrc_industry_end_use, index_col=[0, 1, 2, 3, 4, 5, 6])
    energy_df.columns = energy_df.columns.astype(int).rename('year')
    energy_balances = util.read_tdf(path_to_energy_balances)
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)

    space_heat_demand = (
        energy_df
        .xs("demand")
        .xs("Low enthalpy heat", level="subsection")
        .sum(level=['country_code', 'cat_name', 'unit'])
        .assign(carrier_name="space_heat")
        .set_index("carrier_name", append=True)
    )
    space_heat_consumption = (
        energy_df
        .xs("consumption")
        .xs("Low enthalpy heat", level="subsection")
        .sum(level=['country_code', 'cat_name', 'carrier_name', 'unit'])
        .rename(JRC_TO_EUROSTAT_CARRIER_MAPPING, level="carrier_name")
    )

    industry_energy_balances = (
        energy_balances
        .unstack('year')
        .loc[:, YEAR_RANGE]
        .unstack(['country', "unit"])
        .groupby([
            cat_names.jrc_idees.drop("FC_IND_NE").dropna().to_dict(),
            carrier_names.ind_carrier_name.dropna().to_dict()
        ], level=['cat_code', 'carrier_code']).sum(min_count=1)
        .stack(['country', "unit"])
        .rename_axis(index=["cat_name", "carrier_name", "country_code", "unit"])
        .apply(util.tj_to_ktoe)
        .rename({"TJ": "ktoe"}, level="unit")
    )
    electricity_bau = (
        industry_energy_balances
        .xs("electricity", level="carrier_name")
        .bfill(axis=1)  # missing old data (BIH and MNE) backfilled
        .stack()
        .apply(util.ktoe_to_twh)
        .rename({'ktoe': 'twh'}, level='unit')
        .rename_axis(index={"cat_name": "subsector"})
    )
    industry_energy_balances = (
        industry_energy_balances
        .rename(EUROSTAT_CARRIER_SIMPLIFICATION)
        .sum(level=industry_energy_balances.index.names)
    )
    filled_space_heat_consumption = fill_missing_data(
        industry_energy_balances, space_heat_consumption
    )
    filled_space_heat_demand = fill_missing_data(
        industry_energy_balances, space_heat_demand
    )
    industry_energy_balances_without_heat = (
        industry_energy_balances.subtract(
            filled_space_heat_consumption
            .reorder_levels(industry_energy_balances.index.names)
            .reindex(industry_energy_balances.index)
            .fillna(0)
        )
        .clip(lower=0)
    )
    industry_energy_balances = (
        industry_energy_balances_without_heat.append(
            filled_space_heat_demand.reorder_levels(industry_energy_balances.index.names)
        )
        .apply(util.ktoe_to_twh)
        .rename({'ktoe': 'twh'}, level='unit')
        .rename_axis(index={"cat_name": "subsector", "carrier_name": "carrier"})
        .reorder_levels(['subsector', 'country_code', 'unit', 'carrier'])
        .stack()
    )

    industry_energy_balances.reorder_levels(final_projected_demand_level_order).to_csv(path_to_new_output)
    electricity_bau.reorder_levels(final_bau_demand_level_order).to_csv(path_to_bau_output)


def fill_missing_data(energy_balances, energy_consumption):
    """
    There are 7 countries without relevant data in JRC_IDEES, so we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    2016-2018 data for all 35 countries is filled in based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    # Get annual energy balances
    energy_balances = energy_balances.sum(level=['cat_name', "country_code"], min_count=1)
    JRC_countries = energy_consumption.index.get_level_values("country_code")
    eurostat_countries = energy_balances.index.get_level_values("country_code")
    # If JRC-IDEES data exists, there will be data in 'energy_consumption' for that country
    balances_with_jrc_data = energy_balances.loc[
        eurostat_countries.isin(JRC_countries), energy_consumption.columns
    ]
    # If JRC-IDEES data does not exist, there will only be data for that country in annual energy balances
    balances_without_jrc_data = energy_balances.drop(
        JRC_countries.unique(), level='country_code'
    )
    # Compared to countries with JRC data, those without JRC-IDEES data consume X amount of energy
    # E.g. CH has ~1% of paper and pulp energy consumption compared to all countries with JRC data
    countries_without_jrc_data_contribution = balances_without_jrc_data.div(
        balances_with_jrc_data.sum(level='cat_name')
    )
    # Multiply the JRC data with the contribution from each country
    # So all JRC countries consume 3.43e4ktoe electricity in paper and pulp in 2014,
    # hence CH consumes ~4.34e2ktoe electricity in paper and pulp in 2014
    countries_without_jrc_data_consumption = (
        energy_consumption
        .sum(level=['cat_name', 'unit', 'carrier_name'])
        .align(countries_without_jrc_data_contribution)[0]
        .mul(energy_consumption
             .sum(level=['cat_name', 'unit', 'carrier_name'])
             .align(countries_without_jrc_data_contribution)[1])
    )
    all_euro_calliope_consumption = energy_consumption.append(
        countries_without_jrc_data_consumption
        .reorder_levels(energy_consumption.index.names)
    )
    average_consumption_per_energy_use = (
        all_euro_calliope_consumption
        .div(energy_balances)
        .mean(axis=1)
        .where(lambda x: ~np.isinf(x) & (x > 0))
    )
    # In some instances (e.g. RO: non ferrous metals), JRC IDEES has consumption, but
    # latest Eurostat doesn't, leading to inf values
    average_consumption_per_energy_use = (
        average_consumption_per_energy_use
        .unstack(['cat_name', 'carrier_name', 'unit'])
        .fillna(
            average_consumption_per_energy_use
            .mean(level=['cat_name', 'carrier_name', 'unit'])
        )
        .unstack()
        .reorder_levels(all_euro_calliope_consumption.index.names)
    )

    filled_2016_to_2018 = all_euro_calliope_consumption.where(lambda x: x > 0).fillna(
        energy_balances.mul(average_consumption_per_energy_use, axis=0)
    )
    fill_mid_gaps = filled_2016_to_2018.interpolate(method='linear', limit_area='inside', axis=1)
    all_filled = fill_mid_gaps.T.fillna(fill_mid_gaps.mean(axis=1)).T

    return all_filled


if __name__ == "__main__":
    get_industry_demand(
        path_to_energy_balances=snakemake.input.energy_balances,
        path_to_cat_names=snakemake.input.cat_names,
        path_to_carrier_names=snakemake.input.carrier_names,
        path_to_jrc_industry_end_use=snakemake.input.jrc_industry_end_use,
        path_to_new_output=snakemake.output.new_demand,
        path_to_bau_output=snakemake.output.bau_electricity
    )
