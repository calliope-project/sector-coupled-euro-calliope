import pandas as pd

import util

EFFICIENCY_QUANTILE = 0.25

CARRIERS = {
    'O4652XR5210B': 'petrol',
    'R5210B': 'biofuels',
    'O4671XR5220B': 'diesel',
    'R5220B': 'biofuels',
    'O4630': 'lpg',
    'G3000': 'natural_gas',
    'E7000': 'electricity'
}


def get_transport_demand(
    path_to_energy_balances, path_to_jrc_road_energy, path_to_jrc_road_distance,
    path_to_jrc_rail_energy, path_to_jrc_rail_distance, path_to_distance_output,
    path_to_efficiency_output, path_to_rail_energy_output, path_to_air_energy_output,
    path_to_marine_energy_output, path_to_road_bau_electricity, path_to_rail_bau_electricity
):
    energy_balances = util.read_tdf(path_to_energy_balances)
    road_energy_df = util.read_tdf(path_to_jrc_road_energy)
    road_distance_df = util.read_tdf(path_to_jrc_road_distance)
    rail_energy_df = util.read_tdf(path_to_jrc_rail_energy)
    rail_distance_df = util.read_tdf(path_to_jrc_rail_distance)

    total_road_distance, road_efficiency, road_bau_consumption = get_all_distance_efficiency(
        energy_balances, 'FC_TRA_ROAD_E', road_energy_df, road_distance_df, 'vehicle_subtype'
    )
    # Some cleanup that's specific to road data
    road_efficiency = (
        road_efficiency
        .unstack(0)
        .groupby({
            'Diesel oil engine': 'diesel',
            'Battery electric vehicles': 'electricity',
            'Domestic': 'diesel',
            'International': 'diesel',
            'Gasoline engine': 'petrol',
            'Plug-in hybrid electric': 'petrol'
        }).mean()
        .assign(unit='twh_per_mio_km').set_index('unit', append=True)
        .stack()
    )
    road_electricity_bau = (
        road_bau_consumption
        .sum(level=['carrier', 'vehicle_type', 'country_code', 'year', 'unit'])
        .xs('electricity')
    )

    total_rail_distance, rail_efficiency, rail_bau_consumption = get_all_distance_efficiency(
        energy_balances, 'FC_TRA_RAIL_E', rail_energy_df, rail_distance_df, 'carrier'
    )
    # Some cleanup that's specific to rail data
    rail_energy = (
        total_rail_distance.mul(
            rail_efficiency.xs('electricity', level='carrier'),
            axis=0, level='vehicle_type'
        )
        .sum(level=['vehicle_type', 'country_code'])
        .stack('year')
        .to_frame('twh')
        .rename_axis(columns='unit')
        .stack()
    )
    rail_electricity_bau = (
        rail_bau_consumption
        .sum(level=['carrier', 'section', 'country_code', 'year', 'unit'])
        .xs('electricity')
    )
    air_energy = (
        energy_balances
        .loc[['FC_TRA_DAVI_E', 'INTAVI']]
        .unstack('carrier_code')
        ['O4000XBIO']
        .apply(util.tj_to_twh)
        .sum(level=['country', 'year'])
        .to_frame('twh')
        .rename_axis(columns='unit')
        .stack()
    )
    marine_energy = (
        energy_balances
        .loc[['INTMARB', 'FC_TRA_DNAVI_E']]
        .unstack('carrier_code')
        ['O4000XBIO']
        .apply(util.tj_to_twh)
        .sum(level=['country', 'year'])
        .to_frame('twh')
        .rename_axis(columns='unit')
        .stack()
    )

    total_road_distance.stack('year').to_csv(path_to_distance_output)
    road_efficiency.to_csv(path_to_efficiency_output)
    rail_energy.to_csv(path_to_rail_energy_output)
    air_energy.to_csv(path_to_air_energy_output)
    marine_energy.to_csv(path_to_marine_energy_output)
    road_electricity_bau.to_csv(path_to_road_bau_electricity)
    rail_electricity_bau.to_csv(path_to_rail_bau_electricity)


def get_all_distance_efficiency(
    energy_balances, cat_name, energy_df, distance_df, unique_dim
):

    transport_energy_balance = (
        energy_balances
        .xs(cat_name)
        .unstack('carrier_code')
        .groupby(CARRIERS, axis=1).sum(min_count=1)
        .rename_axis(columns='carrier')
        .stack()
        .apply(util.tj_to_twh)
        .droplevel('unit')
        .rename_axis(index=['country_code', 'year', 'carrier'])
    )
    # contribution of each transport mode to carrier consumption from JRC_IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    carrier_contribution = fill_missing_countries_and_years(
        energy_df.div(energy_df.sum(level=['carrier', 'country_code', 'year']))
    )
    # Energy consumption per transport mode by mapping transport mode
    # carrier contributions to total carrier consumption
    transport_energy_per_mode = pd.merge(
        transport_energy_balance.to_frame('eurostat'),
        carrier_contribution.to_frame('jrc'),
        left_index=True,
        right_on=transport_energy_balance.index.names
    )
    transport_energy_per_mode = (
        transport_energy_per_mode['eurostat'].mul(transport_energy_per_mode['jrc'])
    )
    # Distance per unit energy consumed per transport mode according to JRC IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    transport_efficiency = fill_missing_countries_and_years(
        distance_df
        .where(distance_df > 0)
        .div(
            energy_df
            .where(energy_df > 0)
            .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        )
    )
    # Distance travelled per transport mode, including years 2016-2018,
    # based on JRC IDEES transport efficiency (2015 data for 2016-2018)
    transport_distance_all_years = (
        transport_energy_per_mode
        .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        .mul(transport_efficiency)
    )
    # Use the distance traveled based on Eurostat to fill in blanks in JRC data
    aligned_dfs = (
        transport_distance_all_years
        .reorder_levels(distance_df.index.names)
        .align(distance_df)
    )
    total_transport_distance = (
        aligned_dfs[1]
        .fillna(aligned_dfs[0])
        .sum(level=['vehicle_type', 'country_code', 'unit', 'year'])
        .unstack()
    )
    # 25th percentile of 2015 efficiency in TWh/mio km
    transport_efficiency = (
        (1 / transport_efficiency.xs(2015, level='year'))
        .unstack(['vehicle_type', unique_dim])
        .quantile(EFFICIENCY_QUANTILE)
    )

    return total_transport_distance, transport_efficiency, transport_energy_per_mode


def fill_missing_countries_and_years(jrc_data):
    jrc_data = jrc_data.unstack('country_code')
    balkan_countries = jrc_data[['BG', 'HR', 'HU', 'RO', 'EL']].mean(axis=1)
    nordic_countries = jrc_data[['SE', 'DK']].mean(axis=1)
    ch_neighbours = jrc_data[['DE', 'AT', 'FR', 'IT']].mean(axis=1)
    jrc_data = jrc_data.assign(
        AL=balkan_countries,
        BA=balkan_countries,
        ME=balkan_countries,
        MK=balkan_countries,
        RS=balkan_countries,
        NO=nordic_countries,
        IS=nordic_countries,
        CH=ch_neighbours,
    ).stack().unstack('year')
    jrc_data = jrc_data.assign(
        **{str(i): jrc_data[2015] for i in range(2016, 2019)}
    )
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()


if __name__ == "__main__":
    get_transport_demand(
        path_to_energy_balances=snakemake.input.energy_balances,
        path_to_jrc_road_energy=snakemake.input.jrc_road_energy,
        path_to_jrc_road_distance=snakemake.input.jrc_road_distance,
        path_to_jrc_rail_energy=snakemake.input.jrc_rail_energy,
        path_to_jrc_rail_distance=snakemake.input.jrc_rail_distance,
        path_to_distance_output=snakemake.output.distance,
        path_to_efficiency_output=snakemake.output.efficiency,
        path_to_rail_energy_output=snakemake.output.rail_energy,
        path_to_air_energy_output=snakemake.output.air_energy,
        path_to_marine_energy_output=snakemake.output.marine_energy,
        path_to_road_bau_electricity=snakemake.output.road_bau_electricity,
        path_to_rail_bau_electricity=snakemake.output.rail_bau_electricity,
    )
