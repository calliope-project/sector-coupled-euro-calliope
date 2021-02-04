import math
import pandas as pd

import util


LEVEL_ORDER = ['country_code', 'year', 'unit', 'end_use']
YEAR_RANGE = range(2000, 2019)  # FIXME: this should be in config

idx = pd.IndexSlice


def combine_data(
    building_demand, building_heat_electricity_consumption,
    industry_demand, road_distance, road_vehicles, rail_demand, air_demand, marine_demand,
    road_bau_electricity, rail_bau_electricity, industry_bau_electricity,
    out_path, scaling_factors, industry_config, countries
):
    road_distance_df = util.read_tdf(road_distance)
    road_vehicles_df = util.read_tdf(road_vehicles)
    rail_demand_df = util.read_tdf(rail_demand)
    road_bau_electricity_df = util.read_tdf(road_bau_electricity)
    rail_bau_electricity_df = util.read_tdf(rail_bau_electricity)
    building_demand_df = util.read_tdf(building_demand)
    building_heat_electricity_consumption_df = util.read_tdf(building_heat_electricity_consumption)
    industry_demand_df = util.read_tdf(industry_demand)
    industry_bau_electricity_df = util.read_tdf(industry_bau_electricity)

    df = add_non_industry_demand(
        road_distance_df, road_vehicles_df, rail_demand_df, road_bau_electricity_df,
        rail_bau_electricity_df, building_demand_df, building_heat_electricity_consumption_df
    )
    df = pd.concat([
        df,
        add_industry_demand(
            industry_demand_df, road_distance_df, road_vehicles_df, rail_demand_df,
            air_demand, marine_demand, rail_bau_electricity_df, industry_bau_electricity_df,
            industry_config
        )
    ])

    # FIXME: move this yearrange into config yaml
    df = df[df.index.get_level_values('year').isin(YEAR_RANGE)].rename_axis(index={'country_code': 'id'})

    # Select only countries of interest for the model
    df = df[df.index.get_level_values('id').isin([util.get_alpha3(i) for i in countries])]
    assert len(df.index.get_level_values('id').unique()) == len(countries)

    # FillNA for any missing info
    _df = df.unstack('id').stack(dropna=False)
    print("Filling missing data for {}".format(_df[_df.isna()].sum(level=['dataset', 'cat_name', 'id', 'year']).index.remove_unused_levels()))
    df = _df.fillna(0).groupby(level=['dataset', 'cat_name', 'year', 'id', 'unit', 'end_use']).sum()

    # Scale
    df.loc[df.index.get_level_values('unit') == 'twh'] *= scaling_factors['energy']
    df.loc[df.index.get_level_values('unit') == 'kt'] *= scaling_factors['co2']
    df.loc[df.index.get_level_values('unit') == 'mio_km'] *= scaling_factors['transport']
    df = df.rename(
        {'twh': '{:.2f} twh'.format(1 / scaling_factors['energy']),
         'kt': '{:.2f} kt'.format(1 / scaling_factors['co2']),
         'mio_km': '{:.2f} mio_km'.format(1 / scaling_factors['transport'])},
        level='unit'
    )
    df.to_csv(out_path)


def add_dims(df, keys):
    return pd.concat(
        [df.rename(util.get_alpha3, level='country_code').reorder_levels(LEVEL_ORDER)],
        names=['dataset', 'cat_name'], keys=keys
    )


def get_transport_data(df, vehicle_type, vehicle_abbreviation):
    return df.xs(vehicle_type).to_frame(vehicle_abbreviation).rename_axis(columns='end_use').stack()


def add_non_industry_demand(
    road_distance_df, road_vehicles_df, rail_demand_df, road_bau_electricity_df,
    rail_bau_electricity_df, building_demand_df, building_heat_electricity_consumption_df
):
    df_list = []
    # Household heat, and total current heat met by electricity
    df_list.append(add_dims(
        building_demand_df.xs('household', level='cat_name'),
        [('heat_demand', 'household')]
    ))
    df_list.append(add_dims(
        building_heat_electricity_consumption_df.xs('household', level='cat_name'),
        [('heat_demand', 'household')]
    ))
    # Commercial heat, and total current heat met by electricity
    df_list.append(add_dims(
        building_demand_df.xs('commercial', level='cat_name'),
        [('heat_demand', 'commercial')]
    ))
    df_list.append(add_dims(
        building_heat_electricity_consumption_df.xs('commercial', level='cat_name'),
        [('heat_demand', 'commercial')]
    ))

    def _passenger_vehicle_dfs(df):
        new_names = {
            'Motor coaches, buses and trolley buses': 'bus',
            'Passenger cars': 'passenger_car', 'Powered 2-wheelers': 'motorcycle'
        }
        passenger_df = (
            df.drop(['Heavy duty vehicles', 'Light duty vehicles'], level=0)
        )
        passenger_df.index = passenger_df.index.set_names('end_use', level='vehicle_type')
        return passenger_df.rename(new_names, level='end_use')

    # passenger car distance and # vehicles, and energy thereof met by electricity
    df_list.append(add_dims(
        _passenger_vehicle_dfs(road_distance_df), [('transport_demand', 'road')]
    ))
    df_list.append(add_dims(
        _passenger_vehicle_dfs(road_bau_electricity_df).rename(lambda x: f'{x}_bau_electricity', level='end_use'),
        [('transport_demand', 'road')]
    ))
    df_list.append(add_dims(
        _passenger_vehicle_dfs(road_vehicles_df), [('transport_vehicles', 'road')]
    ))

    # passenger rail distance, and energy thereof met by electricity
    df_list.append(add_dims(
        rail_demand_df
        .drop('Freight', level=0)
        .sum(level=['country_code', 'year', 'unit'])
        .to_frame('electricity')
        .rename_axis(columns='end_use')
        .stack(),
        [('transport_demand', 'rail')]
    ))
    df_list.append(add_dims(
        rail_bau_electricity_df
        .xs('Passenger transport')
        .sum(level=['country_code', 'year', 'unit'])
        .to_frame('bau_electricity')
        .rename_axis(columns='end_use')
        .stack(),
        [('transport_demand', 'rail')]
    ))

    # commercial, light duty car distance and # vehicles, and energy thereof met by electricity
    df_list.append(add_dims(
        get_transport_data(road_distance_df, 'Light duty vehicles', 'ldv'),
        [('transport_demand', 'road')]
    ))
    df_list.append(add_dims(
        get_transport_data(road_vehicles_df, 'Light duty vehicles', 'ldv'),
        [('transport_vehicles', 'road')]
    ))
    df_list.append(add_dims(
        get_transport_data(road_bau_electricity_df, 'Light duty vehicles', 'ldv_bau_electricity'),
        [('transport_demand', 'road')]
    ))

    return pd.concat(df_list)


def add_industry_demand(
    industry_demand_df, road_distance_df, road_vehicles_df, rail_demand_df,
    air_demand, marine_demand, rail_bau_electricity_df, industry_bau_electricity_df,
    industry_config
):
    """
    Some industry subsectors are subdivided based on employment data per NUTS3 region.
    Where this isn't available, we use freight loading data. Both datasets have
    subcategories relating to industry subsectors.
    The more emitting subsectors (e.g. iron & steel) are subdivided using EU-ETS emissions
    data.
    """
    df_list = []
    df_list.append(add_dims(
        industry_demand_df
        .sum(level=['country_code', 'unit', 'year', 'carrier'])
        .rename_axis(index={'carrier': 'end_use'}),
        [('industry_demand', 'industry')]
    ))
    df_list.append(add_dims(
        industry_bau_electricity_df
        .sum(level=['country_code', 'unit', 'year'])
        .to_frame('bau_electricity')
        .rename_axis(columns='end_use')
        .stack(),
        [('industry_demand', 'industry')]
    ))

    df_list.append(add_dims(
        get_transport_data(road_distance_df, 'Heavy duty vehicles', 'hdv'),
        [('industry_demand', 'road')]
    ))
    df_list.append(add_dims(
        get_transport_data(road_vehicles_df, 'Heavy duty vehicles', 'hdv'),
        [('transport_vehicles', 'road')]
    ))

    df_list.append(add_dims(
        get_transport_data(rail_demand_df, 'Freight', 'electricity'),
        [('industry_demand', 'rail')]
    ))

    # BAU electricity
    df_list.append(add_dims(
        get_transport_data(rail_bau_electricity_df, 'Freight transport', 'bau_electricity'),
        [('industry_demand', 'rail')]
    ))

    # Air and Marine energy demand is also distributed according to total industry
    df_list.append(add_dims(
        util.read_tdf(air_demand)
        .rename_axis(index={'country': 'country_code'})
        .to_frame('kerosene')
        .rename_axis(columns='end_use')
        .stack(),
        [('industry_demand', 'air')]
    ))
    df_list.append(add_dims(
        util.read_tdf(marine_demand)
        .rename_axis(index={'country': 'country_code'})
        .to_frame('diesel')
        .rename_axis(columns='end_use')
        .stack(),
        [('industry_demand', 'marine')]
    ))
    return electrify_industry_demand(pd.concat(df_list), industry_config)


def electrify_industry_demand(df, industry_config, level_order=LEVEL_ORDER):

    # Electrify any end uses according to workflow configuration
    # electrification 'efficiencies' are given as 'output/twh_electricity_input'
    to_electrify = [
        i
        for i in industry_config["electrification_efficiency"].keys()
        if i not in industry_config["carriers"]
    ]
    df_to_electrify = df.unstack('end_use').loc[:, to_electrify]
    idx_to_drop = df_to_electrify.stack().index.remove_unused_levels()
    electrification = {
        i: industry_config["electrification_efficiency"][i] if i in to_electrify else 1
        for i in df_to_electrify.columns
    }

    df_to_electrify = (
        df_to_electrify
        .div(electrification)
        .sum(axis=1, min_count=1)
        .where(lambda x: x > 0)
        .dropna()
        .to_frame(('twh', 'electricity'))
        .rename_axis(columns=['unit', 'end_use'])
        .droplevel('unit', axis=0)
        .stack([0, 1])
        .reorder_levels(['dataset', 'cat_name'] + level_order)
    )

    new_df = pd.concat([df.drop(idx_to_drop), df_to_electrify]).groupby(level=df.index.names).sum()
    # Sanity checks
    # only the expected end uses have been removed once electrified
    assert (
        set(df.index.get_level_values('end_use').unique())
        .difference(new_df.index.get_level_values('end_use').unique())
        == set(to_electrify)
    )
    # new electricity demand has increased, (or stayed the same if len(to_electrify) == 0)
    assert new_df.xs('electricity', level='end_use').sum() >= df.xs('electricity', level='end_use').sum()
    # all other demand has remained static

    assert math.isclose(
        new_df.drop('electricity', level='end_use').sum(),
        df.drop(['electricity'] + to_electrify, level='end_use').sum()
    )
    return new_df


def load_eurostat_tsv(path_to_tsv, index_names, slice_idx=None, slice_lvl=None):
    df = pd.read_csv(path_to_tsv, delimiter='\t', index_col=0)
    df.index = df.index.str.split(',', expand=True).rename(index_names)
    if slice_idx is not None:
        df = df.xs(slice_idx, level=slice_lvl)
    df.columns = df.columns.astype(int)
    return df.apply(util.to_numeric)


if __name__ == "__main__":
    combine_data(
        building_demand=snakemake.input.annual_heat_demand,
        building_heat_electricity_consumption=snakemake.input.annual_heat_electricity_consumption,
        industry_demand=snakemake.input.industry_demand,
        road_distance=snakemake.input.road_distance,
        road_vehicles=snakemake.input.road_vehicles,
        rail_demand=snakemake.input.rail_demand,
        air_demand=snakemake.input.air_demand,
        marine_demand=snakemake.input.marine_demand,
        road_bau_electricity=snakemake.input.road_bau_electricity,
        rail_bau_electricity=snakemake.input.rail_bau_electricity,
        industry_bau_electricity=snakemake.input.industry_bau_electricity,
        scaling_factors=snakemake.params.scaling_factors,
        industry_config=snakemake.params.industry_config,
        countries=snakemake.params.countries,
        out_path=snakemake.output.all_annual
    )
