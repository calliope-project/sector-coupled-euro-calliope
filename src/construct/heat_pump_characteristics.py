import pandas as pd
from scipy.interpolate import interp1d
import numpy as np

import util

idx = pd.IndexSlice

SOURCE = {
    'gshp': 'soil',
    'ashp': 'air',
    'hp': 'combined'
}


def get_characteristic(
    weather_pop, hp_characteristics, annual_demand,
    heat_tech_params, characteristic, tech, sink, model_year, out_path
):
    """
    """
    weather_pop_df = (
        pd.read_csv(weather_pop, index_col=0, parse_dates=True, header=[0, 1, 2, 3, 4])
        .filter(regex='temp')
        .rename_axis(columns={'dataset': 'source'})
        .rename(columns={'soil_temp': 'soil', 'air_temp': 'air'})
    )
    hp_data_df = util.read_tdf(hp_characteristics)

    if characteristic == "cop":
        hp_data_df = (
            hp_data_df.xs('COP', level='data_type')
            .groupby(level=['source', 'source_temp', 'sink_temp'])  # average of all products
            .mean()
            .dropna()
        )
        correction = 0.85  # correction factor of 0.85 to go to 'real' system - see when2heat

    elif characteristic == "energy-cap":
        hp_data_df = (
            hp_data_df.xs('Heating capacity', level='data_type')
            .unstack('source')
            .apply(cap_to_cap_ratio, axis=0)
            .stack()
            .groupby(level=['source', 'source_temp', 'sink_temp'])  # average of all products
            .mean()
            .dropna()
        )
        correction = 1

    sink_t_dict = heat_tech_params['heating_temperature']

    site_data = (
        pd.concat(
            [weather_pop_df
             .apply(temperature_to_characteristic, hp_data=hp_data_df, sink_t=v)
             .dropna(axis=1, how='all')
             for v in sink_t_dict.values()],
            keys=list(sink_t_dict.keys()), names=['sink'], axis=1
        )
        .mul(correction)
    )

    def _pop_wavg(x):
        weighted_avg = np.average(x, weights=x.columns.get_level_values('population').astype(float).fillna(0), axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    region_data = site_data.groupby(level=['id', 'sink', 'source'], axis=1).apply(_pop_wavg)

    if "heat_pump_ratio" in heat_tech_params.keys():
        # we create a dataset representing a combination of GSHP and ASHP
        _ratios = pd.Series(heat_tech_params["heat_pump_ratio"]).rename(SOURCE)
        hp_data = region_data.mul(_ratios, axis=1, level='source').sum(level=['id', 'sink'], axis=1)
        hp_data.columns = pd.MultiIndex.from_frame(hp_data.columns.to_frame().assign(source='combined'))
        region_data = region_data.merge(hp_data, left_index=True, right_index=True)

    if sink == 'water-heat':
        region_data.xs(('water', SOURCE[tech]), level=('sink', 'source'), axis=1).to_csv(out_path)
    elif sink == 'space-heat':
        space_heat_data(region_data, heat_tech_params['heat_sink_ratio'], SOURCE[tech]).to_csv(out_path)
    elif sink == 'heat':
        water = region_data.xs(('water', SOURCE[tech]), level=('sink', 'source'), axis=1)
        space = space_heat_data(region_data, heat_tech_params['heat_sink_ratio'], SOURCE[tech])
        annual_demand_df = util.read_tdf(annual_demand).xs(model_year, level="year")
        annual_water = annual_demand_df.xs(('water_heat'), level=('end_use')).sum(level='id')
        annual_space = annual_demand_df.xs(('space_heat'), level=('end_use')).sum(level='id')
        heat = (
            water.mul(annual_water, axis=1)
            .add(space.mul(annual_space, axis=1))
            .div(annual_water.add(annual_space))
        )
        heat.to_csv(out_path)


def space_heat_data(region_data, heat_sink_ratio_dict, source):
    def _heat_sink_wavg(x, ratio_keys, ratio_vals):
        weighted_avg = np.average(x.loc[:, idx[:, ratio_keys]], weights=ratio_vals, axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    return (
        region_data
        .xs(source, level='source', axis=1)
        .groupby('id', axis=1)
        .apply(_heat_sink_wavg,
               ratio_keys=[i for i in heat_sink_ratio_dict.keys()],
               ratio_vals=[i for i in heat_sink_ratio_dict.values()])
    )


def temperature_to_characteristic(x, hp_data, sink_t):
    """
    Take source (timeseries) and sink (static) temperatures and map them to COP
    or energy capacity based on heat pump performance data. temepratures outside points
    given by the performance data are linearly extrapolated/interpolated.
    """
    if x.name[0] == 'air':
        _data = hp_data.xs('air', level='source').unstack()
        source_t = x

    elif x.name[0] == 'soil':
        _data = hp_data.xs('ground', level='source').unstack()
        source_t = x - 5  # soil to brine temperature

    _sink_data = interp1d(
        _data.columns,
        _data.values,
        fill_value="extrapolate"
    )(sink_t)
    characteristic = interp1d(
        _data.dropna().index,
        _sink_data[~np.isnan(_sink_data)],
        fill_value="extrapolate"
    )(source_t)

    return characteristic


def cap_to_cap_ratio(x):
    """
    technology capacity is based on a reference source/sink temperature:
    2C source, 35C sink ashp
    0C source, 35C sink gshp
    We use this information to normalise the heating capacity at other source/sink
    combinations to 1.
    """
    if x.name == 'air':
        return x / x.xs((2, 35), level=('source_temp', 'sink_temp'))
    elif x.name == 'ground':
        return x / x.xs((0, 35), level=('source_temp', 'sink_temp'))


if __name__ == '__main__':
    get_characteristic(
        weather_pop=snakemake.input.weather_pop,
        hp_characteristics=snakemake.input.hp_characteristics,
        annual_demand=snakemake.input.annual_demand,
        heat_tech_params=snakemake.params.heat_tech_params,
        characteristic=snakemake.params.characteristic,
        tech=snakemake.params.tech,
        sink=snakemake.params.sink,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output[0]
    )