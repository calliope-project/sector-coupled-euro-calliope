import xarray as xr
import pandas as pd
from scipy.interpolate import interp1d
import geopandas as gpd
import numpy as np

from hourly_heat_profiles import _map_pop_to_weather, _prep_weather_data
import util

idx = pd.IndexSlice


def get_cop(
    units, path_to_population_tif, soil_temp, air_temp,
    hp_characteristics, heat_tech_params, model_year, characteristic,
    gshp_water_heat_out_path, ashp_water_heat_out_path,
    gshp_space_heat_out_path, ashp_space_heat_out_path
):
    """
    """
    units_gdf = gpd.read_file(units)
    air_temp_ds = xr.open_dataset(air_temp)
    soil_temp_ds = xr.open_dataset(soil_temp)
    coords = air_temp_ds[['lat', 'lon']].to_dataframe()
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

    mapped_pop = _map_pop_to_weather(path_to_population_tif, coords, units_gdf)
    mapped_pop = mapped_pop.assign(
        longitude=mapped_pop.centroid.x, latitude=mapped_pop.centroid.y
    )

    air_temp_df = _prep_weather_data(
        air_temp_ds.temperature, model_year, mapped_pop
    )

    soil_temp_df = _prep_weather_data(
        soil_temp_ds.soil_temperture_5, model_year, mapped_pop
    ) - 273.15  # soil temperature is in K

    temperature = pd.concat([air_temp_df, soil_temp_df], axis=1, keys=['air', 'soil'], names=['source'])
    sink_t_dict = heat_tech_params['heating_temperature']

    site_data = pd.concat(
        [temperature.apply(temperature_to_characteristic, hp_data=hp_data_df, sink_t=v) for v in sink_t_dict.values()],
        keys=list(sink_t_dict.keys()), names=['sink'], axis=1
    ) * correction

    site_pop = mapped_pop.set_index(['latitude', 'longitude']).population

    def _pop_wavg(x, site_pop):
        to_drop = [i for i in x.columns.names if i not in ['latitude', 'longitude']]
        weighted_avg = np.average(x, weights=site_pop.reindex(x.droplevel(to_drop, axis=1).columns).fillna(0), axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    region_data = (
        site_data
        .groupby(level=['country', 'sink', 'source'], axis=1)
        .apply(_pop_wavg, site_pop=site_pop)
        .rename_axis(columns={'country': 'id'})
    )

    region_data.xs(('water', 'soil'), level=('sink', 'source'), axis=1).to_csv(gshp_water_heat_out_path)
    region_data.xs(('water', 'air'), level=('sink', 'source'), axis=1).to_csv(ashp_water_heat_out_path)

    heat_sink_ratio_dict = heat_tech_params['heat_sink_ratio']
    space_heat_data(region_data, heat_sink_ratio_dict, 'soil').to_csv(gshp_space_heat_out_path)
    space_heat_data(region_data, heat_sink_ratio_dict, 'air').to_csv(ashp_space_heat_out_path)


def space_heat_data(region_cop, heat_sink_ratio_dict, source):
    def _heat_sink_wavg(x, ratio_keys, ratio_vals):
        weighted_avg = np.average(x.loc[:, idx[:, ratio_keys]], weights=ratio_vals, axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    return (
        region_cop
        .xs(source, level='source', axis=1)
        .groupby('id', axis=1)
        .apply(_heat_sink_wavg,
               ratio_keys=[i for i in heat_sink_ratio_dict.keys()],
               ratio_vals=[i for i in heat_sink_ratio_dict.values()])
    )


def temperature_to_characteristic(x, hp_data, sink_t):
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
    if x.name == 'air':
        return x / x.xs((2, 35), level=('source_temp', 'sink_temp'))
    elif x.name == 'ground':
        return x / x.xs((0, 35), level=('source_temp', 'sink_temp'))


if __name__ == '__main__':
    get_cop(
        units=snakemake.input.units,
        path_to_population_tif=snakemake.input.path_to_population_tif,
        soil_temp=snakemake.input.soil_temp,
        air_temp=snakemake.input.air_temp,
        hp_characteristics=snakemake.input.hp_characteristics,
        heat_tech_params=snakemake.params.heat_tech_params,
        model_year=snakemake.params.model_year,
        characteristic=snakemake.params.characteristic,
        gshp_water_heat_out_path=snakemake.output.gshp_water_heat,
        ashp_water_heat_out_path=snakemake.output.ashp_water_heat,
        gshp_space_heat_out_path=snakemake.output.gshp_space_heat,
        ashp_space_heat_out_path=snakemake.output.ashp_space_heat
    )
