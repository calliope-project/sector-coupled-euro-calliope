import sys
import os
import math

import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import pycountry
import pytz
from rasterstats import zonal_stats
import rasterio

import util

sys.path.append('data/automatic/when2heat/')
from scripts import demand, read, preprocess, misc


EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
WGS84 = "EPSG:4326"
idx = pd.IndexSlice


def get_heat_profiles(
    population, units, air_temp, wind_10m, model_year, when2heat_path,
    out_path_sh, out_path_wh
):
    """
    Take annual demand per region and map it to hourly profiles, using the methods
    given in When2Heat. We partially use When2Heat methods, but also re-implement
    some since the original implementation does not handle subnational regions.
    """
    units_gdf = gpd.read_file(units)
    air_temp_ds = xr.open_dataset(air_temp)
    wind_speed_ds = xr.open_dataset(wind_10m)
    coords = air_temp_ds[['lat', 'lon']].to_dataframe()

    mapped_pop = _map_pop_to_weather(population, coords, units_gdf)
    mapped_pop = mapped_pop.assign(
        longitude=mapped_pop.centroid.x, latitude=mapped_pop.centroid.y
    )

    air_temp_df = _prep_weather_data(
        air_temp_ds.temperature, model_year, mapped_pop
    ) + 273.15  # to Kelvin, since that's what When2Heat expects

    # Only need site-wide mean wind speed for this analysis
    wind_speed_df = _prep_weather_data(
        wind_speed_ds.wind_speed, model_year, mapped_pop
    )
    wind = wind_speed_df.mean()

    # ref temp is based on the average temp from 3 days prior to each day in the timeseries
    reference_temperature = demand.reference_temperature(air_temp_df)

    # Parameters and how to apply them is based on [BDEW:2015]
    daily_params = read.daily_parameters(os.path.join(when2heat_path, 'input/'))
    hourly_params = read.hourly_parameters(os.path.join(when2heat_path, 'input/'))

    # Get daily demand
    daily_heat = demand.daily_heat(reference_temperature, wind, daily_params)
    daily_water = demand.daily_water(reference_temperature, wind, daily_params)

    # Map profiles to daily demand
    hourly_space, hourly_water = get_demand_profiles(
        reference_temperature, daily_heat, daily_water, hourly_params
    )

    # combine profiles within a region based on a population weighted average
    hourly_space, hourly_water = regional_profiles(
        hourly_space, hourly_water, mapped_pop
    )

    # Save
    hourly_space.to_csv(out_path_sh)
    hourly_water.to_csv(out_path_wh)


def get_demand_profiles(reference_temperature, daily_heat, daily_water, hourly_params):
    """
    Using the data from When2Heat, get demand profile shapes
    """

    def _hourly_profile(x, _ref, hourly_params, hr_idx, hr_weekday_idx):
        building_type = x.name[0]
        # Commercial profiles differ on weekday and weekends
        _idx = hr_weekday_idx if building_type == 'COM' else hr_idx
        return x.mul(hourly_params[building_type].lookup(_idx, _ref[x.name[1:]].values))

    heat_ref = misc.upsample_df(
        (np.ceil(((reference_temperature - 273.15) / 5).astype('float64')) * 5)
        .clip(lower=-15, upper=30),  # get temperature in 5C increments between -15C and +30C
        '60min'
    )
    water_ref = misc.upsample_df(
        pd.DataFrame(
            30,
            index=reference_temperature.index,
            columns=reference_temperature.columns
        ),  # hot water is based on the profile at 30C increment
        '60min'
    )

    upsampled_daily_heat = misc.upsample_df(daily_heat, '60min')
    upsampled_daily_water = misc.upsample_df(daily_water, '60min')

    for v in hourly_params.values():  # param DFs need numeric columns for later lookup
        v.columns = v.columns.astype(int)

    # Prepare lookup indexes
    hr_idx = upsampled_daily_heat.index.map(lambda i: i.strftime('%H:%M'))
    hr_weekday_idx = pd.MultiIndex.from_arrays(
        [upsampled_daily_heat.index.dayofweek, hr_idx],
        names=['weekday', 'time']
    )

    hourly_heat = upsampled_daily_heat.apply(
        _hourly_profile, _ref=heat_ref,
        hourly_params=hourly_params, hr_idx=hr_idx, hr_weekday_idx=hr_weekday_idx
    )
    hourly_water = upsampled_daily_water.apply(
        _hourly_profile, _ref=water_ref,
        hourly_params=hourly_params, hr_idx=hr_idx, hr_weekday_idx=hr_weekday_idx
    )
    hourly_space = (hourly_heat - hourly_water).clip(lower=0)

    return hourly_space, hourly_water


def _map_pop_to_weather(population, coords, units):
    # Get population per merra2 site id and country

    points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(coords.lon.values, coords.lat.values),
        index=coords.index, crs=WGS84
    )
    points_m = points.to_crs(EPSG_3035_PROJ4)
    polys_m = gpd.GeoDataFrame(
        index=points_m.index, geometry=points_m.buffer(25000).envelope  # 50km horizontal resolution
    )
    polys_eu = gpd.overlay(polys_m.to_crs(WGS84).reset_index(), units.to_crs(WGS84))

    with rasterio.open(population) as src:
        array = src.read(1)
        crs = src.crs
        affine = src.transform

        pop_polys = zonal_stats(polys_eu.to_crs(crs), array, affine=affine, stats='sum', nodata=0)
        polys_eu['population'] = [i['sum'] for i in pop_polys]

        pop_eu = zonal_stats(units.to_crs(crs), array, affine=affine, stats='sum', nodata=0)
        units['population'] = [i['sum'] for i in pop_eu]

    # confirm that the total population is valid (i.e. we haven't picked up or lost regions)
    assert math.isclose(
        units.population.sum(), polys_eu.population.sum(), abs_tol=10**3
    )

    return polys_eu.set_index('site').drop(columns=['name', 'type', 'proper'], errors='ignore')


def _prep_weather_data(dataarray, model_year, mapped_pop):
    """
    xarray dataarray to pandas df with id, lat and lon as index.
    'id' is renamed to 'country' to work with when2heat
    """
    _df = dataarray.loc[{'time': str(model_year)}].to_pandas()
    _df = _df.loc[mapped_pop.index].assign(
        latitude=mapped_pop.latitude,
        longitude=mapped_pop.longitude,
        country=mapped_pop.id
    ).set_index(['country', 'latitude', 'longitude'])
    _df.columns = pd.to_datetime(_df.columns)
    _df = _df[(_df.fillna(0) != 0).any(axis=1)]  # remove sites with no data (seems to be logged as all zeros)

    return _df.T


def regional_profiles(hourly_space, hourly_water, mapped_pop):

    site_pop = mapped_pop.set_index(['latitude', 'longitude']).population

    def _wavg(x):
        weights = site_pop.reindex(x.droplevel((0, 1), axis=1).columns).fillna(0)
        if weights.sum() == 0:
            weighted_avg = 0
        else:
            weighted_avg = np.average(x, weights=weights, axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    regional_results = []
    for df in [hourly_space, hourly_water]:
        regional_results.append(
            df.groupby(level=['country', 'building'], axis=1).apply(_wavg)
        )

    return regional_results


if __name__ == "__main__":
    get_heat_profiles(
        population=snakemake.input.population,
        units=snakemake.input.units,
        air_temp=snakemake.input.air_temp,
        wind_10m=snakemake.input.wind_10m,
        model_year=snakemake.params.model_year,
        when2heat_path=snakemake.input.when2heat,
        out_path_sh=snakemake.output.space_heat,
        out_path_wh=snakemake.output.water_heat
    )
