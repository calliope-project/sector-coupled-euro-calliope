import math

import pandas as pd
import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rasterio

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
WGS84 = "EPSG:4326"
idx = pd.IndexSlice


def gridded_weather_population(
    population, units, air_temp, wind_10m, soil_temp, model_year,
    weather_out_path, regions_out_path
):
    units_gdf = gpd.read_file(units)
    air_temp_ds = xr.open_dataset(air_temp)
    wind_speed_ds = xr.open_dataset(wind_10m)
    soil_temp_ds = xr.open_dataset(soil_temp)
    coords = air_temp_ds[['lat', 'lon']].to_dataframe()

    mapped_pop = _map_pop_to_weather(population, coords, units_gdf)
    mapped_pop = mapped_pop.assign(
        longitude=mapped_pop.centroid.x, latitude=mapped_pop.centroid.y
    )

    air_temp_df = _prep_weather_data(
        air_temp_ds.temperature, model_year, mapped_pop
    )

    # Only need site-wide mean wind speed for this analysis
    wind_speed_df = _prep_weather_data(
        wind_speed_ds.wind_speed, model_year, mapped_pop
    )

    soil_temp_df = _prep_weather_data(
        soil_temp_ds.soil_temperture_5, model_year, mapped_pop
    ) - 273.15  # soil temperature is in K

    weather = pd.concat(
        [air_temp_df, soil_temp_df, wind_speed_df],
        keys=['air_temp', 'soil_temp', 'wind_speed'],
        names=['dataset'], axis=1
    )

    units_gdf.set_index("id").country_code.to_csv(regions_out_path)
    weather.to_csv(weather_out_path)


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
    xarray dataarray to pandas df with id, lat, lon, population as columns
    and hourly timeseries as index.
    """
    _df = dataarray.loc[{'time': str(model_year)}].to_pandas()
    _df = _df.loc[mapped_pop.index].assign(
        latitude=mapped_pop.latitude,
        longitude=mapped_pop.longitude,
        id=mapped_pop.id,
        population=mapped_pop.population
    ).set_index(['id', 'latitude', 'longitude', 'population'])
    _df.columns = pd.to_datetime(_df.columns)
    _df = _df[(_df.fillna(0) != 0).any(axis=1)]  # remove sites with no data (seems to be logged as all zeros)

    return _df.T


if __name__ == "__main__":
    gridded_weather_population(
        population=snakemake.input.population,
        units=snakemake.input.units,
        air_temp=snakemake.input.air_temp,
        wind_10m=snakemake.input.wind_10m,
        soil_temp=snakemake.input.soil_temp,
        model_year=snakemake.params.model_year,
        weather_out_path=snakemake.output.weather_pop,
        regions_out_path=snakemake.output.regions,
    )
