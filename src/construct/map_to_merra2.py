import geopandas as gpd
import pandas as pd
import xarray as xr
from rasterstats import zonal_stats
import rasterio

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
WGS84 = "EPSG:4326"


def map_to_merra2(population, air_temp, units, temperature_threshold, out_path):
    # Get population per merra-2 site id and country
    europe_shape = gpd.read_file(units)
    if isinstance(europe_shape.crs, dict):
        europe_shape.crs = europe_shape.crs['init']  # {'init': '...'} is deprecated

    air_temp_ds = xr.open_dataset(air_temp)
    coords = air_temp_ds[['lat', 'lon']].to_dataframe()
    temperatures = air_temp_ds['temperature']
    annual_hdd = get_hdd(temperatures, temperature_threshold)

    points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(coords.lon.values, coords.lat.values),
        index=coords.index, crs=WGS84
    )
    points_m = points.to_crs(EPSG_3035_PROJ4)
    polys_m = gpd.GeoDataFrame(
        index=points_m.index, geometry=points_m.buffer(25000).envelope  # 50km horizontal resolution
    )
    polys_eu = gpd.overlay(polys_m.to_crs(WGS84).reset_index(), europe_shape.to_crs(WGS84))

    with rasterio.open(population) as src:
        array = src.read(1)
        crs = src.crs
        affine = src.transform

        pop_polys = zonal_stats(polys_eu.to_crs(crs), array, affine=affine, stats='sum', nodata=0)
        polys_eu['population'] = [i['sum'] for i in pop_polys]

        pop_eu = zonal_stats(europe_shape.to_crs(crs), array, affine=affine, stats='sum', nodata=0)
        europe_shape['population'] = [i['sum'] for i in pop_eu]

    sig_fig = lambda i: f'{float(f"{i:.4g}"):g}'
    assert (
        sig_fig(polys_eu['population'].sum()) ==
        sig_fig(europe_shape['population'].sum())
    )
    pop_polys_with_hdd = pd.concat(
        [polys_eu.set_index('site')[['country_code', 'population']], annual_hdd],
        axis=1
    )
    pop_polys_with_hdd.to_csv(out_path)


def get_hdd(air_temp, temperature_threshold):
    hdh = air_temp - temperature_threshold  # hdh = heating degree hour
    hdd = hdh.where(hdh >= 0, other=0) / 24  # div by 24 for heating degree hour -> heating degree day

    annual_hdd = hdd.groupby('time.year').sum('time').to_pandas()

    return annual_hdd


if __name__ == '__main__':
    map_to_merra2(
        population=snakemake.input.population,
        air_temp=snakemake.input.air_temp,
        units=snakemake.input.units,
        temperature_threshold=snakemake.params.temperature_threshold,
        out_path=snakemake.output

    )