import geopandas as gpd
import xarray as xr
from rasterstats import zonal_stats
import rasterio

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
WGS84 = "EPSG:4326"


def population_map(population, merra_coords, units, out_path):
    # Get population per merra-2 site id and country
    europe_shape = gpd.read_file(units)
    merra_coords = xr.open_dataset(merra_coords)
    coords = merra_coords[['lat', 'lon']].to_dataframe()
    points = gpd.GeoDataFrame(
        geometry=[shapely.geometry.Point(xy) for xy in zip(coords.lon.values, coords.lat.values)],
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
    polys_eu[['site', 'country_code', 'population']].to_csv(out_path)


if __name__ == '__main__':
    population_map(
        population=snakemake.inputs.population,
        merra_coords=snakemake.inputs.soil_temp,
        units=snakemake.inputs.units,
        out_path=snakemake.output
    )