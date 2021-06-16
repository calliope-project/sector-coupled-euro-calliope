import geopandas as gpd


def regions_gdf_to_csv(inpath, outpath):
    units_gdf = gpd.read_file(inpath)
    units_gdf.set_index("id").country_code.to_csv(outpath)


if __name__ == "__main__":
    regions_gdf_to_csv(snakemake.input.units, snakemake.output[0])
