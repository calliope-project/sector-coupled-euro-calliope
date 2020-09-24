import pandas as pd
import geopandas as gpd


def get_hourly_ev_profiles(
    units, ev_profiles, model_year, out_path
):
    """
    Fill empty countries and map EV profiles to national subregions.
    Profiles are already normalised relative to total number of vehicles in
    the fleet.
    """
    units_gdf = gpd.read_file(units)
    ev_profiles_df = pd.read_csv(ev_profiles, index_col=0, parse_dates=0)

    def _fill_empty_country(df, country_neighbour_dict):
        df = df.unstack('country_code')
        for country, neighbours in country_neighbour_dict.items():
            df[country] = df[neighbours].mean(axis=1)
        return df.stack().rename('cooking_profiles')

    id_to_country = units_gdf.set_index('country_code').id

    # Fill missing countries based on nearest neighbours in the same timezone
    ev_profiles_df = _fill_empty_country(
        ev_profiles_df,
        {'ALB': ['SRB'], 'MKD': ['SRB'], 'GRC': ['BGR'], 'CYP': ['BGR'],
         'BIH': ['SRB', 'HRV'], 'MNE': ['SRB', 'HRV'], 'ISL': ['GBR']}
    )
    ev_profiles_df = ev_profiles_df.reindex(id_to_country, axis=1)

    ev_profiles_df.to_csv(out_path)


if __name__ == "__main__":
    get_hourly_ev_profiles(
        units=snakemake.input.units,
        ev_profiles=snakemake.input.ev_profiles,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output[0],
    )
