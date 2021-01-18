import pandas as pd


def get_hourly_ev_profiles(
    regions_path, ev_profiles_path, model_year, out_path
):
    """
    Fill empty countries and map EV profiles to national subregions.
    Profiles are already normalised relative to total number of vehicles in
    the fleet.
    """
    regions_df = pd.read_csv(regions_path).set_index(['id', 'country_code'])
    ev_profiles_df = (
        pd.read_csv(ev_profiles_path, index_col=[0, 1, 2], parse_dates=[0], squeeze=True)
        .xs(model_year, level='year')
    )

    def _fill_empty_country(df, country_neighbour_dict):
        df = df.unstack('country_code')
        for country, neighbours in country_neighbour_dict.items():
            df[country] = df[neighbours].mean(axis=1)
        return df.stack().rename('ev_profiles')

    # Fill missing countries based on nearest neighbours in the same timezone
    ev_profiles_df = _fill_empty_country(
        ev_profiles_df,
        {'ALB': ['HRV'], 'MKD': ['HRV'], 'GRC': ['ROU'], 'CYP': ['ROU'], 'BGR': ['ROU'],
         'BIH': ['HRV', 'HUN'], 'MNE': ['HRV'], 'ISL': ['GBR'], 'SRB': ['HUN']}
    )
    ev_profiles_df = (
        ev_profiles_df.align(regions_df)[0]
        .droplevel('country_code')
        .unstack('id')
    )
    # to naive timezone, to match all other CSVs in the model
    ev_profiles_df.index = ev_profiles_df.index.tz_localize(None)
    ev_profiles_df.to_csv(out_path)


if __name__ == "__main__":
    get_hourly_ev_profiles(
        regions_path=snakemake.input.regions,
        ev_profiles_path=snakemake.input.ev_profiles,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output[0],
    )
