import pandas as pd


def get_hourly_ev_profiles(
    regions_path, ev_profiles_path, dataset_name, demand_range, first_year, final_year, out_path
):
    """
    Fill empty countries and map EV profiles to national subregions.
    Profiles are already normalised relative to total number of vehicles in
    the fleet.
    """
    regions_df = pd.read_csv(regions_path).set_index(['id', 'country_code'])
    ev_profiles_df = pd.read_csv(
        ev_profiles_path, index_col=[0, 1, 2], parse_dates=[0], squeeze=True
    )
    profiles = []
    for year in range(first_year, final_year + 1):
        profiles.append(get_one_year_hourly_ev_profiles(
            dataset_name, demand_range, regions_df,
            ev_profiles_df.xs(year, level='year')
        ))
    pd.concat(profiles, sort=True).sort_index().to_csv(out_path)


def get_one_year_hourly_ev_profiles(
    dataset_name, demand_range, regions_df, ev_profiles_df
):

    # Demand is normalised to just get the fluctuations in demand rather than absolute values
    # We create two demand profiles, to create a min/max range of total demand that can be met
    # in a timestep (or collection of timesteps)
    if "demand" in dataset_name:
        if "light" in dataset_name:
            ev_profile = (
                ev_profiles_df
                .div(ev_profiles_df.sum(level='country_code'))
                .mul(demand_range[dataset_name.split("-")[1]])
            )
        elif "heavy" in dataset_name:  # assume demand is equal across all timesteps for freight and buses
            ev_profile = (
                ev_profiles_df
                .clip(upper=1, lower=1)
                .div(len(ev_profiles_df.index.get_level_values('datetime').unique()))
                .mul(demand_range[dataset_name.split("-")[1]])
            )
    # % plugged-in EVs is already normalised
    elif dataset_name == 'plugin':
        ev_profile = ev_profiles_df

    def _fill_empty_country(df, country_neighbour_dict):
        df = df.unstack('country_code')
        for country, neighbours in country_neighbour_dict.items():
            df[country] = df[neighbours].mean(axis=1)
        return df.stack().rename('ev_profiles')

    # Fill missing countries based on nearest neighbours in the same timezone
    ev_profile = _fill_empty_country(
        ev_profile,
        {'ALB': ['HRV'], 'MKD': ['HRV'], 'GRC': ['ROU'], 'CYP': ['ROU'], 'BGR': ['ROU'],
         'BIH': ['HRV', 'HUN'], 'MNE': ['HRV'], 'ISL': ['GBR'], 'SRB': ['HUN']}
    )
    ev_profile = (
        ev_profile.align(regions_df)[0]
        .droplevel('country_code')
        .unstack('id')
    )
    # to naive timezone, to match all other CSVs in the model
    ev_profile.index = ev_profile.index.tz_localize(None)
    return ev_profile


if __name__ == "__main__":
    get_hourly_ev_profiles(
        regions_path=snakemake.input.regions,
        ev_profiles_path=snakemake.input.ev_profiles,
        dataset_name=snakemake.wildcards.dataset_name,
        demand_range=snakemake.params.demand_range,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        out_path=snakemake.output[0],
    )
