import pandas as pd
import util


def get_transport_profiles(
    regions, annual_demand, rail_profiles, model_year, out_path
):
    """
    Take DeStInEe rail profiles for weekday, saturday, and sunday and create rail demand
    profiles for all eurospores regions, based on specific year and country timezone.
    We assume that these profiles are valid for rail and electric buses(?).
    """
    annual_demand_df = (
        util.read_tdf(annual_demand)
        .unstack('id')
        .droplevel('unit')
        .xs(('transport_demand', 'rail', model_year))
    )

    regions_df = pd.read_csv(regions, index_col=0, squeeze=True)
    profiles = pd.read_csv(rail_profiles, index_col=0, header=0)
    profiles.columns = profiles.columns.astype(float)
    profiles.index = profiles.index.astype(float)

    annual_profile = pd.Series(index=pd.date_range(str(model_year), str(model_year + 1), freq='1H'))[str(model_year)]
    annual_profile = pd.merge(
        annual_profile
        .groupby([annual_profile.index.date, annual_profile.index.hour, annual_profile.index.dayofweek]).sum()
        .rename_axis(index=['date', 'hour', 'dayofweek'])
        .rename('drop_me'),
        profiles.dropna(axis=1).stack()
        .rename_axis(index=['hour', 'dayofweek'])
        .rename('profiles'),
        left_index=True, right_index=True
    ).drop('drop_me', axis=1)
    annual_profile.index = (
        pd.to_datetime(annual_profile.index.get_level_values('date')) +
        pd.to_timedelta(annual_profile.index.get_level_values('hour').astype(int), unit='h')
    )
    annual_profile = annual_profile.sort_index()
    annual_profile = pd.concat(
        [annual_profile.div(annual_profile.sum())['profiles'] for i in regions_df.index],
        axis=1, keys=regions_df.index, names=['id']
    )

    def _scale_and_shift(x):
        return (-1) * util.update_timeseries_timezone(
            x.mul(annual_demand_df.loc['electricity', x.name]),
            util.get_alpha2(regions_df.loc[x.name], eurostat=False),
            model_year
        )

    annual_profile = annual_profile.apply(_scale_and_shift)

    # Save
    annual_profile.to_csv(out_path)


if __name__ == "__main__":
    get_transport_profiles(
        regions=snakemake.input.regions,
        annual_demand=snakemake.input.annual_demand,
        rail_profiles=snakemake.input.rail_profiles,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output[0]
    )
