import pandas as pd
import pycountry

import util

idx = pd.IndexSlice


def get_heat_demand(
    path_to_annual_demand, path_to_dwelling_ratio, path_to_profile, model_year, key, out_path
):
    """
    Take annual demand per region and use it to scale hourly profiles
    (generated from When2Heat for water heating and space heating, RAMP-cooking for cooking).
    """

    annual_demand_df = util.read_tdf(path_to_annual_demand)
    annual_demand_df = annual_demand_df.xs(model_year, level='year').droplevel('unit')
    dwelling_ratio_df = pd.read_csv(path_to_dwelling_ratio, index_col=0, squeeze=True)
    hourly_profile = pd.read_csv(path_to_profile, index_col=0, parse_dates=True, header=[0, 1])

    # Scale profiles to annual demand
    # (inc. shifting according to timezones and combining SFH and MFH profiles)
    if key in ["heat", "heat-bau-electricity"]:
        _key = [f"space_{key}", f"water_{key}"]
    else:
        _key = key.replace("-", "_")
    scaled_profile = _scale_profiles(
        hourly_profile, annual_demand_df, dwelling_ratio_df, model_year, _key
    )
    # Add industry space heating demand as a flat demand in all hours
    if key in ['space-heat', 'heat']:
        industry_space_heat = annual_demand_df.xs(
            ('industry_demand', 'industry', 'space_heat'),
            level=('dataset', 'cat_name', 'end_use')
        )
        scaled_profile -= industry_space_heat.div(len(scaled_profile.index))

    util.verify_profiles(scaled_profile, _key, annual_demand_df)

    # Save
    scaled_profile.to_csv(out_path)


def _scale_profiles(hourly_profile, annual_demand, dwelling_ratio, model_year, key):
    """
    Scale profile shapes according to total demand for space heating and hot water
    """

    def _shift_profiles(x):
        """
        Shift profiles forward/backward in time based on region timezones
        """
        country = pycountry.countries.lookup(x.name.split('_')[0]).alpha_2
        x = util.update_timeseries_timezone(x, country, model_year)

        return x

    def _get_profiles(df, end_use):
        def _scale_to_annual(df, cat_name):
            demand = annual_demand.loc[('heat_demand', cat_name, df.name, end_use)]
            if isinstance(demand, pd.Series):
                demand = demand.sum()
            return df.div(df.sum()).mul(demand)

        residential_hourly = (
            (df.xs('SFH', level='building', axis=1).mul(1 - dwelling_ratio))
            .add(df.xs('MFH', level='building', axis=1).mul(dwelling_ratio))
            .dropna(axis=1, how='all')
            .apply(_scale_to_annual, cat_name='household')
        )
        commercial_hourly = (
            df.xs('COM', level='building', axis=1)
            .apply(_scale_to_annual, cat_name='commercial')
        )
        print(df)
        print(residential_hourly.add(commercial_hourly))
        return residential_hourly.add(commercial_hourly).apply(_shift_profiles)

    # Calliope expects demand to be negative
    scaled_profile = (-1) * _get_profiles(hourly_profile, key)

    return scaled_profile


if __name__ == "__main__":
    get_heat_demand(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_dwelling_ratio=snakemake.input.dwelling_ratio,
        path_to_profile=snakemake.input.profile,
        model_year=snakemake.params.model_year,
        key=snakemake.params.key,
        out_path=snakemake.output[0],
    )
