import pandas as pd
import numpy as np

import util


def update_hourly_electricity(
    space_heat, water_heat, cooking, public_transport, annual_demand,
    hourly_electricity, path_to_demand_scales, path_to_regions,
    first_year, final_year, out_path, scale_demand, demand_scale_scenario, projection_year

):
    profiles = {  # all profiles are negative values
        k: pd.read_csv(v, index_col=0, parse_dates=True)
        for k, v in {
            'space_heat': space_heat, 'water_heat': water_heat,
            'cooking': cooking, 'transport': public_transport,
        }.items()
    }
    hourly_electricity_df = pd.read_csv(hourly_electricity, index_col=0, parse_dates=True)
    print("There are {} missing data points, which will be backfilled".format(hourly_electricity_df.isna().sum().sum()))
    hourly_electricity_df = hourly_electricity_df.bfill()

    annual_demand_df = util.read_tdf(annual_demand)

    new_hourly_electricity_df = pd.concat([
        update_electricity_per_year(
            annual_demand_df.xs(year, level='year').droplevel('unit'),
            hourly_electricity_df.loc[str(year)],
            {k: v.loc[str(year)] for k, v in profiles.items()},
            scale_demand, demand_scale_scenario, path_to_demand_scales,
            projection_year, path_to_regions
        )
        for year in range(first_year, final_year + 1)
    ]).sort_index()

    # Communicate any timesteps which have ended up being positive
    # (i.e. BAU electricity is greater than demand for that hour...)
    # FIXME: this shouldn't happen (mostly occurs in Cyprus, e.g. 2% of its timesteps in 2016)
    positive_timesteps = new_hourly_electricity_df.where(lambda x: x > 0).dropna(axis=1, how='all').count()
    print(f"Following regions have some positive timesteps which will be clipped to zero:\n{positive_timesteps}")

    # Save
    new_hourly_electricity_df.clip(upper=0).to_csv(out_path)


def update_electricity_per_year(
    annual_demand_df, hourly_electricity_df, profiles,
    scale_demand, demand_scale_scenario, path_to_demand_scales, projection_year, path_to_regions
):
    new_hourly_electricity_df = hourly_electricity_df.copy()

    bau_electricity = (
        annual_demand_df
        .filter(regex='bau_electricity')
        .rename(lambda x: x.replace('_bau_electricity', '').replace('bau_', ''), level='end_use')
    )

    # Subtract old demand
    # Industry is a flatline and is 'added' since hourly electricity is negative
    new_hourly_electricity_df = new_hourly_electricity_df.add(
        bau_electricity
        .xs('industry_demand')
        .sum(level='id')
        .div(len(new_hourly_electricity_df.index))
    )
    # space heating, water heating, cooking, and non-freight transport is based on profiles
    # We use rail profiles for all non-freight transport, which is a simplification in light
    # of the small relative contribution of non-rail transport to electricity demand
    for profile_name, profile in profiles.items():
        if profile_name == 'transport':
            scaled_profile = profile.apply(lambda x: x.div(x.abs().sum()))
            # This profile is 'subtracted' as it is also negative
            new_hourly_electricity_df = new_hourly_electricity_df.sub(
                scaled_profile
                # fillna for regions with no rail transport, but with other transport electricity demand
                .fillna({i: scaled_profile.mean(axis=1) for i in profile.columns})
                .mul(bau_electricity.xs('transport_demand').sum(level='id'))
            )
        else:
            new_hourly_electricity_df = new_hourly_electricity_df.sub(profile)

    # Scale "current" demand to future demand from commercial+residential buildings, if requested
    if scale_demand:
        print(f"Scaling demands with scenario {demand_scale_scenario} and projection year {projection_year}")
        demand_scales = (
            util.read_tdf(path_to_demand_scales)
            .xs((demand_scale_scenario, projection_year), level=("scenario", "year"))
            .rename_axis(index={"id": "country_code"})
            .align(pd.read_csv(path_to_regions).set_index(["id", "country_code"]))[0]
            .droplevel("country_code")
        )

        new_hourly_electricity_df = new_hourly_electricity_df.mul(demand_scales, axis=1).dropna()
    else:
        demand_scales = 1

    # Add back in new demand
    # Industry, this time 'subtracted' to make hourly electricity more negative
    new_hourly_electricity_df = (
        new_hourly_electricity_df.sub(
            annual_demand_df
            .xs('industry_demand')
            .xs('electricity', level='end_use')
            .sum(level='id')
            .div(len(new_hourly_electricity_df.index))
        )
    )
    # Rail transport (profile is negative, so we add it)
    new_hourly_electricity_df = new_hourly_electricity_df.add(profiles['transport'])

    # Verify on an annual level
    assert np.allclose(
        hourly_electricity_df.sum(),
        new_hourly_electricity_df.sum()
        .add(annual_demand_df.xs('electricity', level='end_use').sum(level='id'))
        .div(demand_scales)
        .sub(bau_electricity.sum(level='id'))
        .reindex(hourly_electricity_df.columns)
    )
    return new_hourly_electricity_df


if __name__ == "__main__":
    update_hourly_electricity(
        space_heat=snakemake.input.space_heat,
        water_heat=snakemake.input.water_heat,
        cooking=snakemake.input.cooking,
        public_transport=snakemake.input.public_transport,
        annual_demand=snakemake.input.annual_demand,
        hourly_electricity=snakemake.input.hourly_electricity,
        path_to_demand_scales=snakemake.input.demand_scales,
        path_to_regions=snakemake.input.regions,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        scale_demand=snakemake.params.scale_demand,
        demand_scale_scenario=snakemake.params.demand_scale_scenario,
        projection_year=int(snakemake.params.projection_year),
        out_path=snakemake.output[0]
    )
