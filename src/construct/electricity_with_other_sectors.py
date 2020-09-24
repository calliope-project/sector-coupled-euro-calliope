import pandas as pd
import numpy as np

import util


def update_hourly_electricity(
    space_heating, water_heating, cooking, public_transport, annual_demand,
    hourly_electricity, model_year, out_path
):
    profiles = {
        k: pd.read_csv(v, index_col=0, parse_dates=True)
        for k, v in {
            'space_heating': space_heating, 'water_heating': water_heating,
            'cooking': cooking, 'transport': public_transport,

        }.items()
    }
    hourly_electricity_df = pd.read_csv(hourly_electricity, index_col=0, parse_dates=True)
    new_hourly_electricity_df = hourly_electricity_df.copy()
    annual_demand_df = util.read_tdf(annual_demand)
    annual_demand_df = annual_demand_df.xs(model_year, level='year').droplevel('unit')

    bau_electricity = (
        annual_demand_df
        .filter(regex='bau_electricity')
        .rename(lambda x: x.replace('_bau_electricity', '').replace('bau_', ''), level='end_use')
    )

    # Subtract old demand
    # Industry is a flatline and is 'added' since hourly electricity is negative
    new_hourly_electricity_df = new_hourly_electricity_df.add(
        bau_electricity.sum(level='id').div(len(new_hourly_electricity_df.index))
    )
    # space heating, water heating, cooking, and non-freight transport is based on profiles
    # We use rail profiles for all non-freight transport, which is a simplification in light
    # of the small relative constribution of non-rail transport to electricity demand
    for profile_name, profile in profiles.items():
        if profile_name == 'transport':
            profile_shape = profile.apply(_scale_to_annual)
            profile = profile_shape.mul(bau_electricity.xs('transport_demand').sum(level='id'))
        new_hourly_electricity_df = new_hourly_electricity_df.sub(profile)

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
    # Rail transport
    new_hourly_electricity_df = new_hourly_electricity_df.add(profiles['transport'])

    # Verify on an annual level
    assert np.allclose(
        hourly_electricity_df.sum(),
        new_hourly_electricity_df.sum() -
        bau_electricity +
        annual_demand_df.xs('electricity', level='end_use').sum(level='id')
    )

    # Save
    new_hourly_electricity_df.to_csv(out_path)


def _scale_to_annual(x):
    return x.div(x.abs().sum())


if __name__ == "__main__":
    update_hourly_electricity(
        space_heating=snakemake.input.space_heating,
        water_heating=snakemake.input.water_heating,
        cooking=snakemake.input.cooking,
        public_transport=snakemake.input.public_transport,
        annual_demand=snakemake.input.annual_electricity,
        hourly_electricity=snakemake.input.hourly_electricity,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output[0]
    )
