import sys
import os

import pandas as pd
import numpy as np

sys.path.append('data/automatic/when2heat/')
from scripts import demand, read, misc


def get_heat_profiles(
    weather_pop, when2heat_path, out_path_sh, out_path_wh, out_path_h
):
    """
    Take annual demand per region and map it to hourly profiles, using the methods
    given in When2Heat. We partially use When2Heat methods, but also re-implement
    some since the original implementation does not handle subnational regions.
    """
    weather_pop_df = pd.read_csv(
        weather_pop, index_col=0, parse_dates=True, header=[0, 1, 2, 3, 4]
    ).rename_axis(columns={'id': 'country'})  # When2Heat classifies regions with 'country'

    population = (
        weather_pop_df
        .columns
        .to_frame()
        .population
        .astype(float)
        .groupby(['country', 'latitude', 'longitude']).mean()
    )
    weather_df = weather_pop_df.droplevel('population', axis=1)

    # to Kelvin, since that's what When2Heat expects
    air_temp_df = weather_df.xs('air_temp', level='dataset', axis=1) + 273.15

    # Only need site-wide mean wind speed for this analysis
    wind = weather_df.xs('wind_speed', level='dataset', axis=1).mean()

    # ref temp is based on the average temp from 3 days prior to each day in the timeseries
    reference_temperature = demand.reference_temperature(air_temp_df)

    # Parameters and how to apply them is based on [BDEW:2015]
    daily_params = read.daily_parameters(os.path.join(when2heat_path, 'input/'))
    hourly_params = read.hourly_parameters(os.path.join(when2heat_path, 'input/'))

    # Get daily demand
    daily_heat = demand.daily_heat(reference_temperature, wind, daily_params)
    daily_water = demand.daily_water(reference_temperature, wind, daily_params)

    # Map profiles to daily demand
    hourly_space, hourly_water, hourly_heat = get_demand_profiles(
        reference_temperature, daily_heat, daily_water, hourly_params
    )

    def _wavg(x):
        weights = population.reindex(x.columns.droplevel('building')).fillna(0)
        if weights.sum() == 0:
            weighted_avg = 0
        else:
            weighted_avg = np.average(x, weights=weights, axis=1)
        return pd.Series(data=weighted_avg, index=x.index)

    # Regionalise based on population and save
    for out_path, df in {
        out_path_sh: hourly_space, out_path_wh: hourly_water, out_path_h: hourly_heat
    }.items():
        df.groupby(level=['country', 'building'], axis=1).apply(_wavg).to_csv(out_path)


def get_demand_profiles(reference_temperature, daily_heat, daily_water, hourly_params):
    """
    Using the data from When2Heat, get demand profile shapes
    """

    def _hourly_profile(x, _ref, hourly_params, hr_idx, hr_weekday_idx):
        building_type = x.name[0]
        # Commercial profiles differ on weekday and weekends
        _idx = hr_weekday_idx if building_type == 'COM' else hr_idx
        return x.mul(hourly_params[building_type].lookup(_idx, _ref[x.name[1:]].values))

    heat_ref = misc.upsample_df(
        (np.ceil(((reference_temperature - 273.15) / 5).astype('float64')) * 5)
        .clip(lower=-15, upper=30),  # get temperature in 5C increments between -15C and +30C
        '60min'
    )
    water_ref = misc.upsample_df(
        pd.DataFrame(
            30,
            index=reference_temperature.index,
            columns=reference_temperature.columns
        ),  # hot water is based on the profile at 30C increment
        '60min'
    )

    upsampled_daily_heat = misc.upsample_df(daily_heat, '60min')
    upsampled_daily_water = misc.upsample_df(daily_water, '60min')

    for v in hourly_params.values():  # param DFs need numeric columns for later lookup
        v.columns = v.columns.astype(int)

    # Prepare lookup indexes
    hr_idx = upsampled_daily_heat.index.map(lambda i: i.strftime('%H:%M'))
    hr_weekday_idx = pd.MultiIndex.from_arrays(
        [upsampled_daily_heat.index.dayofweek, hr_idx],
        names=['weekday', 'time']
    )

    hourly_heat = upsampled_daily_heat.apply(
        _hourly_profile, _ref=heat_ref,
        hourly_params=hourly_params, hr_idx=hr_idx, hr_weekday_idx=hr_weekday_idx
    )
    hourly_water = upsampled_daily_water.apply(
        _hourly_profile, _ref=water_ref,
        hourly_params=hourly_params, hr_idx=hr_idx, hr_weekday_idx=hr_weekday_idx
    )
    hourly_space = (hourly_heat - hourly_water).clip(lower=0)

    return hourly_space, hourly_water, hourly_heat


if __name__ == "__main__":
    get_heat_profiles(
        weather_pop=snakemake.input.weather_pop,
        when2heat_path=snakemake.input.when2heat,
        out_path_sh=snakemake.output.space_heat,
        out_path_wh=snakemake.output.water_heat,
        out_path_h=snakemake.output.heat
    )
