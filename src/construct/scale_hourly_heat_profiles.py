import pandas as pd
import numpy as np
import geopandas as gpd
import pycountry

import util

idx = pd.IndexSlice


def get_heat_demand(
    units, annual_demand, dwellings, model_year,
    path_to_sh_profiles, path_to_wh_profiles, path_to_c_profiles,
    nuts_to_regions, sh_key, wh_key, c_key,
    out_path_sh, out_path_wh, out_path_c
):
    """
    Take annual demand per region and use it to scale hourly profiles
    (generated from When2Heat for water heating and space heating, RAMP-cooking for cooking).
    """
    sh_key = sh_key.replace("-", "_")
    wh_key = wh_key.replace("-", "_")
    c_key = c_key.replace("-", "_")

    units_gdf = gpd.read_file(units)

    annual_demand_df = util.read_tdf(annual_demand)
    annual_demand_df = annual_demand_df.xs(model_year, level='year').droplevel('unit')

    dwelling_ratio = _get_dwelling_ratio(dwellings, units_gdf, nuts_to_regions)

    hourly_space_profiles = pd.read_csv(path_to_sh_profiles, index_col=0, parse_dates=True, header=[0, 1])
    hourly_water_profiles = pd.read_csv(path_to_wh_profiles, index_col=0, parse_dates=True, header=[0, 1])
    hourly_cooking_profiles = pd.read_csv(
        path_to_c_profiles, index_col=[0, 1, 2], parse_dates=[0], squeeze=True
    ).xs(model_year, level='year')

    # Scale profiles to annual demand
    # (inc. shifting according to timezones and combining SFH and MFH profiles)
    scaled_hourly_space, scaled_hourly_water = scale_space_water_profiles(
        hourly_space_profiles, hourly_water_profiles, annual_demand_df, dwelling_ratio,
        model_year, sh_key, wh_key
    )
    # Add industry space heating demand as a flat demand in all hours
    if sh_key == 'space_heat':
        industry_space_heat = annual_demand_df.loc[('industry_demand', 'industry')].unstack()['space_heat']
        scaled_hourly_space -= industry_space_heat.div(len(scaled_hourly_space.index))

    # Get cooking profiles using a different dataset (RAMP-cooking output)
    scaled_hourly_cooking = get_hourly_cooking(
        annual_demand_df, hourly_cooking_profiles, units_gdf, model_year, c_key
    )

    verify_profiles(
        {sh_key: scaled_hourly_space, wh_key: scaled_hourly_water, c_key: scaled_hourly_cooking},
        annual_demand=annual_demand_df
    )

    # Save
    scaled_hourly_space.to_csv(out_path_sh)
    scaled_hourly_water.to_csv(out_path_wh)
    scaled_hourly_cooking.to_csv(out_path_c)


def scale_space_water_profiles(
    hourly_space, hourly_water, annual_demand, dwelling_ratio, model_year, sh_key, wh_key
):
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
            return (
                df
                .div(df.sum())
                .mul(annual_demand.loc[('heat_demand', cat_name, df.name, end_use)])
            )

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
        return residential_hourly.add(commercial_hourly).apply(_shift_profiles)

    # Calliope expects demand to be negative
    demand_space = (-1) * _get_profiles(hourly_space, sh_key)
    demand_water = (-1) * _get_profiles(hourly_water, wh_key)

    return demand_space, demand_water


def _get_dwelling_ratio(path_to_dwellings, units, nuts_to_regions):
    """
    Per region, get the ratio of single family homes to multi-family homes.
    Although the underlying data is valid for (2010?) NUTS regions, we just use
    the national average applied to each region in a country.
    Missing locations (e.g. croatia) get the average across all other datasets
    """

    dwellings = pd.read_csv(path_to_dwellings, delimiter='\t', index_col=0)
    dwellings.index = dwellings.index.str.split(',', expand=True).rename(
        ['housing', 'building_type', 'unit', 'year']
    )
    dwellings.columns = dwellings.columns.str.strip().rename('region')

    dwellings = dwellings.droplevel(['unit', 'year']).apply(util.to_numeric)

    nuts_2010 = pd.read_excel(nuts_to_regions, sheet_name='locations')
    nuts_2010['NUTS3_2010'] = nuts_2010["NUTS3_2010"].fillna(nuts_2010.Country.where(nuts_2010.Source != 'NUTS3'))
    nuts_2010 = nuts_2010.dropna(subset=['NUTS3_2010']).set_index('NUTS3_2010').EuroSPORES

    # group dwelling data into eurospores regions
    dwellings_eurospores = dwellings.groupby(nuts_2010, axis=1).sum()
    SFH_eurospores = dwellings_eurospores.xs(("DW", "RES1"))
    MFH_eurospores = dwellings_eurospores.loc[idx["DW", ["RES2", "RES_GE3"]], :].sum()

    MFH_to_SFH = MFH_eurospores / (SFH_eurospores + MFH_eurospores)
    MFH_to_SFH = MFH_to_SFH.reindex(units.id)

    # fill empties
    print(
        'Filling missing regions {} with average dwelling ratio'
        .format(MFH_to_SFH[MFH_to_SFH.isnull()].index.values)
    )
    MFH_to_SFH = MFH_to_SFH.fillna(MFH_to_SFH.mean())

    return MFH_to_SFH


def get_hourly_cooking(annual_demand, cooking_profiles, units, model_year, c_key):
    def _fill_empty_country(df, country_neighbour_dict):
        df = df.unstack('country_code')
        for country, neighbours in country_neighbour_dict.items():
            df[country] = df[neighbours].mean(axis=1)
        return df.stack().rename('cooking_profiles')

    id_to_country = units.set_index('id').country_code
    # We merge commercial and household cooking demands together
    annual_cooking = (
        annual_demand
        .xs(('heat_demand', c_key), level=('dataset', 'end_use'))
        .sum(level='id')
    )
    annual_cooking = (
        pd.merge(
            annual_cooking.to_frame('annual_cooking'),
            id_to_country.to_frame('country_code'),
            left_index=True, right_index=True
        )
        .set_index('country_code', append=True)
    )
    # Fill missing countries based on nearest neighbours in the same timezone
    cooking_profiles = _fill_empty_country(
        cooking_profiles,
        {'ALB': ['SRB'], 'MKD': ['SRB'], 'GRC': ['BGR'], 'CYP': ['BGR'],
         'BIH': ['SRB', 'HRV'], 'MNE': ['SRB', 'HRV'], 'ISL': ['GBR']}
    )
    df = pd.merge(cooking_profiles, annual_cooking, left_index=True, right_index=True)
    hourly_cooking = df['annual_cooking'].mul(df['cooking_profiles'])
    hourly_cooking = hourly_cooking.droplevel('country_code').unstack()

    return (-1) * hourly_cooking


def verify_profiles(profiles, annual_demand):
    for key, df in profiles.items():
        assert (df <= 0).all().all()
        assert np.allclose(
            -1 * df.sum(),
            annual_demand.xs(key, level='end_use').sum(level='id')
        )


if __name__ == "__main__":
    get_heat_demand(
        units=snakemake.input.units,
        annual_demand=snakemake.input.annual_demand,
        dwellings=snakemake.input.dwellings,
        path_to_sh_profiles=snakemake.input.space_profiles,
        path_to_wh_profiles=snakemake.input.water_profiles,
        path_to_c_profiles=snakemake.input.cooking_profiles,
        nuts_to_regions=snakemake.input.nuts_to_regions,
        model_year=snakemake.params.model_year,
        sh_key=snakemake.params.space_heat_key,
        wh_key=snakemake.params.water_heat_key,
        c_key=snakemake.params.cooking_key,
        out_path_sh=snakemake.output.space_heat,
        out_path_wh=snakemake.output.water_heat,
        out_path_c=snakemake.output.cooking,
    )