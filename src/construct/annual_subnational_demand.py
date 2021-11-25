import pandas as pd
import geopandas as gpd
import numpy as np

from annual_national_demand import electrify_industry_demand
import util


FREIGHT_SECTORS = {
    'GT03': 'Mining and quarrying',
    'GT04': 'Food, beverages and tobacco',
    'GT05': 'Textiles and leather',
    'GT06': 'Wood and wood products',
    'GT11': 'Machinery Equipment',
    'GT12': 'Transport Equipment',
    'GT13': 'Other Industrial Sectors'
}

CH_CANTON_TO_NUTS3 = {
    'AG': 'CH033', 'AI': 'CH054', 'AR': 'CH053', 'BE': 'CH021',
    'BL': 'CH032', 'BS': 'CH031', 'FR': 'CH022', 'GE': 'CH013',
    'GL': 'CH051', 'GR': 'CH056', 'JU': 'CH025', 'LU': 'CH061',
    'NE': 'CH024', 'NW': 'CH065', 'OW': 'CH064', 'SG': 'CH055',
    'SH': 'CH052', 'SO': 'CH023', 'SZ': 'CH063', 'TG': 'CH057',
    'TI': 'CH070', 'UR': 'CH062', 'VD': 'CH011', 'VS': 'CH012',
    'ZG': 'CH066', 'ZH': 'CH040'
}
LEVEL_ORDER = ['country_code', 'year', 'id', 'unit', 'end_use']
YEAR_RANGE = range(2000, 2019)  # FIXME: this should be in config

idx = pd.IndexSlice


def subnationalise_demand(
    population, units, building_demand, building_heat_electricity_consumption,
    industry_demand, road_distance, road_vehicles, rail_demand, air_demand, marine_demand,
    road_bau_electricity, rail_bau_electricity, industry_bau_electricity, emissions,
    freight, employees, gva, ch_gva, nuts_to_regions, industry_activity_codes,
    out_path, scaling_factors, industry_config
):

    nuts_to_regions_df = pd.read_csv(nuts_to_regions)
    units = gpd.read_file(units)
    nuts_dfs = {}
    for i in [2006, 2010, 2016]:
        nuts_dfs[i] = (
            nuts_to_regions_df
            .dropna(subset=[f'NUTS3_{i}'])
            .assign(country_code=nuts_to_regions_df.country_code.apply(util.get_alpha3))
            .set_index(f'NUTS3_{i}')
            [['id', 'country_code']]
            .rename(columns={'id': 'region'})
        )
    road_distance_df = util.read_tdf(road_distance)
    road_vehicles_df = util.read_tdf(road_vehicles)
    rail_demand_df = util.read_tdf(rail_demand)
    road_bau_electricity_df = util.read_tdf(road_bau_electricity)
    rail_bau_electricity_df = util.read_tdf(rail_bau_electricity)

    building_demand_df = util.read_tdf(building_demand)

    building_heat_electricity_consumption_df = util.read_tdf(building_heat_electricity_consumption)

    pop_weighted_df = subnational_pop_weighted_demand(
        units, building_demand_df, building_heat_electricity_consumption_df,
        road_distance_df, road_vehicles_df, rail_demand_df, road_bau_electricity_df, rail_bau_electricity_df,
        population
    )
    commercial_df = subnational_commercial_demand(
        building_demand_df, building_heat_electricity_consumption_df,
        road_distance_df, road_vehicles_df,
        road_bau_electricity_df, units, nuts_dfs[2016], gva, ch_gva
    )
    industry_df = subnational_industry_demand(
        units, industry_demand, road_distance_df, road_vehicles_df, rail_demand_df,
        air_demand, marine_demand, emissions, freight,
        rail_bau_electricity_df, industry_bau_electricity,
        employees, nuts_dfs[2006], industry_activity_codes, industry_config
    )

    all_df = pd.concat([commercial_df, industry_df, pop_weighted_df]).droplevel("country_code")

    # FillNA for any missing info
    _df = all_df.unstack('id').stack(dropna=False)
    print("Filling missing data for {}".format(_df[_df.isna()].sum(level=['dataset', 'cat_name', 'id']).index.remove_unused_levels()))
    all_df = _df.fillna(0).groupby(level=['dataset', 'cat_name', 'year', 'id', 'unit', 'end_use']).sum()

    # Scale
    all_df.loc[all_df.index.get_level_values('unit') == 'twh'] *= scaling_factors['energy']
    all_df.loc[all_df.index.get_level_values('unit') == 't'] *= scaling_factors['co2']
    all_df.loc[all_df.index.get_level_values('unit') == 'mio_km'] *= scaling_factors['transport']
    all_df = all_df.rename(
        {'twh': '{:.2f} twh'.format(1 / scaling_factors['energy']),
         't': '{:.2f} t'.format(1 / scaling_factors['co2']),
         'mio_km': '{:.2f} mio_km'.format(1 / scaling_factors['transport'])},
        level='unit'
    )
    all_df.to_csv(out_path)


def get_population_intensity(path_to_population, units):
    population_df = (
        pd.read_csv(path_to_population, index_col=0)
        .set_index(units.set_index(['id', 'country_code']).index)
    )
    return (
        population_df.div(population_df.sum(level='country_code')).population_sum
    )


def subnational_pop_weighted_demand(
    units, heat_demand, heat_electricity_consumption, road_distance_df, road_vehicles_df,
    rail_demand_df, road_bau_electricity_df, rail_bau_electricity_df, population
):
    concat_dfs = []
    population_intensity = get_population_intensity(population, units)

    concat_dfs.append(align_and_scale(
        heat_demand.xs('household', level='cat_name'), population_intensity, units
    ))

    concat_dfs.append(align_and_scale(
        heat_electricity_consumption.xs('household', level='cat_name'),
        population_intensity, units
    ))
    new_names = {
        'Motor coaches, buses and trolley buses': 'bus',
        'Passenger cars': 'passenger_car', 'Powered 2-wheelers': 'motorcycle'
    }
    for df in [road_distance_df, road_bau_electricity_df, road_vehicles_df]:
        passenger_df = align_and_scale(
            df.rename(util.get_alpha3, level='country_code')
            .drop(['Heavy duty vehicles', 'Light duty vehicles'], level=0),
            population_intensity, units
        )
        passenger_df.index = passenger_df.index.set_names('end_use', level='vehicle_type')
        passenger_df = passenger_df.rename(new_names, level='end_use')
        if 'twh' in passenger_df.index.get_level_values('unit'):  # bau electricity consumption
            passenger_df = passenger_df.rename(lambda x: f'{x}_bau_electricity', level='end_use')
        concat_dfs.append(passenger_df)

    passenger_rail_energy = align_and_scale(
        rail_demand_df.rename(util.get_alpha3, level='country_code')
        .drop('Freight', level=0).sum(level=['country_code', 'year', 'unit']),
        population_intensity, units
    )
    concat_dfs.append(
        passenger_rail_energy.to_frame('electricity').rename_axis(columns='end_use').stack()
    )
    passenger_rail_bau_electricity = align_and_scale(
        rail_bau_electricity_df.rename(util.get_alpha3, level='country_code')
        .xs('Passenger transport').sum(level=['country_code', 'year', 'unit']),
        population_intensity, units
    )
    concat_dfs.append(
        passenger_rail_bau_electricity.to_frame('bau_electricity').rename_axis(columns='end_use').stack()
    )

    household_scaled_df = pd.concat(
        [i.reorder_levels(LEVEL_ORDER) for i in concat_dfs],
        names=['dataset', 'cat_name'],
        keys=[('heat_demand', 'household'), ('heat_demand', 'household'),
              ('transport_demand', 'road'), ('transport_demand', 'road'),
              ('transport_vehicles', 'road'),
              ('transport_demand', 'rail'), ('transport_demand', 'rail')]
    )
    return household_scaled_df


def subnational_commercial_demand(
    heat_demand, heat_electricity_consumption, road_distance, road_vehicles,
    road_bau_electricity, units, nuts_2016, gva, ch_gva
):
    """
    Commercial sector demand is subdivided based on gross value added of economic
    activities which can be attributed to the commercial sector.
    """
    gva_df = util.read_eurostat_tsv(gva, ['unit', 'cat_name', 'region'])

    gva_eu_commercial = (
        gva_df
        .xs('MIO_EUR')
        .loc[['G-J', 'K-N', 'O-U']]  # commercial subsectors
        .sum(level='region', min_count=1)
        .reindex(nuts_2016.index)
    )
    gva_eu_total = (
        gva_df
        .xs('MIO_EUR')
        .loc['TOTAL']
        .sum(level='region', min_count=1)
        .reindex(nuts_2016.index)
    )
    gva_eu = (
        gva_eu_commercial
        .fillna(gva_eu_total)
        .rename_axis(index='nuts3')
        .set_index([nuts_2016.region, nuts_2016.country_code], append=True)
    )

    # Get GVA for Switzerland
    ch_gva_dfs = []
    for canton in pd.ExcelFile(ch_gva).sheet_names:
        ch_gva_df = pd.read_excel(
            ch_gva, sheet_name=canton, index_col=0, skiprows=3, skipfooter=16,
            usecols='A,C:L'
        )
        ch_gva_dfs.append(ch_gva_df.loc[['GHIJ', 'K', 'LMNRS', 'O', 'T']].sum())
    ch_gva_df = pd.concat(
        ch_gva_dfs, names=['canton', 'year'], keys=pd.ExcelFile(ch_gva).sheet_names
    ).unstack('year')
    ch_gva_df.index = ch_gva_df.index.map(CH_CANTON_TO_NUTS3)
    ch_gva_df.columns = util.to_numeric(ch_gva_df.columns).rename('year')
    gva_eu.loc[idx[ch_gva_df.index, :, :], ch_gva_df.columns] = ch_gva_df.values

    gva_intensity = (
        gva_eu
        .sum(level=['region', 'country_code'], min_count=1)
        .div(gva_eu.sum(level='country_code', min_count=1))
    )
    gva_intensity = gva_intensity.reindex(units.set_index(['id', 'country_code']).index)
    gva_intensity[units.set_index(['id', 'country_code'])['type'] == 'country'] = 1
    gva_intensity = (
        gva_intensity
        .fillna({col: gva_intensity.mean(axis=1) for col in gva_intensity.columns})
        .rename_axis(columns='year')
        .stack()
    )

    # scale to make up for missing Greek/Spanish islands
    gva_intensity = gva_intensity.div(gva_intensity.sum(level=['year', 'country_code']))
    commercial_heat_demand = align_and_scale(
        heat_demand.xs('commercial', level='cat_name'), gva_intensity, units
    )

    commercial_heat_electricity_consumption = align_and_scale(
        heat_electricity_consumption.xs('commercial', level='cat_name'),
        gva_intensity, units
    )

    ldv_road_distance = align_and_scale(
        road_distance.xs('Light duty vehicles')
        .rename(util.get_alpha3, level='country_code'),
        gva_intensity, units
    )
    ldv_road_distance = ldv_road_distance.to_frame('ldv').rename_axis(columns='end_use').stack()

    ldv_road_vehicles = align_and_scale(
        road_vehicles.xs('Light duty vehicles')
        .rename(util.get_alpha3, level='country_code'),
        gva_intensity, units
    )
    ldv_road_vehicles = ldv_road_vehicles.to_frame('ldv').rename_axis(columns='end_use').stack()

    ldv_road_bau_electricity = align_and_scale(
        road_bau_electricity.xs('Light duty vehicles')
        .rename(util.get_alpha3, level='country_code'),
        gva_intensity, units
    )

    ldv_road_bau_electricity = ldv_road_bau_electricity.to_frame('ldv_bau_electricity').rename_axis(columns='end_use').stack()

    commercial_scaled_df = pd.concat(
        [commercial_heat_demand.reorder_levels(LEVEL_ORDER),
         commercial_heat_electricity_consumption.reorder_levels(LEVEL_ORDER),
         ldv_road_distance.reorder_levels(LEVEL_ORDER),
         ldv_road_vehicles.reorder_levels(LEVEL_ORDER),
         ldv_road_bau_electricity.reorder_levels(LEVEL_ORDER)],
        names=['dataset', 'cat_name'],
        keys=[('heat_demand', 'commercial'), ('heat_demand', 'commercial'),
              ('transport_demand', 'road'), ('transport_vehicles', 'road'),
              ('transport_demand', 'road')]
    )

    return commercial_scaled_df


def industry_subsector_regional_intensity(
    units, emissions, freight, employees, nuts_2006, industry_activity_codes, industry_demand_df
):
    # Get freight data
    freight_df = util.read_eurostat_tsv(freight, ['subsector', 'unit', 'region'])

    nuts_2006.index = nuts_2006.index.str.replace('GR', 'EL')
    freight_eu = (
        freight_df
        .unstack()
        .groupby(FREIGHT_SECTORS, level=0).sum()
        .where(lambda x: x > 0)
        .loc[:, idx[:, nuts_2006.index]]
        .stack([0, 1])
    )
    freight_eu = (
        freight_eu
        .to_frame('freight')
        .set_index(freight_eu.index.get_level_values('region').str[:-1], append=True)
        .rename_axis(index=['subsector', 'year', 'nuts3', 'nuts2'])
    )

    # Get data on number of employees
    # V11210 == local units,  V16110 == persons employed, V13320 == wages and salaries
    industry_employees = util.read_eurostat_tsv(
        employees, ['cat_code', 'indicator', 'region'],
        slice_idx='V16110', slice_lvl='indicator'
    )

    activity_codes_df = pd.read_csv(
        industry_activity_codes, skipfooter=7, index_col=0, header=0, engine='python'
    ).dropna(subset=['Eurostat sector'])

    industry_employees = (
        industry_employees
        .loc[activity_codes_df['Eurostat sector'].dropna().index]
        .unstack()
        .groupby(activity_codes_df['Eurostat sector'].to_dict()).sum(min_count=1)
        .stack([0, 1])
        .rename_axis(index=['subsector', 'year', 'region'])
    )

    # Combine freight and employee data
    freight_employees = pd.concat([
        freight_eu['freight'].reset_index('nuts3'),
        industry_employees.reindex(freight_eu.droplevel('nuts3').index).to_frame('employees_nuts2')
    ], axis=1).set_index('nuts3', append=True)

    freight_employees = (
        freight_employees
        .set_index([
            freight_employees.index.get_level_values('nuts2').str[:2].rename('nuts0'),
            freight_employees.index.get_level_values('nuts2').str[:3].rename('nuts1'),
            nuts_2006.reindex(freight_employees.index.get_level_values('nuts3')).set_index('region').index
        ], append=True)
    )
    freight_employees['employees_nuts1'] = (
        industry_employees
        .reindex(freight_employees.droplevel(['nuts2', 'nuts3', 'nuts0']).index)
        .values
    )
    freight_employees['employees_nuts0'] = (
        industry_employees
        .reindex(freight_employees.droplevel(['nuts2', 'nuts3', 'nuts1']).index)
        .values
    )

    # Get contribution of eurospores regions to national industry activity in each subsector
    # We do this by going from NUTS2 employment data to NUTS3 by using freight in each subregion
    # Then going from NUTS3 up to Euro-Calliope regions
    # Where no employment data is available at any spatial resolution (e.g. in Switzerland), freight data is used directly
    def nuts3_intensity(lvl):
        return (
            freight_employees['freight']
            .div(freight_employees['freight'].sum(level=['subsector', 'year', lvl]))
        )

    employees_per_region = (
        freight_employees['employees_nuts2'].mul(nuts3_intensity('nuts2'))
        .fillna(freight_employees['employees_nuts1'].mul(nuts3_intensity('nuts1')))
        .fillna(freight_employees['employees_nuts0'].mul(nuts3_intensity('nuts0')))
        .sum(level=['subsector', 'year', 'region', 'nuts0'], min_count=1)
    )

    freight_per_region = (
        freight_employees['freight']
        .sum(level=['subsector', 'year', 'region', 'nuts0'], min_count=1)
    )
    # We use 2014 data, due to data gaps elsewhere
    subsector_employee_freight_intensity = (
        employees_per_region
        .div(employees_per_region.sum(level=['subsector', 'year', 'nuts0']))
        .fillna(freight_per_region.div(freight_per_region.sum(level=['subsector', 'year', 'nuts0'])))
        .unstack('subsector')
        .xs(2014)
        .droplevel('nuts0')
    )

    emissions_gdf = gpd.read_file(emissions)
    # Map industries to Euro-calliope regions
    emissions_eu = gpd.overlay(emissions_gdf, units)

    emissions_intensity_eu = (
        emissions_eu
        .groupby(['id', 'Subsector', 'country_code']).sum()
        .apply(lambda x: x / x.sum(level=['Subsector', 'country_code']))
        .reset_index('country_code', drop=True)
        .loc[:, 'emissions']
        .drop('Other Industrial Sectors', level='Subsector')
        .unstack('Subsector')
    )

    # Combine freight/employee subsector intensities with emissions subsector intensities
    all_intensities = (
        pd.concat([emissions_intensity_eu, subsector_employee_freight_intensity], axis=1)
        .reindex(units.set_index('id').index)
        .where(units.set_index('id')['type'] != 'country', other=1)
        .rename_axis(columns='subsector')
        .reindex(industry_demand_df.index.levels[0], axis=1)
        .set_index(units.set_index('country_code').index.map(util.get_alpha2), append=True)
    )
    # Find subsectors with zero intensity across all subregions in the country
    missing_intensities = (
        all_intensities.sum(level='country_code', min_count=1).isna().stack().where(lambda x: x).dropna()
    )
    # we normalise the average intensity, since it sometimes doesn't quite
    # add up to 1 in a country.
    average_intensity = (
        all_intensities.mean(axis=1).div(all_intensities.mean(axis=1).sum(level='country_code'))
    )
    for i in missing_intensities.index:
        print(
            "Missing industry regional subsector intensities for {}; "
            "filling with average from other subsectors as {}"
            .format(i, average_intensity.loc[idx[:, i[0]]])
        )
        all_intensities.loc[idx[:, i[0]], i[1]] = average_intensity.loc[idx[:, [i[0]]]]

    assert (all_intensities.sum(level='country_code').round(5) == 1).all().all()
    # anything else (e.g. one region in a country with NaN intensity): fill with zero
    return all_intensities.fillna(0).stack()


def subnational_industry_demand(
    units, industry_demand, road_distance_df, road_vehicles_df, rail_demand_df,
    air_demand, marine_demand, emissions, freight,
    rail_bau_electricity_df, industry_bau_electricity,
    employees, nuts_2006, industry_activity_codes, industry_config
):
    """
    Some industry subsectors are subdivided based on employment data per NUTS3 region.
    Where this isn't available, we use freight loading data. Both datasets have
    subcategories relating to industry subsectors.
    The more emitting subsectors (e.g. iron & steel) are subdivided using EU-ETS emissions
    data.
    """

    industry_demand_df = util.read_tdf(industry_demand)

    all_subsectors_industry_intensity = industry_subsector_regional_intensity(
        units, emissions, freight, employees, nuts_2006, industry_activity_codes, industry_demand_df
    )
    # FIXME: move this yearrange into config yaml
    industry_demand_df = industry_demand_df[
        industry_demand_df.index.get_level_values('year').isin(YEAR_RANGE)

    ]
    industry_energy_df = align_and_scale(
        industry_demand_df, all_subsectors_industry_intensity, units
    )
    industry_energy_df = (
        industry_energy_df
        .sum(level=['id', 'country_code', 'unit', 'year', 'carrier'])
        .rename_axis(index={'carrier': 'end_use'})
    )

    # BAU electricity
    industry_bau_electricity_df = util.read_tdf(industry_bau_electricity)

    industry_bau_electricity_df = industry_bau_electricity_df[
        industry_bau_electricity_df.index.get_level_values('year').isin(YEAR_RANGE)
    ]
    industry_bau_electricity_df = align_and_scale(
        industry_bau_electricity_df, all_subsectors_industry_intensity, units
    )
    industry_bau_electricity_df = (
        industry_bau_electricity_df
        .sum(level=['id', 'country_code', 'unit', 'year'])
        .to_frame('bau_electricity')
        .rename_axis(columns='end_use')
        .stack()
    )

    # Next three datasets can't handle industry subsectors, so we create a new
    # intensity df, based on all subsector energy consumption
    industry_energy_intensity = (
        industry_energy_df
        .xs('twh', level='unit')
        .sum(level=['id', 'country_code', 'year'])
        .div(industry_energy_df.xs('twh', level='unit').sum(level=['country_code', 'year']))
    )
    road_distance = align_and_scale(
        road_distance_df.xs('Heavy duty vehicles'), industry_energy_intensity, units
    )
    road_distance = road_distance.to_frame('hdv').rename_axis(columns='end_use').stack()
    road_vehicles = align_and_scale(
        road_vehicles_df.xs('Heavy duty vehicles'), industry_energy_intensity, units
    )
    road_vehicles = road_vehicles.to_frame('hdv').rename_axis(columns='end_use').stack()

    rail_freight = align_and_scale(
        rail_demand_df.xs('Freight'), industry_energy_intensity, units
    )
    rail_freight = rail_freight.to_frame('electricity').rename_axis(columns='end_use').stack()

    # BAU electricity
    rail_freight_bau = align_and_scale(
        rail_bau_electricity_df.xs('Freight transport'), industry_energy_intensity, units
    )
    rail_freight_bau = rail_freight_bau.to_frame('bau_electricity').rename_axis(columns='end_use').stack()

    # Air and Marine energy demand is also distributed according to total industry
    air_energy = align_and_scale(
        util.read_tdf(air_demand)
        .loc[idx[:, [i for i in YEAR_RANGE], :]]
        .rename_axis(index={'country': 'country_code'}),
        industry_energy_intensity,
        units
    ).to_frame('kerosene').rename_axis(columns='end_use').stack()
    marine_energy = align_and_scale(
        util.read_tdf(marine_demand)
        .loc[idx[:, [i for i in YEAR_RANGE], :]]
        .rename_axis(index={'country': 'country_code'}),
        industry_energy_intensity,
        units
    ).to_frame('diesel').rename_axis(columns='end_use').stack()

    industry_scaled_df = pd.concat(
        [industry_energy_df.reorder_levels(LEVEL_ORDER),
         industry_bau_electricity_df.reorder_levels(LEVEL_ORDER),
         road_distance.reorder_levels(LEVEL_ORDER),
         road_vehicles.reorder_levels(LEVEL_ORDER),
         rail_freight.reorder_levels(LEVEL_ORDER),
         rail_freight_bau.reorder_levels(LEVEL_ORDER),
         air_energy.reorder_levels(LEVEL_ORDER),
         marine_energy.reorder_levels(LEVEL_ORDER)],
        names=['dataset', 'cat_name'],
        keys=[('industry_demand', 'industry'), ('industry_demand', 'industry'),
              ('transport_demand', 'road'), ('transport_vehicles', 'road'),
              ('industry_demand', 'rail'),
              ('industry_demand', 'rail'), ('industry_demand', 'air'),
              ('industry_demand', 'marine')]
    )
    industry_scaled_df = electrify_industry_demand(
        industry_scaled_df, industry_config, level_order=LEVEL_ORDER
    )

    return industry_scaled_df


def align_and_scale(orig_df, scaling_df, units):

    aligned_df = scaling_df.align(orig_df)
    scaled_df = aligned_df[0].mul(aligned_df[1]).dropna()

    # make sure numbers add up
    check_scaling_df = scaled_df.sum(level=orig_df.index.names)

    assert np.allclose(
        check_scaling_df.sum(level='country_code'),
        orig_df.reindex(check_scaling_df.index).sum(level='country_code'),
        equal_nan=True
    )

    return scaled_df


if __name__ == "__main__":
    subnationalise_demand(
        population=snakemake.input.population,
        units=snakemake.input.units,
        building_demand=snakemake.input.annual_heat_demand,
        building_heat_electricity_consumption=snakemake.input.annual_heat_electricity_consumption,
        industry_demand=snakemake.input.industry_demand,
        road_distance=snakemake.input.road_distance,
        road_vehicles=snakemake.input.road_vehicles,
        rail_demand=snakemake.input.rail_demand,
        air_demand=snakemake.input.air_demand,
        marine_demand=snakemake.input.marine_demand,
        road_bau_electricity=snakemake.input.road_bau_electricity,
        rail_bau_electricity=snakemake.input.rail_bau_electricity,
        industry_bau_electricity=snakemake.input.industry_bau_electricity,
        emissions=snakemake.input.emissions,
        freight=snakemake.input.freight,
        employees=snakemake.input.employees,
        gva=snakemake.input.gva,
        ch_gva=snakemake.input.ch_gva,
        nuts_to_regions=snakemake.input.nuts_to_regions,
        industry_activity_codes=snakemake.input.industry_activity_codes,
        scaling_factors=snakemake.params.scaling_factors,
        industry_config=snakemake.params.industry_config,
        out_path=snakemake.output.all_annual
    )
