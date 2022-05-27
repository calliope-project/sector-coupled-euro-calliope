import pandas as pd
import numpy as np

import util

H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
RECYCLED_STEEL = 0.5  # 50% H-DRI Iron, 50% scrap steel
H2_TO_STEEL = (1 - RECYCLED_STEEL) * 0.05  # 0.05t_h2/t_steel in H-DRI
HDRI_CONSUMPTION = 0.0116  # H-DRI: 135kWh_e/t = 0.0116ktoe/kt
MWH_PER_T_TO_KTOE_PER_KT = 0.08598  # 1 MWh/t -> 1 GWh/kt -> 0.08598 ktoe/kt
METHANOL_LHV_KTOE = 0.476  # 19.915 MJ/kg LHV -> 19915000MJ/kt -> 0.476ktoe/kt

YEAR_RANGE = slice(2000, 2018)

final_projected_demand_level_order = ["subsector", "country_code", "unit", "carrier", "year"]
final_bau_demand_level_order = ["subsector", "country_code", "unit", "year"]


def get_industry_demand(
    path_to_energy_balances, path_to_cat_names, path_to_carrier_names,
    path_to_jrc_industry_end_use, path_to_jrc_industry_production,
    path_to_demand_scales, scale_demand,
    path_to_new_output, path_to_bau_output
):
    energy_df = pd.read_csv(path_to_jrc_industry_end_use, index_col=[0, 1, 2, 3, 4, 5, 6])
    prod_df = pd.read_csv(path_to_jrc_industry_production, index_col=[0, 1, 2, 3])
    energy_balances = pd.read_csv(path_to_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True)
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)

    demand = energy_df.xs('demand').sum(level=['section', 'subsection', 'country_code', 'cat_name', 'unit'])

    # if it can be met by electricity (exclusively or otherwise),
    # then it's an end-use electricity demand
    electrical_consumption = (
        get_carrier_demand('Electricity', demand, energy_df)
        .assign(carrier='electricity').set_index('carrier', append=True)
    )
    # If it can only be met by natural gas (steam heating) then it's natural gas
    nat_gas_consumption = (
        get_carrier_demand('Natural gas (incl. biogas)', demand, energy_df)
        .drop(electrical_consumption.droplevel('carrier').index, errors='ignore')
        .assign(carrier='methane').set_index('carrier', append=True)
    )
    # If it can only be met by diesel (backup generators) then it's diesel
    diesel_consumption = (
        get_carrier_demand('Diesel oil (incl. biofuels)', demand, energy_df)
        .drop(nat_gas_consumption.droplevel('carrier').index, errors='ignore')
        .drop(electrical_consumption.droplevel('carrier').index, errors='ignore')
        .assign(carrier='diesel').set_index('carrier', append=True)
    )

    steel_energy_consumption = get_steel_energy_consumption(energy_df, prod_df)
    chem_energy_consumption = get_chem_energy_consumption(electrical_consumption, prod_df, demand)

    # All electricity/nat gas consumption elsewhere
    all_other_industry_consumption = pd.concat([
        i
        .drop(['Iron and steel', 'Chemicals Industry'], level='cat_name', errors='ignore')
        .drop('Low enthalpy heat', level='subsection', errors='ignore')
        .sum(level=['cat_name', 'country_code', 'unit', 'carrier'])
        for i in [electrical_consumption, nat_gas_consumption, diesel_consumption]
    ])
    all_other_industry_space_heat = (
        demand
        .drop(['Iron and steel', 'Chemicals Industry'], level='cat_name')
        .xs('Low enthalpy heat', level='subsection')
        .sum(level=['cat_name', 'country_code', 'unit'])
        .assign(carrier='space_heat').set_index('carrier', append=True)
    )

    all_consumption = pd.concat([
        steel_energy_consumption.reorder_levels(all_other_industry_consumption.index.names),
        chem_energy_consumption.reorder_levels(all_other_industry_consumption.index.names),
        all_other_industry_consumption,
        all_other_industry_space_heat
    ])
    all_consumption.columns = all_consumption.columns.astype(int).rename('year')
    all_filled_consumption = fill_missing_data(energy_balances, cat_names, carrier_names, all_consumption)
    units = all_filled_consumption.index.get_level_values('unit')
    all_filled_consumption.loc[units == 'ktoe'] = (
        all_filled_consumption.loc[units == 'ktoe'].apply(util.ktoe_to_twh)
    )
    all_filled_consumption = all_filled_consumption.rename({'ktoe': 'twh'}, level='unit')
    all_filled_consumption.index = all_filled_consumption.index.set_names('subsector', level='cat_name')
    all_filled_consumption = all_filled_consumption.stack()

    if scale_demand:
        print("Scaling demands according to increase in value added from industry subsectors")
        demand_scales = (
            util.read_tdf(path_to_demand_scales)
            .xs(2050, level="year")
            .rename_axis(index={"id": "country_code"})
            .rename(lambda x: util.get_alpha2(x), level="country_code")
        )
        demand_scales = demand_scales.reindex(
            all_filled_consumption.groupby(level=demand_scales.index.names).sum().index
        ).fillna(1)

        all_filled_consumption = all_filled_consumption.mul(demand_scales)
        assert not all_filled_consumption.isna().any()

    all_filled_consumption.reorder_levels(final_projected_demand_level_order).to_csv(path_to_new_output)

    # Also save current electricity consumption, for removal from the ENTSOE hourly electricity profiles
    electricity_bau = (
        energy_balances
        .xs('E7000', level='carrier_code')
        .unstack(['unit', 'country', 'year'])
        .groupby(cat_names.jrc_idees.dropna().to_dict()).sum(min_count=1)
        .rename_axis(index='subsector', columns={'country': 'country_code'})
        .stack(['unit', 'country_code'])
        .bfill(axis=1)  # missing old data (BIH and MNE) backfilled
        .stack()
        .apply(util.tj_to_twh)
        .rename({'TJ': 'twh'}, level='unit')
    )
    electricity_bau.reorder_levels(final_bau_demand_level_order).to_csv(path_to_bau_output)


def get_carrier_demand(carrier, all_demand, energy_df):
    """
    Get demand for a specific carrier, assuming all end use demand that could consume
    that carrier are completely met by that carrier
    """
    energy = energy_df.xs(carrier, level='carrier_name')
    energy_efficiency = energy.xs('demand').div(energy.xs('consumption'))
    # Fill NaNs (where there is demand, but no electrical consumption in that country)
    # with the average efficiency a. from the country, b. from all countries
    energy_efficiency = (
        energy_efficiency
        .T  # have to reorient the array thanks to a NotImplementedError
        .fillna(energy_efficiency.mean(axis=1))
        .T  # have to reorient the array thanks to a NotImplementedError
        .fillna(energy_efficiency.mean())
    )
    return all_demand.reindex(energy_efficiency.index).div(energy_efficiency)


def get_steel_energy_consumption(energy_df, prod_df):
    """
    Calculates energy consumption in the iron and steel industry based on expected
    change in process to avoid fossil feedstocks. All process specific energy consumption
    (energy/t_steel) is based on the Electric Arc process (EAF), except sintering, which
    will be required for iron ore processed using H-DRI, but is not required by EAF.

    This function does the following:
    1. Finds all the specific consumption values by getting
        a. process energy demand / produced steel => specific demand
        b. process electrical demand / electrical consumption => electrical efficiency
        c. specific demand / electricial efficiency => specific electricity consumption
    2. Gets total process specific electricity consumption by adding specific consumptions
    for direct electric processes, EAF, H-DRI, smelting, sintering, refining, and finishing
    3. Gets specific hydrogen consumption for all countries that will process iron ore
    4. Gets specific space heat demand based on demand associated with EAF plants
    5. Gets total demand for electricity, hydrogen, and space heat by multiplying specific
    demand by total steel production (by both EAF and BF-BOF routes).
    """
    def get_specific_electricity_consumption(process, subprocess, energy_df, prod_df):
        consumption = energy_df.xs(('consumption', process, subprocess))
        demand = energy_df.xs(('demand', process, subprocess))
        specific_demand = demand.sum(level=['country_code']).div(prod_df.loc[process].droplevel('unit'))
        efficiency = demand.div(consumption)
        electrical_efficiency = (
            efficiency
            .where(efficiency > 0)
            .xs('Electricity', level='carrier_name')
            .T  # have to reorient the array thanks to a NotImplementedError
            .fillna(efficiency.xs('Electricity', level='carrier_name').mean(axis=1))
            .T  # have to reorient the array thanks to a NotImplementedError
            .fillna(efficiency.xs('Electricity', level='carrier_name').mean())
        )

        specific_consumption = specific_demand.div(electrical_efficiency).rename(index={'ktoe': 'ktoe/kt'})
        assert (
            specific_consumption.fillna(-1).droplevel('unit')
            >= specific_demand.fillna(-1)
        ).all().all()

        return specific_consumption.fillna(0)

    def get_auxiliary_electricity_consumption(process, energy_df, prod_df):
        auxiliaries = ['Lighting', 'Air compressors', 'Motor drives', 'Fans and pumps']
        consumption = (
            energy_df
            .xs(('consumption', process))
            .loc[auxiliaries]
            .sum(level=['country_code', 'unit', 'cat_name'])
        )
        specific_consumption = consumption.div(prod_df.loc[process].droplevel('unit'))
        specific_consumption.index = specific_consumption.index.set_levels(['ktoe/kt'], level='unit')
        return specific_consumption.fillna(0)

    # sintering/pelletising
    sintering_specific_consumption = get_specific_electricity_consumption(
        'Integrated steelworks', 'Steel: Sinter/Pellet making',
        energy_df, prod_df
    )

    # smelters
    eaf_smelting_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Smelters',
        energy_df, prod_df
    )

    # EAF
    eaf_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Electric arc',
        energy_df, prod_df
    )

    # Rolling & refining
    refining_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Furnaces, Refining and Rolling',
        energy_df, prod_df
    )
    finishing_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Products finishing',
        energy_df, prod_df
    )

    # Auxiliaries (lighting, motors, etc.)
    auxiliary_specific_consumption = get_auxiliary_electricity_consumption(
        'Electric arc', energy_df, prod_df
    )

    # Total electricity consumption
    ## If the country produces steel from Iron ore (assuming 50% recycling):
    ### sintering/pelletising * iron_ore_% + smelting * recycled_steel_% + H-DRI + EAF + refining/rolling + finishing + auxiliaries
    ## If the country only recycles steel:
    ### smelting + EAF + refining/rolling + finishing + auxiliaries

    total_specific_consumption = (
        sintering_specific_consumption.mul(1 - RECYCLED_STEEL).add(
            eaf_smelting_specific_consumption
            # if no sintering, this country/year recycles 100% of steel
            .where(sintering_specific_consumption == 0)
            # if there is some sintering, we update smelting consumption to equal our assumed 2050 recycling rate and add weighted H-DRI consumption to process the remaining iron ore
            .fillna(eaf_smelting_specific_consumption.mul(RECYCLED_STEEL).add(HDRI_CONSUMPTION))
        )
        .add(eaf_specific_consumption)
        .add(refining_specific_consumption)
        .add(finishing_specific_consumption)
        .add(auxiliary_specific_consumption)
    )
    # In case our model now says a country does produce steel,
    # we give them the average of energy consumption of all other countries
    total_specific_consumption = (
        total_specific_consumption
        .where(total_specific_consumption > 0)
        .fillna(total_specific_consumption.mean())
        .assign(carrier='electricity')
        .set_index('carrier', append=True)
    )

    # Hydrogen consumption for H-DRI, only for those country/year combinations that handle iron ore
    # and don't recycle all their steel
    h2_specific_consumption = H2_LHV_KTOE * H2_TO_STEEL
    total_specific_h2_consumption = (
        total_specific_consumption
        .where(sintering_specific_consumption > 0)
        .fillna(0)
        .where(lambda x: x == 0)
        .fillna(h2_specific_consumption)
        .rename(index={'electricity': 'hydrogen'})
    )
    total_specific_consumption = total_specific_consumption.append(total_specific_h2_consumption)
    # Space heat
    space_heat_specific_demand = (
        energy_df
        .xs(('demand', 'Electric arc', 'Low enthalpy heat'))
        .div(prod_df.xs('Electric arc').droplevel('unit'))
        .assign(carrier='space_heat').set_index('carrier', append=True)
        .sum(level=total_specific_consumption.index.names)
        .rename(index={'ktoe': 'ktoe/kt'})
    )
    total_specific_consumption = total_specific_consumption.append(space_heat_specific_demand)

    steel_consumption = (
        total_specific_consumption
        .mul(prod_df.xs('Iron and steel', level='cat_name').sum(level='country_code'), level='country_code')
        .rename(index={'ktoe/kt': 'ktoe'})
    )

    return steel_consumption


def get_chem_energy_consumption(electrical_consumption, prod_df, demand):
    """
    We remove feedstock and steam-based processing from "basic chemicals" production,
    which will be replaced by direct H2 provision.
    All other demand (non-basic chemicals & basic chemicals electricity) is directly passed over.
    """
    chem_electricity_consumption = (
        electrical_consumption
        .drop('Low enthalpy heat', level='subsection')
        .xs('Chemicals Industry', level='cat_name')
        .sum(level=['country_code', 'unit', 'carrier'])
    )

    h2_demand = {  # t/t
        'Ethylene': 0,
        'Propylene': 0,
        'BTX': 0,
        'Ammonia': 0.178,
        'Methanol': 0.189,
    }
    methanol_demand = {
        'Ethylene': 2.83,
        'Propylene': 2.83,
        'BTX': 4.3,
        'Ammonia': 0,
        'Methanol': 1
    }
    co2_demand = {  # tCO2/t
        'Ethylene': 0,
        'Propylene': 0,
        'BTX': 0,
        'Ammonia': 0.112,  # inc. 0.35t_urea/t_ammonia
        'Methanol': 1.373,
    }

    energy_demand = {  # MWh/t
        'Ethylene': 1.4,  # MTO process
        'Propylene': 1.4,  # MTO process
        'BTX': 1.4,  # MTA process
        'Ammonia': 2.05,  # inc. 0.35t_urea/t_ammonia
        'Methanol': 1.5,  # auxiliary demand, could just be assumed as already included in existing electricity demand
    }

    mass = {  # Bazzanella and Ausfelder, 2017
        'Ethylene': 21.7,
        'Propylene': 17,
        'BTX': 15.7,
        'Ammonia': 17,
        'Methanol': 2,
    }
    molar_mass = {
        'Ethylene': 28.05,
        'Propylene': 42.08,
        'BTX': 93,
        'Ammonia': 17.01,
        'Methanol': 32.04,
    }
    moles = {k: mass[k] / molar_mass[k] for k in mass.keys()}
    molar_ratio = {k: v / sum(moles.values()) for k, v in moles.items()}

    basic_chemicals_mass = prod_df.xs('Basic chemicals (kt ethylene eq.)').sum(level=['country_code'])
    basic_chemicals_moles = basic_chemicals_mass / molar_mass['Ethylene']
    masses = {
        k: basic_chemicals_moles * molar_ratio[k] * molar_mass[k]
        for k in molar_ratio.keys()
    }

    chem_h2_consumption = sum(masses[k] * h2_demand[k] * H2_LHV_KTOE for k in h2_demand.keys())
    chem_methanol_consumption = sum(masses[k] * methanol_demand[k] * METHANOL_LHV_KTOE for k in methanol_demand.keys())
    chem_co2_consumption = sum(masses[k] * co2_demand[k] for k in co2_demand.keys())
    chem_energy_consumption = sum(masses[k] * energy_demand[k] * MWH_PER_T_TO_KTOE_PER_KT for k in energy_demand.keys())
    # Space heat
    space_heat_demand = (
        demand
        .xs('Low enthalpy heat', level='subsection')
        .xs('Chemicals Industry', level='cat_name')
        .sum(level=['country_code', 'unit'])
    )

    chem_consumption = pd.concat(
        [chem_electricity_consumption.add(chem_energy_consumption).droplevel('carrier'),
         chem_methanol_consumption.assign(unit='ktoe').set_index('unit', append=True),
         chem_h2_consumption.assign(unit='ktoe').set_index('unit', append=True),
         chem_co2_consumption.mul(1e3).assign(unit='t').set_index('unit', append=True),
         space_heat_demand],

        names=['carrier'], keys=['electricity', 'methanol', 'hydrogen', 'co2', 'space_heat']
    )

    chem_consumption = chem_consumption.assign(cat_name='Chemicals Industry').set_index('cat_name', append=True)

    return chem_consumption


def fill_missing_data(energy_balances, cat_names, carrier_names, energy_consumption):
    """
    There are 7 countries without relevant data in JRC_IDEES, so we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    2016-2018 data for all 35 countries is filled in based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    # Get annual energy balances
    industry_energy_balances = (
        energy_balances
        .unstack(['year', 'country'])
        .groupby([
            cat_names.jrc_idees.dropna().to_dict(),
            carrier_names.ind_carrier_name.dropna().to_dict()
        ], level=['cat_code', 'carrier_code']).sum(min_count=1)
        .sum(level='cat_code', min_count=1)
        .stack('country')
        .rename_axis(index=['cat_name', 'country_code'])
        .apply(util.tj_to_ktoe)
    )

    # If JRC-IDEES data exists, there will be data in 'energy_consumption' for that country
    balances_with_jrc_data = industry_energy_balances.loc[
        pd.IndexSlice[:, energy_consumption.index.levels[1]], energy_consumption.columns
    ]
    # If JRC-IDEES data does not exist, there will only be data for that country in annual energy balances
    balances_without_jrc_data = industry_energy_balances.drop(
        energy_consumption.index.levels[1], level='country_code'
    )
    # Compared to countries with JRC data, those without JRC-IDEES data consume X amount of energy
    # E.g. CH has ~1% of paper and pulp energy consumption compared to all countries with JRC data
    countries_without_jrc_data_contribution = balances_without_jrc_data.div(
        balances_with_jrc_data.sum(level='cat_name')
    )
    # Multiply the JRC data with the contribution from each country
    # So all JRC countries consume 3.43e4ktoe electricity in paper and pulp in 2014,
    # hence CH consumes ~3.43e2ktoe electricity in paper and pulp in 2014
    countries_without_jrc_data_consumption = (
        energy_consumption
        .sum(level=['cat_name', 'unit', 'carrier'])
        .align(countries_without_jrc_data_contribution)[0]
        .mul(energy_consumption
             .sum(level=['cat_name', 'unit', 'carrier'])
             .align(countries_without_jrc_data_contribution)[1])
    )
    all_euro_calliope_consumption = energy_consumption.append(
        countries_without_jrc_data_consumption
        .reorder_levels(energy_consumption.index.names)
    )
    average_consumption_per_energy_use = (
        all_euro_calliope_consumption
        .div(industry_energy_balances)
        .mean(axis=1)
        .where(lambda x: ~np.isinf(x) & (x > 0))
    )
    # In some instances (e.g. RO: non ferrous metals), JRC IDEES has consumption, but
    # latest Eurostat doesn't, leading to inf values
    average_consumption_per_energy_use = (
        average_consumption_per_energy_use
        .unstack(['cat_name', 'carrier', 'unit'])
        .fillna(
            average_consumption_per_energy_use
            .mean(level=['cat_name', 'carrier', 'unit'])
        )
        .unstack()
        .reorder_levels(all_euro_calliope_consumption.index.names)
    )
    # Fill data where JRC says there is no consumption of any form in a country's industry sybsector
    # But where the energy balances show consumption (e.g. UK, Wood and wood products)
    _to_fill = all_euro_calliope_consumption.stack().unstack(["carrier", "unit"])
    _filler = industry_energy_balances.stack().mul(average_consumption_per_energy_use, axis=0).unstack(["carrier", "unit"])
    filled_jrc_no_data = _to_fill.where(_to_fill.sum(axis=1) > 0).fillna(_filler)

    # Fill data for all years between 2016 and 2018, since there is no data from JRC
    _to_fill = filled_jrc_no_data.stack(["unit", "carrier"]).unstack("year")
    new_years_filler = _filler.stack(["unit", "carrier"]).unstack("year")
    filled_2016_to_2018 = _to_fill.assign(**{str(_year): new_years_filler.loc[_to_fill.sum(axis=1) > 0, _year] for _year in range(2016, 2019)})
    filled_2016_to_2018.columns = filled_2016_to_2018.columns.astype(int)

    # Fill any remaining missing data (e.g. BA, pre 2013) with backfill (older years) / forwardfill (newer years) / linear interpolation (middle years)
    all_filled = filled_2016_to_2018.interpolate(axis=1, limit_direction="both")

    verify_data(all_filled, industry_energy_balances)

    return all_filled


def verify_data(euro_calliope_consumption, industry_energy_balances):
    """
    Check that all our processing of data leads to sensible results relative to
    industry energy balances. What we mean by this is that > 95% of total euro-calliope
    energy demand should be within 0.5x and 1.5x historical primary energy
    consumption across subsectors. It's a crude estimate, but should catch cases
    where we're completely off / have missed some energy demand.
    """
    # Get the ratio of euro-calliope energy consumption to eurostat energy consumption
    consumption_diff = (
        euro_calliope_consumption
        .xs('ktoe', level='unit')
        .sum(level=['cat_name', 'country_code'], min_count=1)
        .div(industry_energy_balances.where(industry_energy_balances > 0))
        .loc[:, YEAR_RANGE]
        .stack()
    )
    # If industry energy balances have finite values, so does euro-calliope
    assert consumption_diff[consumption_diff == 0].sum() == 0

    # Calculate contribution of outliers in the consumption ratios
    Q1 = consumption_diff.quantile(0.25)
    Q3 = consumption_diff.quantile(0.75)
    IQR = Q3 - Q1
    # We would consider ratios outside 0.5-1.5 to be outliers, so just make sure the
    # threshold for outliers in the data is at least better than that
    assert Q1 - IQR * 3 >= 0.5
    assert Q3 + IQR * 3 <= 1.5
    outliers = consumption_diff[(consumption_diff < 0.5) | (consumption_diff > 1.5)]

    # We're happy with < 10% outliers
    assert len(outliers) / len(consumption_diff) <= 0.10

    # We're happy with < 5% of total demand being within the outliers
    assert (
        industry_energy_balances.stack().reindex(outliers.index).sum() /
        industry_energy_balances.loc[:, YEAR_RANGE].stack().sum() <= 0.05
    )


if __name__ == "__main__":
    get_industry_demand(
        path_to_energy_balances=snakemake.input.energy_balances,
        path_to_cat_names=snakemake.input.cat_names,
        path_to_carrier_names=snakemake.input.carrier_names,
        path_to_jrc_industry_end_use=snakemake.input.jrc_industry_end_use,
        path_to_jrc_industry_production=snakemake.input.jrc_industry_production,
        path_to_demand_scales=snakemake.input.demand_scales,
        scale_demand=snakemake.params.scale_demand,
        path_to_new_output=snakemake.output.new_demand,
        path_to_bau_output=snakemake.output.bau_electricity
    )
