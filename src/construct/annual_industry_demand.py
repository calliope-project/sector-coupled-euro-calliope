import pandas as pd
import numpy as np

import util

H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
RECYCLED_STEEL = 0.5  # 50% H-DRI Iron, 50% scrap steel
H2_TO_STEEL = RECYCLED_STEEL * 0.05  # 0.05t_h2/t_steel in H-DRI
HDRI_CONSUMPTION = 0.0116  # H-DRI: 135kWh/t = 0.0116ktoe/kt
MWH_PER_T_TO_KTOE_PER_KT = 0.08598  # 1 MWh/t -> 1 GWh/kt -> 0.08598 ktoe/kt
METHANOL_LHV_KTOE = 0.476  # 19.915 MJ/kg LHV -> 19915000MJ/kt -> 0.476ktoe/kt


def get_industry_demand(
    path_to_energy_balances, path_to_cat_names,
    path_to_jrc_industry_end_use, path_to_jrc_industry_production,
    path_to_new_output, path_to_bau_output
):
    energy_df = pd.read_csv(path_to_jrc_industry_end_use, index_col=[0, 1, 2, 3, 4, 5, 6])
    prod_df = pd.read_csv(path_to_jrc_industry_production, index_col=[0, 1, 2, 3])
    energy_balances = pd.read_csv(path_to_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True)
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)

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
    # If it can only be met by diesel (backup generators) then it's oil
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

    all_filled_consumption = fill_missing_data(energy_balances, cat_names, all_consumption)
    units = all_filled_consumption.index.get_level_values('unit')
    all_filled_consumption.loc[units == 'ktoe'] = (
        all_filled_consumption.loc[units == 'ktoe'].apply(util.ktoe_to_twh)
    )
    all_filled_consumption = all_filled_consumption.rename({'ktoe': 'twh'}, level='unit')
    all_filled_consumption.index = all_filled_consumption.index.set_names('subsector', level='cat_name')

    all_filled_consumption.stack().to_csv(path_to_new_output)

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
    electricity_bau.to_csv(path_to_bau_output)


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
    change in process to avoid fossil feedstocks
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

        specific_consumption = specific_demand.div(electrical_efficiency)
        specific_consumption.index = specific_consumption.index.set_levels(['ktoe/kt'], level='unit')
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

    ## sintering/pelletising
    sintering_specific_consumption = get_specific_electricity_consumption(
        'Integrated steelworks', 'Steel: Sinter/Pellet making',
        energy_df, prod_df
    )

    ## smelters
    smelting_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Smelters',
        energy_df, prod_df
    )

    ## EAF
    eaf_specific_consumption = get_specific_electricity_consumption(
        'Electric arc', 'Steel: Electric arc',
        energy_df, prod_df
    )

    ## Rolling & refining
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

    ## Total
    ### If the country produces steel from Iron ore:
    #### sintering/pelletising * 0.5 + smelting * 0.5 + H-DRI + EAF + refining/rolling + finishing + auxiliaries
    ### If the country only recycles steel:
    #### smelting + EAF + refining/rolling + finishing + auxiliaries

    total_specific_consumption = (
        sintering_specific_consumption.mul(1 - RECYCLED_STEEL).add(
            smelting_specific_consumption
            .where(sintering_specific_consumption == 0)
            .fillna(smelting_specific_consumption.mul(RECYCLED_STEEL).add(HDRI_CONSUMPTION))
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

    # Hydrogen consumption for H-DRI
    h2_specific_consumption = H2_LHV_KTOE * H2_TO_STEEL
    total_specific_h2_consumption = total_specific_consumption.where(lambda x: x < 0).fillna(h2_specific_consumption)
    total_specific_h2_consumption.index = total_specific_h2_consumption.index.set_levels(['hydrogen'], level='carrier')
    total_specific_consumption = total_specific_consumption.append(total_specific_h2_consumption)
    # Space heat
    space_heat_specific_demand = (
        energy_df
        .xs(('demand', 'Electric arc', 'Low enthalpy heat'))
        .div(prod_df.xs('Electric arc').droplevel('unit'))
        .assign(carrier='space_heat').set_index('carrier', append=True)
        .sum(level=total_specific_consumption.index.names)
    )
    space_heat_specific_demand.index = space_heat_specific_demand.index.set_levels(['ktoe/kt'], level='unit')
    total_specific_consumption = total_specific_consumption.append(space_heat_specific_demand)

    steel_consumption = total_specific_consumption.mul(
        prod_df.xs('Iron and steel', level='cat_name').sum(level='country_code'),
        level='country_code'
    )
    steel_consumption.index = steel_consumption.index.set_levels(['ktoe'], level='unit')

    return steel_consumption


def get_chem_energy_consumption(electrical_consumption, prod_df, demand):
    """
    We remove feedstock and steam-based processing, which will be replaced by direct H2 provision
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
         chem_co2_consumption.assign(unit='kt').set_index('unit', append=True),
         space_heat_demand],

        names=['carrier'], keys=['electricity', 'methanol', 'hydrogen', 'co2', 'space_heat']
    )

    chem_consumption = chem_consumption.assign(cat_name='Chemicals Industry').set_index('cat_name', append=True)

    return chem_consumption


def fill_missing_data(energy_balances, cat_names, energy_consumption):
    """
    There are 7 countries without relevant data in JRC_IDEES, so we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    2016-2018 data for all 35 countries is filled in based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    industry_energy_balances = (
        energy_balances
        .unstack(['year', 'country'])
        .xs('TOTAL', level='carrier_code')
        .groupby(cat_names.jrc_idees.dropna().to_dict(), level=0)
        .sum()
        .stack('country')
        .rename_axis(index=['cat_name', 'country_code'])
    )
    balances_with_jrc_data = industry_energy_balances.loc[
        pd.IndexSlice[:, energy_consumption.index.levels[1]], energy_consumption.columns
    ]
    balances_without_jrc_data = industry_energy_balances.drop(
        energy_consumption.index.levels[1], level='country_code'
    )
    countries_without_jrc_data_contribution = balances_without_jrc_data.div(
        balances_with_jrc_data.sum(level='cat_name')
    )
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
        .where(lambda x: ~np.isinf(x))
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

    filled_2016_to_2018 = all_euro_calliope_consumption.fillna(
        industry_energy_balances.mul(average_consumption_per_energy_use, axis=0)
    )
    all_filled = filled_2016_to_2018.T.fillna(filled_2016_to_2018.mean(axis=1)).T

    return all_filled


if __name__ == "__main__":
    get_industry_demand(
        path_to_energy_balances=snakemake.input.energy_balances,
        path_to_cat_names=snakemake.input.cat_names,
        path_to_jrc_industry_end_use=snakemake.input.jrc_industry_end_use,
        path_to_jrc_industry_production=snakemake.input.jrc_industry_production,
        path_to_new_output=snakemake.output.new_demand,
        path_to_bau_output=snakemake.output.bau_electricity
    )
