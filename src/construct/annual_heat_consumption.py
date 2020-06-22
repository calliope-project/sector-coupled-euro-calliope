import glob

import pandas as pd
import numpy as np
import xarray as xr

import util

idx = pd.IndexSlice

END_USE_CAT_NAMES = {
    'FC_OTH_HH_E': 'all',
    'FC_OTH_HH_E_CK': 'cooking',
    'FC_OTH_HH_E_SH': 'space_heating',
    'FC_OTH_HH_E_WH': 'water_heating'
}

HH_END_USE_CARRIERS_ALL = {  # we keep these carriers, since we don't want fossils subcategories but we do want renewables subcategories
    'E7000': 'electricity',
    'G3000': 'natural_gas',
    'H8000': 'heat',
    'O4000': 'oil',
    'R5110-5150_W6000RI': 'biomass',
    'R5300': 'biogas',
    'RA410': 'solar_thermal',
    'RA600': 'ambient_heat',
    'SFF_P1000_S2000': 'coal',
    'TOTAL': 'total'
}
HH_END_USE_RENEWABLES = [  # These are the carriers that fall into the renewables subcategory
    'R5110-5150_W6000RI', 'R5300', 'RA410', 'RA600'
]

CH_ENERGY_CARRIER_TRANSLATION = {
    'Heizöl': 'oil',
    'Erdgas': 'natural_gas',
    'El. Widerstandsheizungen': 'direct_electric',
    'El. Wärmepumpen 1)': 'heat_pump',
    "El. Ohm'sche Anlagen": 'direct_electric',
    'El. Wärmepumpen': 'heat_pump',
    'Elektrizität': 'electricity',
    'Holz': 'biomass',
    'Kohle': 'coal',
    'Fernwärme': 'heat',
    'Umweltwärme': 'ambient_heat',
    'Solar': 'solar_thermal',
}

CH_HH_END_USE_TRANSLATION = {
    'Raumwärme': 'space_heating',
    'Warmwasser': 'water_heating',
    #'Klima, Lüftung, HT': 'Space cooling',
    #'Unterhaltung, I&K': 'Lighting and appliances',
    'Kochen / Geschirrspülen': 'cooking',
    #'Beleuchtung': 'Lighting and appliances',
    #'Waschen & Trocknen': 'Lighting and appliances',
    #'Gefrieren & Kühlen': 'Lighting and appliances',
    #'sonstige Elektrogeräte': 'Lighting and appliances',
}


def generate_annual_energy_balance_nc(
    hh_end_use, ch_end_use, countries, path_to_output
):
    """
    Combine annual energy balance information with data on household
    end use energy consumption
    """

    country_codes = [util.get_alpha2(i, eurostat=True) for i in countries]

    hh_end_use_df = pd.read_csv(hh_end_use, delimiter='\t', index_col=0)
    hh_end_use_df.index = (
        hh_end_use_df.index.str.split(',', expand=True)
        .rename(['cat_code', 'carrier_code', 'unit', 'country'])
    )
    hh_end_use_df.columns = hh_end_use_df.columns.astype(int).rename('year')
    hh_end_use_df = (
        hh_end_use_df
        .loc[idx[[i for i in END_USE_CAT_NAMES.keys()], :, 'TJ', country_codes], :]
        .transform(util.to_numeric)
        .astype(float)
        .dropna(how='all')
    )
    hh_end_use_df.index = (
        hh_end_use_df.index
        .remove_unused_levels()
        .set_levels([END_USE_CAT_NAMES.values(), HH_END_USE_CARRIERS_ALL.values()],
                    level=['cat_code', 'carrier_code'])
        .rename(['cat_name', 'carrier_name'], level=['cat_code', 'carrier_code'])
    )

    ch_hh_end_use_df_tot = get_ch_sheet(  # maybe remove - do we need total hh energy consumption info?
        ch_end_use, 'Tabelle 15', skipfooter=5, translation=CH_HH_END_USE_TRANSLATION
    )
    ch_hh_end_use_df_sh = get_ch_sheet(
        ch_end_use, 'Tabelle 18', skipfooter=8, translation=CH_ENERGY_CARRIER_TRANSLATION
    )
    ch_hh_end_use_df_hw = get_ch_sheet(
        ch_end_use, 'Tabelle 20', skipfooter=5, translation=CH_ENERGY_CARRIER_TRANSLATION
    )
    ch_hh_end_use_df_c = get_ch_sheet(
        ch_end_use, 'Tabelle 21', skipfooter=4, translation=CH_ENERGY_CARRIER_TRANSLATION
    )

    ch_hh_end_use_df = (
        pd.concat(
            [ch_hh_end_use_df_sh, ch_hh_end_use_df_hw, ch_hh_end_use_df_c],
            keys=('space_heating', 'water_heating', 'cooking'),
            names=['cat_name', 'carrier_name']
        )
        .append(
            ch_hh_end_use_df_tot  # maybe remove - do we need total hh energy consumption info?
            .assign(carrier_name='total')
            .set_index('carrier_name', append=True)
        )
        .assign(unit='TJ', country_code='CH')
        .set_index(['unit', 'country_code'], append=True)
        .drop('Δ ’00 – ’18', axis=1)
    )
    ch_hh_end_use_df.columns = ch_hh_end_use_df.columns.astype(int)

    hh_end_use_df = (
        hh_end_use_df
        .append(ch_hh_end_use_df, sort=True)
        .sort_index()
        .where(hh_end_use_df > 0)
        .dropna(how='all')
    )

    hh_end_use_df.to_csv(path_to_output)


def get_ch_sheet(path_to_excel, sheet, skipfooter, translation):
    _df = (
        pd.read_excel(
        path_to_excel, sheet_name=sheet, skiprows=9, skipfooter=skipfooter, index_col=1)
        .drop('Unnamed: 0', axis=1)
    )
    return _df.groupby(translation).sum()


if __name__ == "__main__":
    generate_override(
        hh_end_use=snakemake.input.hh_end_use
        ch_end_use=snakemake.input.ch_end_use
        countries=snakemake.params.countries
        path_to_output=snakemake.output
    )