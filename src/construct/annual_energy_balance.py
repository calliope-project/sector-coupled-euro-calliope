import glob
from string import digits

import click
import pandas as pd
import numpy as np
import xarray as xr

import util

@click.command()
@click.argument("path_to_input")
@click.argument("path_to_cat_names")
@click.argument("path_to_carrier_names")
@click.argument("path_to_ch_excel")
@click.argument("path_to_output")
@click.argument("countries")
def generate_annual_energy_balance_nc(
    path_to_input, path_to_cat_names, path_to_carrier_names, path_to_ch_excel,
    path_to_result, countries
):
    """
    Open a TSV file and reprocess it into a xarray dataset, including long names for
    Eurostat codes.
    Switzerland is not included in Eurostat, so we splicei n data from their govt.
    statistics.
    """
    # Names for each consumption category/sub-category and carriers have been prepared by hand
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)
    country_codes = [util.get_alpha2(i, eurostat=True) for i in countries]

    df = pd.read_csv(path_to_input, delimiter='\t', index_col=0)
    df.index = (
        df.index.str.split(',', expand=True)
        .rename(['cat_code', 'carrier_code', 'unit', 'country'])  # comes as 'nrg_bal,siec,unit,geo\\time'
    )
    df.columns = df.columns.astype(int).rename('year')
    df = df.transform(util.to_numeric)
    df = df.reorder_levels([
        'cat_code', 'carrier_code', 'unit', 'country'
    ])
    df = df.loc[
        [cat_names.index, carrier_names.index, 'TJ', country_codes], :
    ].dropna(how='all')

    demand_ds = df.stack().unstack('country').to_xarray()
    demand_ds = demand_ds.assign_coords(cat_names.to_xarray())
    demand_ds = demand_ds.assign_coords(carrier_names.to_xarray())

    demand_ds.attrs['units'] = 'TJ'


    # Add CH energy use (only covers a subset of sectors and carriers, but should be enough)
    demand_ds = demand_ds.merge(ch_energy_use_ds)

    demand_ds.to_netcdf(path_to_result)


def add_ch_energy_balance(path_to_ch_excel):
    ch_hh_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, 'T17a', skipfooter=9).assign(cat_code='FC_OTH_HH_E'
    ).reset_index()
    ch_ind_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, 'T17b', skipfooter=12).assign(cat_code='FC_IND_E'
    ).reset_index()
    ch_ser_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, 'T17c', skipfooter=12).assign(cat_code='FC_OTH_CP_E'
    ).reset_index()

    ch_energy_use_ds = (
        pd.concat([ch_hh_energy_use, ch_ind_energy_use, ch_ser_energy_use])
        .set_index(['year', 'cat_code']).stack().to_frame(name='CH').to_xarray()
    ).transpose('cat_code', 'carrier_code', 'year')

    return ch_energy_use_ds


def get_ch_energy_balance_sheet(path_to_excel, sheet, skipfooter):
    ch_energy_carriers = {
        'Erdölprodukte': 'O4000XBIO',
        'Elektrizität': 'E7000',
        'Gas': 'G3000',
        'Kohle': 'C0000X0350-0370',
        'Holzenergie': 'R5110-5150_W6000RI',
        'Fernwärme': 'H8000',
        'Industrieabfälle': 'W6100_6220',
        'Übrige erneuerbare Energien': 'RA000',
        'Total\n= %': 'TOTAL',
    }
    # Footnote labels lead to some strings randomly ending in numbers; we remove them here
    remove_digits = str.maketrans('', '', digits)
    _df = (
        pd.read_excel(path_to_excel, skiprows=6, skipfooter=skipfooter,
                      index_col=0, sheet_name=sheet, header=[0, 1, 2, 3, 4])
        # Ignore columns giving % use
        .loc[slice(2010, None), idx[:, :, :, :, 'TJ']]
        # ignore the column giving subset of oil use which is light oil
        .iloc[:, [i for i in range(10) if i != 1]]
    )
    _df.columns = (
        _df.columns.get_level_values(0).str.translate(remove_digits)
        .map(ch_energy_carriers).rename('carrier_code')
    )
    _df.index.rename('year', inplace=True)

    _df = _df.transform(util.to_numeric).astype('float')

    return _df


if __name__ == "__main__":
    generate_annual_energy_balance_nc(
        path_to_input=snakemake.input.energy_balance
        path_to_ch_excel=snakemake.input.ch_energy_balance
        path_to_cat_names=snakemake.input.cat_names
        path_to_carrier_names=snakemake.input.carrier_names
        countries=snakemake.params.countries
        path_to_result=snakemake.output
    )