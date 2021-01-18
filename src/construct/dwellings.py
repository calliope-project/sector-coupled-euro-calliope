import pandas as pd

import util

idx = pd.IndexSlice


def get_dwelling_ratio(path_to_dwellings, regions, nuts_to_regions, out_path):
    """
    Per region, get the ratio of single family homes to multi-family homes.
    Although the underlying data is valid for (2010?) NUTS regions, we just use
    the national average applied to each region in a country.
    Missing locations (e.g. croatia) get the average across all other datasets
    """
    regions_df = pd.read_csv(regions, index_col=0)
    dwellings = pd.read_csv(path_to_dwellings, delimiter='\t', index_col=0)
    dwellings.index = dwellings.index.str.split(',', expand=True).rename(
        ['housing', 'building_type', 'unit', 'year']
    )
    dwellings.columns = dwellings.columns.str.strip().rename('id')

    dwellings = dwellings.droplevel(['unit', 'year']).apply(util.to_numeric)
    if "national" in out_path:
        grouper = regions_df.index.apply(util.get_alpha2)
    else:
        nuts_2010 = pd.read_csv(nuts_to_regions)
        nuts_2010['NUTS3_2010'] = nuts_2010["NUTS3_2010"].fillna(
            nuts_2010.country_code.where(nuts_2010.Source != 'NUTS3')
        )
        grouper = nuts_2010.dropna(subset=['NUTS3_2010']).set_index('NUTS3_2010').id

    # group dwelling data into regions
    regional_dwellings = dwellings.groupby(grouper, axis=1).sum()
    regional_dwellings_SFH = regional_dwellings.xs(("DW", "RES1"))
    regional_dwellings_MFH = regional_dwellings.loc[idx["DW", ["RES2", "RES_GE3"]], :].sum()

    MFH_to_SFH = regional_dwellings_MFH / (regional_dwellings_SFH + regional_dwellings_MFH)
    MFH_to_SFH = MFH_to_SFH.reindex(regions_df.index)

    # fill empties
    print(
        'Filling missing regions {} with average dwelling ratio'
        .format(MFH_to_SFH[MFH_to_SFH.isnull()].index.values)
    )
    MFH_to_SFH = MFH_to_SFH.fillna(MFH_to_SFH.mean())

    MFH_to_SFH.to_csv(out_path)


if __name__ == "__main__":
    get_dwelling_ratio(
        path_to_dwellings=snakemake.input.dwellings,
        regions=snakemake.input.regions,
        nuts_to_regions=snakemake.input.nuts_to_regions,
        out_path=snakemake.output[0]
    )