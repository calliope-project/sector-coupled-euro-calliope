import geopandas as gpd

import util


def nuts_to_nuts(nuts_orig, nuts_per_region, units):
    """
    Go from the NUTS regions used in euro-calliope (e.g. 2006) to the regions
    used in a Eurostat dataset (e.g. 2010).
    `nuts_per_region` should be a pandas Seires mapping from NUTS regions to a
    larger region set (e.g. new Euro-Calliope regions).
    `nuts_orig` is the raw regions coming from Eurostat, as a geodataframe.
    `units` is the shapes of the regions for which we have `nuts_per_region` information
    """

    # Quicker than 'overlay'
    nuts_orig = nuts_orig[
        (nuts_orig.bounds.minx >= units.total_bounds[0]) &
        (nuts_orig.bounds.miny >= units.total_bounds[1]) &
        (nuts_orig.bounds.maxx <= units.total_bounds[2]) &
        (nuts_orig.bounds.maxy <= units.total_bounds[3]) &
        nuts_orig.CNTR_CODE.isin(units.country_code.apply(util.get_alpha2))
    ]

    # First just map across all those with the same ID
    nuts3_valid_in_both = nuts_orig.merge(
        nuts_per_region.rename('regions'), left_on='id', right_index=True, how='outer'
    ).set_index('id')

    # Then go through the task of overlaying all remaining shapes, ensuring
    # we only match with regions in the correct country, and picking the region
    # to which the greatest share of nuts_orig shape area can be assigned.
    missing_ids = nuts3_valid_in_both[nuts3_valid_in_both.regions.isna()]
    overlayed_units = gpd.overlay(missing_ids, units, keep_geom_type=False).set_index('FID')
    overlayed_units = overlayed_units[
        overlayed_units.CNTR_CODE.apply(util.get_alpha3) == overlayed_units.country_code
    ]

    for i in missing_ids.index:
        if i not in overlayed_units.index:
            print(f'skipping {i} entirely')
            continue

        relevant_units = overlayed_units.loc[[i]]
        if len(relevant_units) == 1:
            nuts3_valid_in_both.loc[i, 'regions'] = relevant_units.id.item()
        else:
            relevant_units['area_weight'] = relevant_units.area / relevant_units.area.sum()
            print(
                'Choosing {} with {}% of {} area'
                .format(relevant_units.sort_values('area_weight').iloc[-1].id,
                        relevant_units['area_weight'].max() * 100, i)
            )
            nuts3_valid_in_both.loc[i, 'regions'] = relevant_units.sort_values('area_weight').iloc[-1].id

    return nuts3_valid_in_both
