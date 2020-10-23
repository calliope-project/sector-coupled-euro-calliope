import jinja2
import numpy as np
import pandas as pd
import pycountry

TEMPLATE = """
overrides:
    outer-countries:
        techs.outer_country_supply:
            essentials:
                carrier_out: electricity
                parent: supply
                name: Neighbouring country supply
            constraints:
                resource: .inf
                energy_cap_max: .inf

        locations:
        {%- for row_id in gtcs.Link_from.unique() %}
            {{ row_id }}:
                techs: {outer_country_supply}
                coordinates: {{ coords[row_id.split('_')[0]] }}
        {%- endfor %}

        links:
        {%- for row_id, row in gtcs.iterrows() %}
            {{ row.Link_from }},{{ row.Link_to }}.techs:
                ac_transmission:
                    constraints:
                        energy_cap_equals: {{ row.ac }} # [{{ pow_scaling_factor }} MW]
                    distance: {{ row.Length }}  # [{{ dist_scaling_factor }} m]
        {%- endfor %}
"""

# these are not the centroids of the respective countries, but locations closer to Europe, for better viz
OUTER_COUNTRY_COORDS = {
    'MAR': {'lat': 33.693266, 'lon': -5.026300},
    'DZA': {'lat': 34.209077, 'lon': 3.092407},
    'TUN': {'lat': 35.792865, 'lon': 9.460798},
    'LBY': {'lat': 29.558167, 'lon': 17.075974},
    'EGY': {'lat': 29.710955, 'lon': 31.151270},
    'TUR': {'lat': 40.867002, 'lon': 29.831711},
    'RUS': {'lat': 55.593247, 'lon': 37.306013},
    'UKR': {'lat': 50.362284, 'lon': 30.620247}
}


def generate_override(path_to_gtc, scaling_factors, path_to_result, resolution):
    """Generate a file that represents an override in Calliope."""
    gtcs = _read_gtcs(path_to_gtc, scaling_factors, resolution)

    override = jinja2.Template(TEMPLATE).render(
        coords=OUTER_COUNTRY_COORDS,
        gtcs=gtcs,
        pow_scaling_factor=1 / scaling_factors['power'],
        dist_scaling_factor=1 / scaling_factors['distance']
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(override)


def _read_gtcs(path_to_gtc, scaling_factors, resolution):
    def _get_alpha3(alpha2):
        if alpha2 == 'UK':
            alpha2 = 'GB'
        return pycountry.countries.get(alpha_2=alpha2).alpha_3

    gtcs = pd.read_excel(path_to_gtc, sheet_name="links_external", header=0)

    locations = pd.read_excel(path_to_gtc, sheet_name="locations", header=0, index_col=0)
    locations.index = locations.index.str.lower()
    locations['Country'] = locations.Country.map(_get_alpha3)

    if resolution == 'national':
        gtcs['Link_from'] = gtcs.Link_from.replace(locations['Country'])
        gtcs['Link_to'] = gtcs.Link_to.replace(locations['Country'])
    elif resolution == 'eurospores':
        gtcs['Link_from'] = gtcs.Link_from.replace(locations['EuroSPORES'])
        gtcs['Link_to'] = gtcs.Link_to.replace(locations['EuroSPORES'])
    gtcs = (
        gtcs[['Link_from', 'Link_to', 'Length', '100% RES']]
        .groupby(['Link_from', 'Link_to']).sum()
    ).rename(columns={"100% RES": 'ac'})
    # Remove reference to internal links
    gtcs = gtcs.loc[gtcs.index.get_level_values(0) != gtcs.index.get_level_values(1)].astype(np.float32)

    gtcs['ac'] = gtcs['ac'].mul(scaling_factors['power'])
    gtcs['Length'] = gtcs['Length'].mul(scaling_factors['distance'])

    return gtcs.reset_index()


if __name__ == "__main__":
    generate_override(
        path_to_gtc=snakemake.input.gtc,
        path_to_result=snakemake.output[0],
        scaling_factors=snakemake.params.scaling_factors,
        resolution=snakemake.wildcards[0]
    )

