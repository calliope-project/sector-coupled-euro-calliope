import jinja2
import pandas as pd
import pycountry

TEMPLATE = """
links:
    {% for row_id, row in gtcs.iterrows() %}
    {{ row.Link_from }},{{ row.Link_to }}.techs:
    {%- if row.ac > 0 %}
        ac_transmission:
            constraints:
                energy_cap_equals: {{ row.ac }} # [{{ pow_scaling_factor }} MW]
    {%- endif -%}
    {% if row.dc > 0 %}
        dc_transmission:
            constraints:
                energy_cap_equals: {{ row.dc }} # [{{ pow_scaling_factor }} MW]
    {%- endif -%}
    {% endfor %}
"""
COUNTRY_CODE_COLUMN = "country_code"



def generate_links(path_to_gtc, scaling_factor, path_to_result, resolution):
    """Generate a file that represents links in Calliope."""
    gtcs = _read_gtcs(path_to_gtc, scaling_factor, resolution)
    links = jinja2.Template(TEMPLATE).render(
        gtcs=gtcs,
        pow_scaling_factor=1 / scaling_factor
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _read_gtcs(path_to_gtc, scaling_factor, resolution):
    def _get_alpha3(alpha2):
        if alpha2 == 'UK':
            alpha2 = 'GB'
        return pycountry.countries.get(alpha_2=alpha2).alpha_3

    gtcs = pd.read_excel(path_to_gtc, sheet_name="links_internal", header=0)

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
        gtcs[['Link_from', 'Link_to', "('GTC', 'AC')", "('GTC', 'DC')"]]
        .groupby(['Link_from', 'Link_to']).sum()
    ).rename(columns={"('GTC', 'AC')": 'ac', "('GTC', 'DC')": 'dc'})
    # Remove reference to internal links
    gtcs = gtcs.loc[gtcs.index.get_level_values(0) != gtcs.index.get_level_values(1)]
    gtcs
    return gtcs.astype(pd.np.float32).mul(scaling_factor).reset_index()


if __name__ == "__main__":
    generate_links(
        path_to_gtc=snakemake.input.gtc,
        path_to_result=snakemake.output[0],
        scaling_factor=snakemake.params.scaling_factor,
        resolution=snakemake.wildcards[0]
    )

