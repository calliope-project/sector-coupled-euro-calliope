import sys

import jinja2
import numpy as np
import pandas as pd

sys.path.append('euro-calliope/src/')
import filters


TEMPLATE = """
links:
    {% for row_id, row in gtcs.iterrows() %}
    {{ row_id[0] }},{{ row_id[1] }}.techs:
        {% for link, data in row.iteritems() %}
        {% if data > 0 %}
        {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission:
            constraints:
                energy_cap_min: {{ data * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
            costs.monetary.energy_cap: {{ costs["{}-{}".format(link[0], link[1].lower())] * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR/MW") }}
        {% if link[2] == "one-way" %}
            one_way: true
        {% endif %}
        {% endif %}
        {% endfor %}
    {% endfor %}
"""


def generate_links(path_to_gtc, scaling_factors, path_to_result, resolution, costs):
    """Generate a file that represents links in Calliope."""
    gtcs = _read_gtcs(path_to_gtc, resolution)

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    links = env.from_string(TEMPLATE).render(
        gtcs=gtcs,
        scaling_factors=scaling_factors,
        costs=costs
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _read_gtcs(path_to_gtc, resolution):

    gtcs = (
        pd.read_csv(path_to_gtc, header=0)
        .set_index(['from', 'to', 'current', 'type', 'direction'])
        .drop('comment', axis=1)
        .squeeze()
        .unstack(['current', 'type', 'direction'])
    )
    if resolution == 'national':
        gtcs.index = gtcs.index.map(lambda x: tuple(i.split('_')[0] for i in x))
        # Remove internal links and sum links across the same borders
        gtcs = gtcs.groupby(level=[0, 1]).sum().drop([i for i in gtcs.index.values if i[0] == i[1]])

    return gtcs


if __name__ == "__main__":
    generate_links(
        path_to_gtc=snakemake.input.gtc,
        path_to_result=snakemake.output[0],
        scaling_factors=snakemake.params.scaling_factors,
        costs=snakemake.params.costs,
        resolution=snakemake.wildcards[0]
    )
