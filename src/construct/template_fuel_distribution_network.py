import pandas as pd

import jinja2

TEMPLATE = """
overrides:
    synfuel_transmission:
        techs:
            diesel_distribution_export:
                essentials:
                    carrier: diesel
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            methane_distribution_export:
                essentials:
                    carrier: methane
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            diesel_distribution_import:
                essentials:
                    carrier: diesel
                    parent: supply
                constraints:
                    resource: .inf

            methane_distribution_import:
                essentials:
                    carrier: methane
                    parent: supply
                constraints:
                    resource: .inf

        locations:
            {{ regions|join(",") }}:
                techs:
                    diesel_distribution_export:
                    methane_distribution_export:
                    diesel_distribution_import:
                    methane_distribution_import:


"""

def parametrise_template(path_to_regions, path_to_output):
    regions = pd.read_csv(path_to_regions, index_col=0).index
    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)

    rendered = env.from_string(TEMPLATE).render(
        regions=regions
    )
    with open(path_to_output, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parametrise_template(
        path_to_regions=snakemake.input.regions,
        path_to_output=snakemake.output[0]
    )
