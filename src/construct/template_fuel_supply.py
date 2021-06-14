"""Applies config parameters to 2030 basis template files."""
from pathlib import Path

import numpy as np
import pandas as pd

import jinja2

from eurocalliopelib import filters


def parameterise_template(path_to_template, path_to_fuel_costs, scaling_factors, fuel_cost_source, fuel_cost_year, path_to_result):
    """Applies config parameters to template files."""

    fuel_cost_df = pd.read_excel(path_to_fuel_costs, sheet_name="Appendix 1", header=0, index_col=[0, 1, 2, 3])
    fuel_cost_series = fuel_cost_df.xs((fuel_cost_source, fuel_cost_year, "2015 REF"), level=(1, 2, 3)).loc[:, "EU"]
    assert isinstance(fuel_cost_series, pd.Series)
    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    path_to_template = Path(path_to_template)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(path_to_template.parent), lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit
    env.globals["mean"] = np.mean
    rendered = env.get_template(path_to_template.name).render(
        scaling_factors=scaling_factors,
        fuel_cost=fuel_cost_series
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        path_to_fuel_costs=snakemake.input.fuel_costs,
        scaling_factors=snakemake.params["scaling_factors"],
        fuel_cost_source=snakemake.params["fuel_cost_source"],
        fuel_cost_year=snakemake.params["fuel_cost_year"],
        path_to_result=snakemake.output[0]
    )
