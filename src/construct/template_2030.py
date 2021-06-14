"""Applies config parameters to 2030 basis template files."""
from pathlib import Path

import numpy as np
import pandas as pd

import jinja2

from eurocalliopelib import filters


def parameterise_template(path_to_template, scaling_factors, heat, path_to_result):
    """Applies config parameters to template files."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    path_to_template = Path(path_to_template)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(path_to_template.parent), lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit
    env.globals["mean"] = np.mean
    rendered = env.get_template(path_to_template.name).render(
        scaling_factors=scaling_factors,
        heat=heat,
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        scaling_factors=snakemake.params["scaling_factors"],
        heat=snakemake.params["heat"],
        path_to_result=snakemake.output[0]
    )
