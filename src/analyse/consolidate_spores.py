from pathlib import Path
import os

import pandas as pd

from friendly_calliope.io import write_dpkg


def consolidate_spores(
    path_to_cost_optimal_metrics, paths_to_spores_metrics, output_dir
):
    print(paths_to_spores_metrics)
    cost_opt_metrics = get_metrics(path_to_cost_optimal_metrics)
    spores_metrics = {
        Path(spore): get_metrics(spore) for spore in paths_to_spores_metrics
    }
    processed_spores = pd.Series(list(spores_metrics.keys()))

    spores_metrics_grouped = {
        metric: pd.concat(
            [spore[metric] for spore in spores_metrics.values()],
            keys=processed_spores.index, names=["scenario"]
        )
        for metric in cost_opt_metrics.keys()
    }
    _metrics_to_dpkg(spores_metrics_grouped, output_dir, "spores")
    _metrics_to_dpkg(cost_opt_metrics, output_dir, "cost_opt")

    processed_spores.to_csv(os.path.join(output_dir, "processed_spores.csv"))


def get_metrics(path_to_metrics):
    metrics = {}
    for metric_path in Path(path_to_metrics).glob("*.csv"):
        metric_df = pd.read_csv(metric_path, index_col=False)
        metric_series = metric_df.set_index(list(metric_df.columns[:-1])).squeeze()
        assert isinstance(metric_series, pd.Series)
        metrics[metric_series.name] = metric_series
    return metrics


def _metrics_to_dpkg(metrics, output_dir, spore_or_cost_opt):
    if spore_or_cost_opt == "spores":
        meta_name = "euro-calliope-spores-results"
        meta_description = "Calliope Euro-SPORES output dataset"
        extra_keywords = []
    elif spore_or_cost_opt == "cost_opt":
        meta_name = "euro-calliope-cost-opt-data"
        meta_description = "Calliope SPORES cost optimal baseline output dataset"
        extra_keywords = ["cost-optimal"]

    # FIXME: remove when friendly_data can handle this.
    for metric_name, metric_series in metrics.items():
        if not set(metric_series.index.names).intersection(["locs", "techs", "scenario"]):
            metrics[metric_name] = metrics[metric_name].to_frame("foo").rename_axis(columns="scenario").stack()

    meta = {
        "name": meta_name,
        "description": meta_description,
        "keywords": ["calliope", "SPORES", "2h_resolution"] + extra_keywords,
        "license": "CC-BY-4.0"
    }
    write_dpkg(metrics, os.path.join(output_dir, spore_or_cost_opt), meta)


if __name__ == '__main__':
    consolidate_spores(
        path_to_cost_optimal_metrics=snakemake.input.cost_optimal_model,
        paths_to_spores_metrics=snakemake.input.spores,
        output_dir=snakemake.output[0]
    )
