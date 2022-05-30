import argparse
import os


import xarray as xr

import calliope

EXCL_SCORE_COST = 100
COST_OPT_VAL = 1003.5068965694

from spores_exclude import rerun_model, setup_cost_opt_model, update_run_config_for_exclusion_only

def run_updated_demand(
    dir_path, path_to_old_demand_model, path_to_demand_scales,
    n_spores, slack, techs=None, tech_groups=None, first_exclude_only=False
):
    print(dir_path)
    calliope.set_log_verbosity()
    if "spore_0.nc" in os.listdir(dir_path):
        print("Running SPORES from cost optimal model")
        rerun_model(
            dir_path, n_spores, slack, techs, tech_groups, first_exclude_only
        )
    elif first_exclude_only:
        print("Running only first exclude run, without SPORES scores")
        new_demand_model = prep_new_demand_model(
            dir_path, path_to_old_demand_model, path_to_demand_scales, n_spores, slack
        )
        new_demand_model._model_data["cost"].loc[{"costs":"monetary"}] *= (
            COST_OPT_VAL / new_demand_model._model_data["cost"].loc[{"costs":"monetary"}].sum()
        )
        new_demand_model = setup_cost_opt_model(new_demand_model, techs, tech_groups)
        new_demand_model = update_run_config_for_exclusion_only(new_demand_model)
        new_demand_model.run(force_rerun=True)
        excl_path_string = ",".join([i for i in [techs, tech_groups] if i is not None])
        new_demand_model.to_netcdf(dir_path + "spore_excl-{}-0.nc".format(excl_path_string))

    else:
        print("Running SPORES from scratch (inc. optimising for cost first)")
        run_from_scratch(
            dir_path, path_to_old_demand_model, path_to_demand_scales, n_spores, slack
        )

def run_from_scratch(
    dir_path, path_to_old_demand_model, path_to_demand_scales, n_spores, slack
):

    new_demand_model = prep_new_demand_model(
        dir_path, path_to_old_demand_model, path_to_demand_scales, n_spores, slack
    )
    if hasattr(new_demand_model, "results"):
        to_delete = list(new_demand_model.results.data_vars.keys())
        new_demand_model._model_data = new_demand_model._model_data.drop_vars(to_delete)
        new_demand_model = calliope.Model(config=None, model_data=new_demand_model._model_data)

    assert not hasattr(new_demand_model, "results")
    print("Running prepared new demand model from scratch")
    print(new_demand_model.run_config)

    new_demand_model.run(force_rerun=True)


def prep_new_demand_model(dir_path, path_to_old_demand_model, path_to_demand_scales, n_spores, slack):
    old_demand_model = calliope.read_netcdf(path_to_old_demand_model)
    demand_scales = xr.open_dataset(path_to_demand_scales)

    new_demand_vars = old_demand_model.inputs * demand_scales.reindex_like(old_demand_model.inputs).fillna(1)
    new_demand_model = update_run_config_for_spores(old_demand_model, n_spores, slack, dir_path)

    for k, v in new_demand_vars.items():
        v.attrs = new_demand_model._model_data[k].attrs
        new_demand_model._model_data[k] = v

    if "objective_cost_class" in new_demand_model._model_data.data_vars:
        new_demand_model._model_data = new_demand_model._model_data.drop_vars(["objective_cost_class"])
        new_demand_model = calliope.Model(config=None, model_data=new_demand_model._model_data)

    print("Prepared new demand model")

    new_demand_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}] = (
        new_demand_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}]
        .clip(min=1, max=1)
    )

    return new_demand_model


def update_run_config_for_spores(model, n_spores, slack, dir_path):
    model.run_config["objective_options"]["cost_class"] = {"spores_score": 0, "excl_score": 0, "monetary": 1}
    model.run_config["mode"] = "spores"
    model.run_config["solver_io"] = "python"
    model.run_config["solver"] = "gurobi_persistent"
    model.run_config["solver_options"]["BarHomogeneous"] = 1
    model.run_config["spores_options"]["slack"] = slack
    model.run_config["spores_options"]["spores_number"] = n_spores
    model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_{}.nc"

    model._model_data.attrs.pop("objective_function_value", None)
    model._model_data.attrs.pop("termination_condition", None)

    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spore_dir", help="Directory to load cost optimal model and save new spores")
    parser.add_argument("path_to_old_demand_model", help=".nc file built to today's end-use service demands")
    parser.add_argument("path_to_demand_scales", help=".nc file with ratios of demand change to apply")
    parser.add_argument("--n_spores", help="Number of exclusion spores to run", default=10, type=int)
    parser.add_argument("--slack", help="Cost slack", default=0.1, type=float)
    parser.add_argument("--techs", help="Tech to exclude")
    parser.add_argument("--tech_groups", help="Tech group to exclude")
    parser.add_argument("--first_exclude_only", help="If true, just run an exclusion score minimisation", action="store_true")

    parser.set_defaults(first_exclude_only=False)
    args = parser.parse_args()

    run_updated_demand(
        args.spore_dir, args.path_to_old_demand_model, args.path_to_demand_scales,
        args.n_spores, args.slack,
        techs=args.techs, tech_groups=args.tech_groups,
        first_exclude_only=args.first_exclude_only
    )
