import argparse
import glob

import numpy as np
import xarray as xr

import calliope

EXCL_SCORE_COST = 100


def rerun_model(dir_path, n_spores, techs=None, tech_groups=None):
    calliope.set_log_verbosity()
    excl_path_string = ",".join([i for i in [techs, tech_groups] if i is not None])
    cost_op_model = calliope.read_netcdf(dir_path + "spore_0.nc")
    print("Loaded cost optimal model")
    exclusion_model_paths = glob.glob(dir_path + "spore_excl-{}-*.nc".format(excl_path_string))

    if exclusion_model_paths:
        latest_spore = sorted([int(i.split("-")[-1].replace(".nc", "")) for i in exclusion_model_paths])[-1]
        latest_results = calliope.read_netcdf(dir_path + "spore_excl-{}-{}.nc".format(excl_path_string, latest_spore))
        print("Loaded latest exclusion SPORE model: spore_excl-{}-{}.nc".format(excl_path_string, latest_spore))
    else:
        latest_spore = 0

    cost_op_model = update_run_config(cost_op_model, n_spores, dir_path, excl_path_string)
    cost_op_model = setup_cost_opt_model(cost_op_model, techs, tech_groups)

    if latest_spore == 0:
        new_scores = _get_new_scores(cost_op_model._model_data)
    if latest_spore > 0:
        cost_op_model._model_data["cost_energy_cap"] = latest_results._model_data["cost_energy_cap"]
        new_scores = _get_new_scores(latest_results._model_data)

    if "spores" not in cost_op_model._model_data.coords:
        cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [latest_spore]))
    else:
        cost_op_model._model_data.coords["spores"] = [latest_spore]

    if "objective_cost_class" in cost_op_model._model_data.data_vars:
        cost_op_model._model_data = cost_op_model._model_data.drop_vars(["objective_cost_class"])

    cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score"}] += new_scores

    cost_op_model.run(force_rerun=True)


def _get_new_scores(model_data):
    return xr.where(model_data.energy_cap > 1e-3, 100, 0).loc[model_data.loc_techs_investment_cost].drop_vars("loc_techs")


def setup_cost_opt_model(cost_op_model, techs, tech_groups):
    cost_op_model._model_data["group_cost_max"].loc[{"costs": "monetary"}] = (
        cost_op_model._model_data.cost.sum().item() *
        (1 + cost_op_model.run_config["spores_options"]["slack"])
    )
    if tech_groups is not None and techs is None:
        tech_groups = tech_groups.split(",")
        techs = np.unique(np.concatenate([
            cost_op_model._model_data.techs.loc[cost_op_model._model_data.inheritance.str.find(tech_group) > -1].values
            for tech_group in tech_groups
        ]))
    elif techs is not None:
        techs = techs.split(",")

    loc_techs = [
        f"{loc}::{tech}"
        for loc in cost_op_model._model_data.locs.values
        for tech in techs
    ]
    valid_loc_techs = cost_op_model._model_data.loc_techs_investment_cost.to_index().intersection(loc_techs)
    cost_op_model._model_data.cost_energy_cap.loc[{"costs": "excl_score", "loc_techs_investment_cost": valid_loc_techs}] = EXCL_SCORE_COST

    cost_op_model._model_data.attrs.pop("objective_function_value", None)
    cost_op_model._model_data.attrs.pop("termination_condition", None)

    print(f"Running with exclusion scores activated for techs {techs} and a total of {len(loc_techs)} loc::techs")
    return cost_op_model


def update_run_config(cost_op_model, n_spores, dir_path, excl_path_string):
    cost_op_model.run_config["spores_options"]["skip_cost_op"] = True
    cost_op_model.run_config["objective_options"]["cost_class"] = {"spores_score": 1, "excl_score": 10, "monetary": 0}
    cost_op_model.run_config["spores_options"]["spores_number"] = n_spores
    cost_op_model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_excl-{}-{{}}.nc".format(excl_path_string)
    return cost_op_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spore_dir", help="Directory to load cost optimal model and save new spores")
    parser.add_argument("--n_spores", help="Number of exclusion spores to run", default=10, type=int)
    parser.add_argument("--techs", help="Tech to exclude")
    parser.add_argument("--tech_groups", help="Tech group to exclude")
    args = parser.parse_args()

    rerun_model(
        args.spore_dir, args.n_spores,
        techs=args.techs, tech_groups=args.tech_groups
    )
