import argparse
import glob

import numpy as np
import xarray as xr

import calliope

EXCL_SCORE_COST = 100


def rerun_model(
    dir_path, n_spores, slack, techs=None, tech_groups=None, first_exclude_only=False
):
    calliope.set_log_verbosity()
    excl_path_string = ",".join([i for i in [techs, tech_groups] if i is not None])
    cost_op_model = calliope.read_netcdf(dir_path + "spore_0.nc")
    print("Loaded cost optimal model")
    exclusion_model_paths = glob.glob(dir_path + "spore_excl-{}-*.nc".format(excl_path_string))

    if exclusion_model_paths and not first_exclude_only:
        latest_spore = sorted([int(i.split("-")[-1].replace(".nc", "")) for i in exclusion_model_paths])[-1]
        if latest_spore != 0:
            latest_results = calliope.read_netcdf(dir_path + "spore_excl-{}-{}.nc".format(excl_path_string, latest_spore))
            print("Loaded latest exclusion SPORE model: spore_excl-{}-{}.nc".format(excl_path_string, latest_spore))
    else:
        latest_spore = 0

    cost_op_model = update_run_config_for_spores(cost_op_model, n_spores, slack, dir_path, excl_path_string)
    cost_op_model = setup_cost_opt_model(cost_op_model, techs, tech_groups)

    if latest_spore == 0:
        if first_exclude_only:
            print("Running only first exclude run, without SPORES scores")
            cost_op_model = update_run_config_for_exclusion_only(cost_op_model)
            if "objective_cost_class" in cost_op_model._model_data.data_vars:
                cost_op_model._model_data = (
                    cost_op_model._model_data.drop_vars(["objective_cost_class"])
                )
            cost_op_model.run(force_rerun=True)
            cost_op_model.to_netcdf(dir_path + "spore_excl-{}-0.nc".format(excl_path_string))
            return None  # i.e. we break out of this script
        new_scores = _get_new_scores(cost_op_model._model_data, cost_op_model)
    if latest_spore > 0:
        cost_op_model._model_data["cost_energy_cap"] = latest_results._model_data["cost_energy_cap"]
        new_scores = _get_new_scores(latest_results._model_data, cost_op_model)

    if "spores" not in cost_op_model._model_data.coords:
        cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [latest_spore]))
    else:
        cost_op_model._model_data.coords["spores"] = [latest_spore]

    if "objective_cost_class" in cost_op_model._model_data.data_vars:
        cost_op_model._model_data = cost_op_model._model_data.drop_vars(["objective_cost_class"])

    relevant_loc_techs = _get_relevant_loc_techs(cost_op_model)
    cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score", "loc_techs_investment_cost": relevant_loc_techs}] += new_scores

    assert len((cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score"}].dropna("loc_techs_investment_cost"))) == len(relevant_loc_techs)
    print(f"Running SPORES with slack of {cost_op_model.run_config['spores_options']['slack']}")
    print(f"Running SPORES with options {cost_op_model._model_data.attrs['run_config']}")

    cost_op_model.run(force_rerun=True)


def _get_new_scores(model_data, cost_op_model):
    relevant_loc_techs = _get_relevant_loc_techs(cost_op_model)
    # First, we ensure that we don't erroneously score techs that have to be there no matter what.
    if "energy_cap_min" in cost_op_model._model_data.data_vars.keys():
        print("I am min")
        min_cap = cost_op_model._model_data.energy_cap_min.fillna(0)
    else:
        print("I am not min")
        min_cap = 0

    if "energy_cap_equals" in cost_op_model._model_data.data_vars.keys():
        forced_energy_cap = cost_op_model._model_data.energy_cap_equals.fillna(min_cap)
    else:
        forced_energy_cap = min_cap

    return xr.where(model_data.energy_cap - forced_energy_cap > 1e-3, 100, 0).loc[relevant_loc_techs].drop_vars("loc_techs")


def _get_relevant_loc_techs(cost_op_model):
    return (
        cost_op_model
        ._model_data
        .cost_energy_cap
        .loc[{"costs": "spores_score"}]
        .dropna("loc_techs_investment_cost")
        .loc_techs_investment_cost
    )


def setup_cost_opt_model(cost_op_model, techs, tech_groups):
    cost_op_model._model_data["group_cost_max"].loc[{"costs": "monetary"}] = (
        cost_op_model._model_data.cost.sum().item() *
        (1 + cost_op_model.run_config["spores_options"]["slack"])
    )
    all_techs = set()
    if tech_groups is not None:
        tech_groups = tech_groups.split(",")
        all_techs.update(np.concatenate([
            cost_op_model._model_data.techs.loc[cost_op_model._model_data.inheritance.str.find(tech_group) > -1].values
            for tech_group in tech_groups
        ]))
    if techs is not None:
        all_techs.update(techs.split(","))

    loc_techs = [
        f"{loc}::{tech}"
        for loc in cost_op_model._model_data.locs.values
        for tech in all_techs
    ]
    valid_loc_techs = cost_op_model._model_data.loc_techs_investment_cost.to_index().intersection(loc_techs)
    cost_op_model._model_data.cost_energy_cap.loc[{"costs": "excl_score", "loc_techs_investment_cost": valid_loc_techs}] = EXCL_SCORE_COST
    cost_op_model._model_data.cost_depreciation_rate.loc[{"costs": "excl_score", "loc_techs_investment_cost": valid_loc_techs}] = 1
    cost_op_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}] = (
        cost_op_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}]
        .clip(min=1, max=1)
    )

    cost_op_model._model_data.attrs.pop("objective_function_value", None)
    cost_op_model._model_data.attrs.pop("termination_condition", None)

    print(f"Running with exclusion scores activated for techs {all_techs} and a total of {len(valid_loc_techs)} loc::techs")
    return cost_op_model


def update_run_config_for_spores(cost_op_model, n_spores, slack, dir_path, excl_path_string):
    cost_op_model.run_config["spores_options"]["skip_cost_op"] = True
    cost_op_model.run_config["objective_options"]["cost_class"] = {"spores_score": 1, "excl_score": 10, "monetary": 0}
    cost_op_model.run_config["mode"] = "spores"
    cost_op_model.run_config["solver_io"] = "python"
    cost_op_model.run_config["solver"] = "gurobi_persistent"
    cost_op_model.run_config["spores_options"]["slack"] = slack
    cost_op_model.run_config["spores_options"]["spores_number"] = n_spores
    cost_op_model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_excl-{}-{{}}.nc".format(excl_path_string)
    return cost_op_model


def update_run_config_for_exclusion_only(cost_op_model):
    cost_op_model.run_config["mode"] = "plan"
    cost_op_model.run_config["objective_options"]["cost_class"] = {"excl_score": 1, "monetary": 0, "spores_score": 0}
    return cost_op_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spore_dir", help="Directory to load cost optimal model and save new spores")
    parser.add_argument("--n_spores", help="Number of exclusion spores to run", default=10, type=int)
    parser.add_argument("--slack", help="Cost slack", default=0.1, type=float)
    parser.add_argument("--techs", help="Tech to exclude")
    parser.add_argument("--tech_groups", help="Tech group to exclude")
    parser.add_argument("--first_exclude_only", help="If true, just run an exclusion score minimisation", action="store_true")

    parser.set_defaults(first_exclude_only=False)
    args = parser.parse_args()

    rerun_model(
        args.spore_dir, args.n_spores, args.slack,
        techs=args.techs, tech_groups=args.tech_groups,
        first_exclude_only=args.first_exclude_only
    )
