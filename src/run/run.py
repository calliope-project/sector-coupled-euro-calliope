import os

import pyomo.core as po
import numpy as np
import calliope

from calliope.backend.pyomo.util import (
    get_param,
    split_comma_list,
    get_timestep_weight,
    invalid,
)


def run_model(path_to_model, path_to_output):
    current_path = os.path.abspath(os.curdir)
    os.chdir(os.environ["TMPDIR"])

    calliope.set_log_verbosity("info", include_solver_output=True, capture_warnings=True)
    model = calliope.read_netcdf(os.path.join(current_path, path_to_model))

    model.run(build_only=True)

    add_eurocalliope_constraints(model)
    new_model = model.backend.rerun()

    if new_model.results.attrs.get('termination_condition', None) not in ['optimal', 'feasible']:
        calliope.exceptions.BackendError("Problem is non optimal, not saving anything.")

    new_model.to_netcdf(os.path.join(current_path, path_to_output))


def add_eurocalliope_constraints(model):
    backend_model = model._backend_model
    if "energy_cap_max_time_varying" in model._model_data.data_vars:
        print("Building production_max_time_varying constraint")
        add_production_max_time_varying_constraint(model, backend_model)
    if "cb" in model._model_data.data_vars:
        print("Building chp_cb constraint")
        add_chp_cb_constraint(model, backend_model)
    if "cv" in model._model_data.data_vars:
        print("Building chp_cv constraint")
        add_chp_cv_constraint(model, backend_model)
    if "energy_cap_ratio" in model._model_data.data_vars and "cb" in model._model_data.data_vars:
        print("Building chp_p2h constraint")
        add_chp_p2h_constraint(model, backend_model)
    if any(var.startswith("capacity_factor") for var in model._model_data.data_vars):
        print("Building capacity_factor constraint")
        add_capacity_factor_constraints(model, backend_model)
    if any(var.startswith("carrier_prod_per_month") for var in model._model_data.data_vars):
        print("Building carrier_prod_per_month constraint")
        add_carrier_prod_per_month_constraints(model, backend_model)


def equalizer(lhs, rhs, sign):
    if sign == "max":
        return lhs <= rhs
    elif sign == "min":
        return lhs >= rhs
    elif sign == "equals":
        return lhs == rhs
    else:
        raise ValueError("Invalid sign: {}".format(sign))


def add_production_max_time_varying_constraint(model, backend_model):

    def _carrier_production_max_time_varying_constraint_rule(
        backend_model, loc_tech, timestep
    ):
        """
        Set maximum carrier production for technologies with time varying maximum capacity
        """
        energy_cap_max = backend_model.energy_cap_max_time_varying[loc_tech, timestep]
        if invalid(energy_cap_max):
            return po.Constraint.Skip
        model_data_dict = backend_model.__calliope_model_data["data"]
        timestep_resolution = backend_model.timestep_resolution[timestep]
        loc_tech_carriers_out = split_comma_list(
            model_data_dict["lookup_loc_techs_conversion_plus"]["out", loc_tech]
        )

        carrier_prod = sum(
            backend_model.carrier_prod[loc_tech_carrier, timestep]
            for loc_tech_carrier in loc_tech_carriers_out
        )
        return carrier_prod <= (
            backend_model.energy_cap[loc_tech] * timestep_resolution * energy_cap_max
        )

    backend_model.loc_tech_carrier_production_max_time_varying_constraint = po.Set(
        initialize=[
            loc_tech for loc_tech in backend_model.loc_techs_conversion_plus
            if model.inputs.energy_cap_max_time_varying.loc[{"loc_techs": loc_tech}].notnull().all()
        ],
        ordered=True
    )
    model.backend.add_constraint(
        "carrier_production_max_time_varying_constraint",
        ["loc_tech_carrier_production_max_time_varying_constraint", "timesteps"],
        _carrier_production_max_time_varying_constraint_rule,
    )


def add_chp_cb_constraint(model, backend_model):
    def _chp_extraction_cb_constraint_rule(backend_model, loc_tech, timestep):
        """
        Set backpressure line for CHP plants with extraction/condensing turbine
        """
        model_data_dict = backend_model.__calliope_model_data
        loc_tech_carrier_out = model_data_dict["data"]["lookup_loc_techs_conversion_plus"][
            ("out", loc_tech)
        ]
        loc_tech_carrier_out_2 = model_data_dict["data"][
            "lookup_loc_techs_conversion_plus"
        ][("out_2", loc_tech)]

        power_to_heat_ratio = get_param(backend_model, "cb", (loc_tech))

        return backend_model.carrier_prod[loc_tech_carrier_out, timestep] >= (
            backend_model.carrier_prod[loc_tech_carrier_out_2, timestep]
            * power_to_heat_ratio
        )
    backend_model.loc_techs_chp_extraction_cb_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs_conversion_plus
            if model._model_data.cb.loc[{"loc_techs": loc_tech}].notnull()
            and (
                "energy_cap_ratio" not in model._model_data.data_vars or
                model._model_data.energy_cap_ratio.loc[{"loc_techs": loc_tech}].isnull()
            )
        ],
        ordered=True
    )
    model.backend.add_constraint(
        "chp_extraction_cb_constraint",
        ["loc_techs_chp_extraction_cb_constraint", "timesteps"],
        _chp_extraction_cb_constraint_rule,
    )


def add_chp_cv_constraint(model, backend_model):
    def _chp_extraction_cv_constraint_rule(backend_model, loc_tech, timestep):
        """
        Set extraction line for CHP plants with extraction/condensing turbine
        """
        model_data_dict = backend_model.__calliope_model_data
        loc_tech_carrier_out = model_data_dict["data"]["lookup_loc_techs_conversion_plus"][
            ("out", loc_tech)
        ]
        loc_tech_carrier_out_2 = model_data_dict["data"][
            "lookup_loc_techs_conversion_plus"
        ][("out_2", loc_tech)]

        power_loss_factor = get_param(backend_model, "cv", (loc_tech))

        return backend_model.carrier_prod[loc_tech_carrier_out, timestep] <= (
            backend_model.energy_cap[loc_tech]
            - backend_model.carrier_prod[loc_tech_carrier_out_2, timestep]
            * power_loss_factor
        )
    backend_model.loc_techs_chp_extraction_cv_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs_conversion_plus
            if model._model_data.cv.loc[{"loc_techs": loc_tech}].notnull()
        ],
        ordered=True
    )
    model.backend.add_constraint(
        "chp_extraction_cv_constraint",
        ["loc_techs_chp_extraction_cv_constraint", "timesteps"],
        _chp_extraction_cv_constraint_rule,
    )


def add_chp_p2h_constraint(model, backend_model):
    def _chp_extraction_p2h_constraint_rule(backend_model, loc_tech, timestep):
        """
        Set power-to-heat tail for CHPs that allow trading off power output for heat
        """
        model_data_dict = backend_model.__calliope_model_data
        loc_tech_carrier_out = model_data_dict["data"]["lookup_loc_techs_conversion_plus"][
            ("out", loc_tech)
        ]
        loc_tech_carrier_out_2 = model_data_dict["data"][
            "lookup_loc_techs_conversion_plus"
        ][("out_2", loc_tech)]

        power_to_heat_ratio = get_param(backend_model, "cb", loc_tech)
        energy_cap_ratio = get_param(
            backend_model, "energy_cap_ratio", ("out_2", loc_tech_carrier_out_2)
        )
        slope = power_to_heat_ratio / (energy_cap_ratio - 1)
        return backend_model.carrier_prod[loc_tech_carrier_out, timestep] <= (
            slope
            * (
                backend_model.energy_cap[loc_tech] * energy_cap_ratio
                - backend_model.carrier_prod[loc_tech_carrier_out_2, timestep]
            )
        )
    backend_model.loc_techs_chp_extraction_p2h_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs_conversion_plus
            if model._model_data.cb.loc[{"loc_techs": loc_tech}].notnull()
            and (
                "energy_cap_ratio" in model._model_data.data_vars and
                model._model_data.energy_cap_ratio.loc[{"loc_techs": loc_tech}].notnull()
            )
        ],
        ordered=True
    )

    model.backend.add_constraint(
        "chp_extraction_p2h_constraint",
        ["loc_techs_chp_extraction_p2h_constraint", "timesteps"],
        _chp_extraction_p2h_constraint_rule,
    )


def add_capacity_factor_constraints(model, backend_model):

    def _capacity_factor_min_constraint_rule(backend_model, loc_tech):
        """
        If there is capacity of a technology, force the annual capacity factor to be
        at least a certain amount
        """
        return _capacity_factor_constraint_rule_factory(backend_model, loc_tech, "min")

    def _capacity_factor_max_constraint_rule(backend_model, loc_tech):
        """
        If there is capacity of a technology, force the annual capacity factor to be
        at most a certain amount
        """
        return _capacity_factor_constraint_rule_factory(backend_model, loc_tech, "max")

    def _capacity_factor_constraint_rule_factory(backend_model, loc_tech, sense):
        """
        If there is capacity of a technology, force the annual capacity factor to be
        at most (sense="max") or at least (sense="min") a certain amount
        """
        capacity_factor = get_param(backend_model, f"capacity_factor_{sense}", (loc_tech))
        if invalid(capacity_factor):
            return po.Constraint.Skip
        model_data_dict = backend_model.__calliope_model_data
        loc_tech_carrier = model_data_dict["data"]["lookup_loc_techs"][loc_tech]
        lhs = sum(
            backend_model.carrier_prod[loc_tech_carrier, timestep]
            * backend_model.timestep_weights[timestep]
            for timestep in backend_model.timesteps
        )
        rhs = (
            backend_model.energy_cap[loc_tech]
            * capacity_factor
            * get_timestep_weight(backend_model)
            * 8760
        )
        return equalizer(lhs, rhs, sense)

    backend_model.loc_tech_carriers_capacity_factor_min_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs
            if model._model_data.capacity_factor_min.loc[{"loc_techs": loc_tech}].notnull()
        ],
        ordered=True
    )
    backend_model.loc_tech_carriers_capacity_factor_max_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs
            if model._model_data.capacity_factor_max.loc[{"loc_techs": loc_tech}].notnull()
        ],
        ordered=True
    )
    model.backend.add_constraint(
        "capacity_factor_min_constraint",
        ["loc_tech_carriers_capacity_factor_min_constraint"],
        _capacity_factor_min_constraint_rule,
    )
    model.backend.add_constraint(
        "capacity_factor_max_constraint",
        ["loc_tech_carriers_capacity_factor_max_constraint"],
        _capacity_factor_max_constraint_rule,
    )


def add_carrier_prod_per_month_constraints(model, backend_model):

    def _carrier_prod_per_month_constraint_rule_generator(sense):
        def __carrier_prod_per_month_constraint_rule(backend_model, loc_tech, month):
            """
            Set the min/max amount of carrier consumption (relative to annual consumption)
            for a specific loc tech that must take place in a given calender month in the model
            """
            model_data_dict = backend_model.__calliope_model_data
            loc_tech_carrier = model_data_dict["data"]["lookup_loc_techs_conversion"][("out", loc_tech)]

            prod = backend_model.carrier_prod
            prod_total = sum(
                prod[loc_tech_carrier, timestep]
                for timestep in backend_model.timesteps
            )
            prod_month = sum(
                prod[loc_tech_carrier, timestep]
                for timestep in backend_model.timesteps
                if backend_model.month_numbers[timestep].value == month
            )
            if "timesteps" in [
                i.name
                for i in getattr(
                    backend_model, f"carrier_prod_per_month_{sense}_time_varying"
                )._index.subsets()
            ]:
                prod_fraction = sum(
                    get_param(
                        backend_model,
                        f"carrier_prod_per_month_{sense}_time_varying",
                        (loc_tech, timestep),
                    )
                    * backend_model.timestep_resolution[timestep]
                    for timestep in backend_model.timesteps
                    if backend_model.month_numbers[timestep].value == month
                )
            else:
                prod_fraction = get_param(
                    backend_model, f"carrier_prod_per_month_{sense}", (loc_tech)
                )

            return equalizer(prod_month, prod_total * prod_fraction, sense)
        return __carrier_prod_per_month_constraint_rule

    backend_model.months = po.Set(
        initialize=np.unique(model._model_data.timesteps.dt.month.values), ordered=True
    )
    month_numbers = model._model_data.timesteps.dt.month.to_series()
    month_numbers.index = month_numbers.index.strftime("%Y-%m-%d %H:%M")

    backend_model.month_numbers = po.Param(
        backend_model.timesteps,
        initialize=month_numbers.to_dict(),
        mutable=True,
        within=po.Reals,
    )
    backend_model.__calliope_datetime_data.add(("data_vars", "month_numbers"))

    for sense in ["min", "max", "equals"]:
        if f"carrier_prod_per_month_{sense}_time_varying" in model._model_data.data_vars:
            setattr(
                backend_model,
                f"loc_techs_carrier_prod_per_month_{sense}",
                po.Set(
                    initialize=[
                        loc_tech
                        for loc_tech in backend_model.loc_techs
                        if (
                            model
                            ._model_data[f"carrier_prod_per_month_{sense}_time_varying"]
                            .loc[{"loc_techs": loc_tech}]
                            .notnull()
                            .all()
                        )
                    ],
                    ordered=True
                )
            )
            model.backend.add_constraint(
                f"carrier_prod_per_month_{sense}_constraint",
                [f"loc_techs_carrier_prod_per_month_{sense}", "months"],
                _carrier_prod_per_month_constraint_rule_generator(sense),
            )


if __name__ == "__main__":
    run_model(
        path_to_model=snakemake.input.model,
        path_to_output=snakemake.output[0]
    )
