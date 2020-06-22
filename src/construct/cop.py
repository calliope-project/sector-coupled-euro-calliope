import xarray as xr
import pandas as pd

idx = pd.IndexSlice


def get_cop(
    population, soil_temp, air_temp, annual_consumption,
    heat_tech_params, cop_out_path, annual_consumption_out_path
):
    """
    heat_tech_params include:
        ratio of carnot efficiency to actual efficiency: 'carnot_performance',
        tech efficiencies: 'natural_gas_eff', 'oil_eff', 'coal_eff', 'wood_eff', 'solar_thermal_eff', 'direct_electric_eff',
        working temperatures for COP: 'space_heating_temp', 'water_heating_temp'
    """
    sh_T = heat_tech_params['space_heating_temp'] + 273.15  # give in Kelvin
    wh_T = heat_tech_params['water_heating_temp'] + 273.15  # give in Kelvin

    annual_consumption_df = pd.read_csv(annual_consumption)
    soil_temp_ds = xr.open_dataset(soil_temp)
    air_temp_ds = xr.open_dataset(air_temp)

    cop_dict = {}
    for country in population.country_code.unique():
        sites = population[population.country_code == country].site.values
        soil_temp_df = soil_temp_ds.soil_temperature_5.loc[{'site': sites}].to_pandas()
        air_temp_df = air_temp_ds.temperature.loc[{'site': sites}].to_pandas()
        air_temp_df = air_temp_df.where(
            air_temp_df >= heat_tech_params['degree_day_threshhold']
        )

        population_weights = (
            population.loc[sites, 'population'].transform(lambda x: x / x.sum())
        )
        cop_dict[country] = get_national_cop(
            air_temp_df, soil_temp_df, heat_tech_params['carnot_performance'],
            population_weights, heat_tech_params['sh_T'], heat_tech_params['wh_T']
        )

    # index = [country, year], columns=[cat_name, tech]
    cop_df = pd.concat(cop_dict.values(), keys=cop_dict.keys())
    cop_df = get_air_ground_heat_pump_ratio(annual_consumption, cop_df)
    cop_df.to_csv(cop_out_path)

    # Update annual_consumption to include guessed heat pump electricity consumption
    annual_consumption_df = get_heat_pump_elec_consumption(annual_consumption_df, cop_df)
    annual_consumption_df.to_csv(annual_consumption_out_path)


def get_national_cop(air_temp, soil_temp, carnot_factor, population_weights, sh_T, wh_T):
    """For COP method see [Nouvel_2015]; then, we population weight it"""
    def cop(t_source, t_out, carnot_factor):
        """Returns annual average COP for a country, based on population weighting"""
        return (
            (carnot_factor * t_source / (t_source - t_out))
            .mul(population_weights, axis=0)
            .sum()
            .groupby(t_out.columns.year, axis=1).mean()
        )

    cop_sh_ashp = cop(sh_T, air_temp, carnot_factor)
    cop_sh_gshp = cop(sh_T, soil_temp, carnot_factor)

    cop_wh_ashp = cop(wh_T, air_temp, carnot_factor)
    cop_wh_gshp = cop(wh_T, soil_temp, carnot_factor)

    cop = pd.concat(
        [cop_sh_ashp, cop_sh_gshp, cop_wh_ashp, cop_wh_gshp],
        axis=1, names=[
            ('space_heating', 'ashp'), ('space_heating', 'gshp'),
            ('water_heating', 'ashp'), ('water_heating', 'gshp')
        ]
    )

    cop.columns = pd.MultiIndex.from_tuples(cop.columns).rename(['cat_name', 'tech'])

    return cop


def get_national_cop(air_temp, soil_temp, carnot_factor, population_weights, sh_T, wh_T):
    """For COP method see [Nouvel_2015]; then, we population weight it"""
    def cop(t_source, t_out, carnot_factor):
        """Returns annual average COP for a country, based on population weighting"""
        return (
            (carnot_factor * t_source / (t_source - t_out))
            .mul(population_weights, axis=0)
            .sum()
            .groupby(t_out.columns.year, axis=1).mean()
        )

    cop_sh_ashp = cop(sh_T, air_temp, carnot_factor)
    cop_sh_gshp = cop(sh_T, soil_temp, carnot_factor)

    cop_wh_ashp = cop(wh_T, air_temp, carnot_factor)
    cop_wh_gshp = cop(wh_T, soil_temp, carnot_factor)

    cop = pd.concat(
        [cop_sh_ashp, cop_sh_gshp, cop_wh_ashp, cop_wh_gshp],
        axis=1, names=[
            ('space_heating', 'ashp'), ('space_heating', 'gshp'),
            ('water_heating', 'ashp'), ('water_heating', 'gshp')
        ]
    )

    cop.columns = pd.MultiIndex.from_tuples(cop.columns).rename(['cat_name', 'tech'])

    return cop


def get_air_ground_heat_pump_ratio(annual_consumption, cop):
    """
    Use known Swiss COP to guess the ratio of ASHP to GSHP (since they have different COPs)
    COP = delivered heat / input electricity. delivered heat = input electricity + ambient heat
    COP = 1 + ambient heat / input electricity
    """
    # change to have years and category in index, carrier in columns
    df = annual_consumption.loc[idx[:, :, 'CH', ['ambient_heat', 'heat_pump']], :].stack().unstack('carrier_name')
    # get COP from carrier info, put category name in columns
    known_cop = (df.ambient_heat.div(df.heat_pump, axis=0) + 1).unstack('cat_name')
    # get modelled COP in the form: index = years, columns=[cat_name, tech]
    modelled_cop = cop.loc[idx['CH', :], :].reindex(known_cop.index)

    for cat_name in ['space_heating', 'water_heating']:
        ground_factor = (  # ground_factor = 1 - air_factor; COP = air_factor * cop_ashp + ground_factor * cop_gshp
            known_cop.loc[:, cat_name] -
            modelled_cop.loc[:, idx[cat_name, 'ashp']]
        ).div(
            modelled_cop.loc[:, idx[cat_name, 'gshp']] -
            modelled_cop.loc[:, idx[cat_name, 'ashp']]
        )
        # Give all countries the same gshp and ashp factors
        cop[(cat_name, 'gshp_factor')] = ground_factor.loc[cop.index.get_level_values('year')].values
        cop[(cat_name, 'ashp_factor')] = (1 - ground_factor).loc[cop.index.get_level_values('year')].values

    return cop


def get_heat_pump_elec_consumption(df, cop):
    elec = df.loc[:, 'electricity']
    heat = df.loc[:, 'ambient_heat']
    cop = cop.stack('cat_name').reorder_levels(df.index.names)
    ave_cop = cop.loc[:, 'ashp'] * cop.loc[:, 'ashp_factor'] + cop.loc[:, 'gshp'] * cop.loc[:, 'gshp_factor']
    hp_elec = heat.reindex(df.index).fillna(0).div(ave_cop.reindex(df.index) - 1, axis=0)
    direct_heat_elec = elec.sub(hp_elec)

    df['direct_electric'] = direct_heat_elec
    df['heat_pump'] = hp_elec

    return df


if __name__ == '__main__':
    get_cop(
        population = snakemake.inputs.population,
        soil_temp = snakemake.inputs.soil_temp,
        air_temp = snakemake.inputs.air_temp,
        annual_consumption = snakemake.inputs.annual_consumption,
        heat_tech_params = snakemake.params.heat_tech_params,
        cop_out_path = snakemake.output.cop,
        annual_consumption_out_path = snakemake.output.annual_consumption
    )