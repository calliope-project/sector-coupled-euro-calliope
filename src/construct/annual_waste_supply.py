import pandas as pd
import geopandas as gpd

import util
from annual_subnational_demand import align_and_scale, get_population_intensity

idx = pd.IndexSlice


def annual_waste_supply(
    path_to_annual_energy_balances, path_to_population, path_to_units,
    renewable_out_path, total_out_path
):
    """
    Get total waste
    """
    energy_balance = util.read_tdf(path_to_annual_energy_balances)
    units = gpd.read_file(path_to_units)
    population_intensity = get_population_intensity(path_to_population, units)

    renewable_waste = get_disaggregated_waste_consumption_for_chp(
        energy_balance, population_intensity, units, "W6210"
    )
    non_renewable_waste = get_disaggregated_waste_consumption_for_chp(
        energy_balance, population_intensity, units, "W6100_6220"
    )
    total_waste_consumption = non_renewable_waste.add(renewable_waste, fill_value=0)

    renewable_waste.to_csv(renewable_out_path)
    total_waste_consumption.to_csv(total_out_path)


def get_disaggregated_waste_consumption_for_chp(
    energy_balance, population_intensity, units, carrier_code
):
    waste_consumption = (
        energy_balance.loc[idx["TI_EHG_E", carrier_code, "TJ", :, :]]
        .sum(level=['country', 'year', 'unit'])
        .apply(util.tj_to_twh)  # TJ to TWh as energy units
        .rename_axis(index={"country": "country_code"})
        .rename(util.get_alpha3, level='country_code')  # e.g. GB to GBR
        .rename({"TJ": "twh"}, level='unit')
    )
    return align_and_scale(
        waste_consumption, population_intensity, units
    )


if __name__ == "__main__":
    annual_waste_supply(
        path_to_annual_energy_balances=snakemake.input.energy_balance,
        path_to_population=snakemake.input.population,
        path_to_units=snakemake.input.units,
        renewable_out_path=snakemake.output.renewable,
        total_out_path=snakemake.output.total,
    )
