import pandas as pd
import geopandas as gpd

import util
from annual_subnational_demand import align_and_scale, get_population_intensity

idx = pd.IndexSlice


def annual_waste_supply(
    path_to_annual_energy_balances, path_to_population, path_to_units, out_path
):
    """
    Get total waste
    """
    energy_balance = util.read_tdf(path_to_annual_energy_balances)
    waste_consumption = (
        energy_balance.loc[idx["TI_EHG_E", ["W6210", "W6100_6220"], "TJ", :, :]]
        .sum(level=['country', 'year', 'unit'])
        .apply(util.tj_to_twh)  # TJ to TWh as energy units
        .rename_axis(index={"country": "country_code"})
        .rename(util.get_alpha3, level='country_code')  # e.g. GB to GBR
        .rename({"TJ": "twh"}, level='unit')
    )
    units = gpd.read_file(path_to_units)
    population_intensity = get_population_intensity(path_to_population, units)

    disaggregated_waste_supply = align_and_scale(
        waste_consumption, population_intensity, units
    )
    disaggregated_waste_supply.to_csv(out_path)


if __name__ == "__main__":
    annual_waste_supply(
        path_to_annual_energy_balances=snakemake.input.energy_balance,
        path_to_population=snakemake.input.population,
        path_to_units=snakemake.input.units,
        out_path=snakemake.output[0]
    )
