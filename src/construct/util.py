
import pycountry
import pandas as pd

def get_alpha2(country, eurostat=True):
    if country in ["United Kingdom", "GB", "GBR"] and eurostat is True:
        return "UK"
    elif country in ["Greece", "GR", "GRC"] and eurostat is True:
        return "EL"
    else:
        return pycountry.countries.lookup(country).alpha_2


def get_alpha3(country):
    if country == "UK":
        country = "GB"
    elif country == "EL":
        country = "GR"
    return pycountry.countries.lookup(country).alpha_3


def to_numeric(series):
    return pd.to_numeric(series, errors='coerce')


def tj_to_twh(array):
    array /= 3600
    return array