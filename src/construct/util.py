import pycountry
import pandas as pd
import pytz
import numpy as np

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
    series = series.astype(str).str.extract('(\-*\d+\.*\d*)')[0]
    return pd.to_numeric(series, errors='coerce')


def pj_to_twh(array):
    return array / 3.6


def tj_to_twh(array):
    return pj_to_twh(array) * 1e-3


def ktoe_to_twh(array):
    return array * 1.163e-2


def update_timeseries_timezone(x, country, model_year):
    """
    Shift a generic profile forward/backward in time based on a country's timezone
    """
    if country == 'UK':
        country = 'GB'
    elif country == 'EL':
        country = 'GR'
    tz = pytz.country_timezones[country][0]
    try:
        idx = x.index.tz_localize(tz, nonexistent='shift_forward').tz_convert('UTC')
    except pytz.AmbiguousTimeError as err:
        idx = x.index.tz_localize(
            tz, ambiguous=x.index != err.args[0], nonexistent='shift_forward'
        ).tz_convert('UTC')
    shift = len(idx[idx.year > model_year]) - len(idx[idx.year < model_year])

    x = np.roll(x, shift=shift)

    return x


def read_tdf(filename):
    df = pd.read_csv(filename, header=0)
    tdf = df.set_index([i for i in df.columns[:-1]]).squeeze()
    return tdf


def get_timedelta(model_time, model_year):
    model_timedelta = (pd.to_datetime(model_time[1]) - pd.to_datetime(model_time[0])).days + 1
    if pd.to_datetime(model_year).is_leap_year:
        model_timedelta = model_timedelta / 366
    else:
        model_timedelta = model_timedelta / 365
    return model_timedelta
