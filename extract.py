#!/usr/bin/env python2

import os, calendar, datetime, pandas as pd, numpy as np
from scipy.io import netcdf
from copy import deepcopy
from IPython import embed as IP

class WindSite():
    pass

def day_in_a_year(this_day):
    d0 = datetime.datetime(this_day.year, 1, 1)
    return (this_day - d0).days + 1

def read_netcdf(fnc):
    # Given a WIND toolkit netCDF file, return a dummy object includeing all
    # components
    nc = netcdf.netcdf_file(fnc, 'r')
    s = WindSite()
    for k in nc.variables.keys():
        setattr( s, k, deepcopy(nc.variables[k][:]) )
    nc.close()
    return s

def offset_WIND_data(original, utc=0, forecast=False):
    # Fix the time shift issue and return the correct array.
    # utc is the UTC time zone, defaut to 0, for example, if the site is in
    # Colorado (Time zone: UTC -7), you should set utc = -7.
    # Note the standard time zone should be used, not the daylight saving time.
    # The array original should be the complete data that is from 2007 to 2013.
    if utc is 0:
        return original
    else:
        offset = abs(utc)
    new      = np.array([])
    nndays = 0 # Cumulative number of days
    # Valid year from 2007 to 2012, 6 years in total
    for y in range(2007, 2013):
        year_end  = datetime.datetime(y, 12, 31)
        ndays     = day_in_a_year(year_end)
        if forecast:
            i_start   = nndays*24
            nndays    = nndays + ndays
            i_end     = nndays*24
            this_year = original[i_start: i_end]
            head      = this_year[0: offset]
            tail      = this_year[offset: ]
        else:
            i_start   = nndays*24*12
            nndays    = nndays + ndays
            i_end     = nndays*24*12
            this_year = original[i_start: i_end]
            head      = this_year[0: offset*12]
            tail      = this_year[offset*12: ]
        # print y
        # print i_start, i_end
        # print head.shape[0], tail.shape[0]
        this_year_new = np.append( tail, head )
        new = np.append(new, this_year_new)
    return new

################################################################################

def extract_day(trace, year, month, day, forecast=False):
    # Trace is a 631296 x 1 array, if forecast = True, we will extract forecast
    # data. Default value is actual data (forecast = False).
    # In WIND toolkit, actual data is in 5-min resolution and forecasts are hourly.
    nday = 0
    for y in range(2007, year):
        nday += day_in_a_year( datetime.datetime(y, 12, 31) )
    nday += day_in_a_year( datetime.datetime(year, month, day) )
    if forecast:
        i_start = (nday-1)*24
        i_end   = nday*24
    else:
        i_start = (nday-1)*24*12
        i_end   = nday*24*12
    return trace[i_start: i_end]

def write_csv_file(fnc_actual, fnc_frcst, date_list=None):
    # This script read strings of filename fnc_actual and fnc_crcst, and extract
    # HOURLY AVERAGE wind power, both actual and day-ahead forecast from, of the
    # date specified in date_list and write to a csvfile. If date_list is empty
    # then all date will be extracted, which is from 2007-01-01 to 2012-12-31.
    # Include actual data and day ahead forecasted data.
    directory, fname = os.path.split(fnc_actual)
    str_sid = fname.split('.')[0]

    if date_list:
        pass
    else:
        date_start = datetime.datetime(2007, 1, 1)
        date_end   = datetime.datetime(2012, 12, 31)
        delta = date_end - date_start
        date_list = list()
        for i in range(0, delta.days+1):
            date_list.append(date_start + datetime.timedelta(days=i))

    s_actual     = read_netcdf(fnc_actual)
    power_actual = offset_WIND_data(s_actual.power, -6, False) # Texas is in UTC -6
    s_frcst      = read_netcdf(fnc_frcst)
    power_frcst = dict()
    item_forecasted = [
        'day_ahead_power',
        '4_hour_ahead_power',
        '6_hour_ahead_power',
        'hour_ahead_power',
    ]
    for i in item_forecasted:
        power_frcst[i] = offset_WIND_data(getattr(s_frcst, i), -6, True) # Texas is in UTC -6
    output = {
        'year':        list(),
        'month':       list(),
        'day':         list(),
        'hour':        list(),
        'minute':      list(),
        'actual':      np.array([]),
    }
    p2 = dict() # Forecasted
    for i in item_forecasted:
        output[i] = np.array([])
        p2[i] = np.array([])
    p1 = np.array([]) # Actual

    for d in date_list:
        p1_day = extract_day(power_actual, d.year, d.month, d.day)
        p1_day = np.reshape(p1_day, (12, 24), order='F')
        p1_day_hourly = np.mean(p1_day, axis=0) # Hourly average
        p1 = np.append(p1, p1_day_hourly)
        for i in item_forecasted:
            tmp = extract_day(power_frcst[i], d.year, d.month, d.day, True)
            p2[i] = np.append(p2[i], tmp)
        output['year']   += [d.year]*24
        output['month']  += [d.month]*24
        output['day']    += [d.day]*24
        output['hour']   += range(0, 24)
        output['minute'] +=[0]*24
    output['actual']    = np.append(output['actual'], p1)
    for i in item_forecasted:
        output[i] = np.append(output[i], p2[i])
    df_output = pd.DataFrame(output)
    df_output.to_csv(
        str_sid + '.csv',
        index=False,
        columns=[
            'year',
            'month',
            'day',
            'hour',
            'minute',
            'actual',
            # 'day-ahead',
        ] + item_forecasted,
    )

def extract_texas2000():
    # Extract RT and forecasting wind power of all wind farms from the TEXAS2000
    # system and store them to csv files, in format as specified by Mucun.
    # Note that this function calls write_csv_file(), which extracts hourly data.
    dir_home = os.getcwd()
    dir_save = 'test'
    os.chdir(dir_save)
    path_actual = 'D:\\WINDTOOKIT\\api\\met_data'
    path_frcst  = 'D:\\WINDTOOKIT\\api\\fcst_data'
    df_wind_mapping   = pd.DataFrame.from_csv('C:\\Users\\bxl180002\\git\\FlexibleRampSCUC\\TEXAS2k_B\\wind_generator_data.csv')
    # df_wind_meta_file = pd.DataFrame.from_csv('C:\\Users\\bxl180002\\git\\FlexibleRampSCUC\\TEXAS2k_B\\WIND\\wtk_site_metadata.csv')
    df_wind_meta_file = pd.DataFrame.from_csv('C:\\Users\\bxl180002\\OneDrive\\Tmp_RampWind\\code_before_github\\WIND\\wtk_site_metadata.csv')
    map_wind2site = df_wind_mapping['SITE_ID'].to_dict()
    map_site2dir  = df_wind_meta_file['full_timeseries_directory'].to_dict()
    wind_farms = set()
    s_done = set()
    for w in map_wind2site:
        s = map_wind2site[w]
        d = map_site2dir[s]
        f_actual = os.path.sep.join( [path_actual, str(d), str(s)+'.nc'] )
        f_frcst  = os.path.sep.join( [path_frcst,  str(d), str(s)+'.nc'] )
        if s not in s_done:
            write_csv_file(f_actual, f_frcst)
        s_done.add(s)
    os.chdir(dir_home)

def extract_site_actual(s):
    # This is extracting one element by one element, the basic implementation,
    # we can improve its time consumption by block extraction.
    dir_home = os.getcwd()
    dir_save = 'test'
    path_actual = 'D:\\WINDTOOKIT\\api\\met_data'
    path_frcst  = 'D:\\WINDTOOKIT\\api\\fcst_data'
    os.chdir(dir_save)
    df_wind_meta_file = pd.DataFrame.from_csv('C:\\Users\\bxl180002\\OneDrive\\Tmp_RampWind\\code_before_github\\WIND\\wtk_site_metadata.csv')
    map_site2dir  = df_wind_meta_file['full_timeseries_directory'].to_dict()
    d = map_site2dir[s]
    f_actual = os.path.sep.join( [path_actual, str(d), str(s)+'.nc'] )
    f_frcst  = os.path.sep.join( [path_frcst,  str(d), str(s)+'.nc'] )
    s_actual     = read_netcdf(f_actual)
    power_actual = offset_WIND_data(s_actual.power, utc=-6, forecast=False)
    speed_actual = offset_WIND_data(s_actual.wind_speed, utc=-6, forecast=False)
    output = {
        'year':        list(),
        'month':       list(),
        'day':         list(),
        'hour':        list(),
        'minute':      list(),
        'actual':      list(),
        'speed (m/s)': list(),
        # 'day-ahead':   power_frcst,
    }
    directory, fname = os.path.split(f_actual)
    str_sid = fname.split('.')[0]
    for y in range(2007, 2013):
        for m in range(1, 13):
            mrange = calendar.monthrange(y, m)
            eom = mrange[1] # End of month
            for d in range(1, eom+1):
                for h in range(0, 24): # Let's start from 0 o'clock
                    for m5 in range(0, 12): # 1 hour has 12 5-min interval
                        d0 = datetime.datetime(2007, 1, 1)
                        d_now = datetime.datetime(y, m, d)
                        i = ((d_now-d0).days*24 + h)*6 + m5
                        output['year'].append(y)
                        output['month'].append(m)
                        output['day'].append(d)
                        output['hour'].append(h)
                        output['minute'].append(5*m5)
                        output['actual'].append(power_actual[i])
                        output['speed (m/s)'].append(speed_actual[i])
    df_output = pd.DataFrame(output)
    df_output.to_csv(
        str_sid + '.csv',
        index=False,
        columns=[
            'year',
            'month',
            'day',
            'hour',
            'minute',
            'actual',
            'speed (m/s)',
            # 'day-ahead',
        ],
    )
    os.chdir(dir_home)


if __name__ == "__main__":
    extract_texas2000()
    # for s in [20481, 21055, 22235]:
    #     extract_site_actual(s)

