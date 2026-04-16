#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
from scipy import constants
import math
import argparse
import shutil
import json
import requests
import sys

def dist_lonlat(lon1,lat1,lon2,lat2,coordtype):

    """
    Find dist (km) btw 2 points given their lon,lat
    """

    R = 6371.0

    if coordtype == 'degree':
        lon1 = lon1 * constants.pi / 180
        lat1 = lat1 * constants.pi / 180
        lon2 = lon2 * constants.pi / 180
        lat2 = lat2 * constants.pi / 180

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
    c = 2 * np.arcsin(np.sqrt(a))

    dist_lonlat = R * c

    return dist_lonlat

def convert_accel_units(acceleration, from_, to_='cm/s/s'):  # noqa
    """
    Converts acceleration from/to different units

    :param acceleration: the acceleration (numeric or numpy array)
    :param from_: unit of `acceleration`: string in "g", "m/s/s", "m/s**2",
        "m/s^2", "cm/s/s", "cm/s**2" or "cm/s^2"
    :param to_: new unit of `acceleration`: string in "g", "m/s/s", "m/s**2",
        "m/s^2", "cm/s/s", "cm/s**2" or "cm/s^2". When missing, it defaults
        to "cm/s/s"

    :return: acceleration converted to the given units (by default, 'cm/s/s')
    """
    m_sec_square = ("m/s/s", "m/s**2", "m/s^2")
    cm_sec_square = ("cm/s/s", "cm/s**2", "cm/s^2")
    acceleration = np.asarray(acceleration)
    if from_ == 'g':
        if to_ == 'g':
            return acceleration
        if to_ in m_sec_square:
            return acceleration * g
        if to_ in cm_sec_square:
            return acceleration * (100 * g)
    elif from_ in m_sec_square:
        if to_ == 'g':
            return acceleration / g
        if to_ in m_sec_square:
            return acceleration
        if to_ in cm_sec_square:
            return acceleration * 100
    elif from_ in cm_sec_square:
        if to_ == 'g':
            return acceleration / (100 * g)
        if to_ in m_sec_square:
            return acceleration / 100
        if to_ in cm_sec_square:
            return acceleration

    raise ValueError("Unrecognised time history units. "
                     "Should take either ''g'', ''m/s/s'' or ''cm/s/s''")

#%%

parser = argparse.ArgumentParser()
parser.add_argument("--imt", type=str, required=True, help="Intensity Measure Type (e.g., SA, PGA)")
parser.add_argument("--T", type=str, required=False, default="1_000", help="Period (e.g., 1_000)")
parser.add_argument("--ev_id", type=str, required=True, help="Event ESM ID (e.g., 20160903_0000009)")
parser.add_argument("--fault_multiplier", type=int, required=False, default=1, help="Fault multiplier (default=1)")
parser.add_argument("--max_dist", type=float, required=False, default=20, help="Stations max distance (default=20 km)")

args = parser.parse_args()

imt = args.imt
T = args.T
ev_id = args.ev_id
fault_multiplier = args.fault_multiplier
max_dist = args.max_dist

out_path = os.path.join(os.getcwd(), f"INPUT_FILES/DATA/{ev_id}")
if not os.path.exists(out_path):
    os.makedirs(out_path)

#%%

# Construct the query URL from ESM website
fields = [ 'u_pga', 'v_pga', 'rotd50_pga', f"u_t{T}", f"v_t{T}", f"rotd50_t{T}", 
          'mw', 'network_code', 'station_code', 'st_latitude', 'st_longitude', 
          'event_time', 'ev_latitude', 'ev_longitude', 'ev_depth_km'] 

url_base = "https://esm-db.eu/esmws/flatfile/1/query" 

query = ( f"{url_base}?eventid={ev_id}" f"&include-fields={','.join(fields)}" "&output-format=csv" ) 
print("Query ESM URL:", query) 

response = requests.get(query) 

esm_query_file = os.path.join(out_path, "esm_query_result.csv") 

with open(esm_query_file, "w", encoding='utf-8') as f: 
    f.write(response.text) 
try: 
    df = pd.read_csv(esm_query_file, delimiter=';', low_memory=False) 
    if df.empty: 
        print(f"No ESM data available for event {ev_id}") 
    else: 
        print(f"Loaded {len(df)} records from ESM")
except pd.errors.EmptyDataError: 
    print(f"No ESM data available for event {ev_id} (empty CSV).") 
    sys.exit(1)

distances = dist_lonlat(df['st_longitude'], df['st_latitude'], df['ev_longitude'], df['ev_latitude'], 'degree')
df['distance_epicenter'] = distances

# Get event metadata
coordinates = [
    float(df['ev_longitude'].iloc[0]),
    float(df['ev_latitude'].iloc[0]),
    float(df['ev_depth_km'].iloc[0])
]
time_str = df['event_time'].iloc[0]
Mw = float(df['mw'].iloc[0])

# Derive rupture length using Mai and Beroza (2000)

print("Mw = ", Mw)

seismic_moment = 10 ** (3 / 2 * (Mw + 10.7)) * 10 ** (-7)
fault_length = 10 ** (-5.20 + 0.35 * math.log10(seismic_moment)) # km
print("Fault length = ", round(fault_length,2), "km")

if fault_length < max_dist:
    dmax = max_dist
else:
    dmax = fault_multiplier * fault_length
print("Stations Max Distance = ", dmax, "km")

# Filter data within a multiple of the fault length (or within [user input] km if fault_length < [user input] km)
df = df[df['distance_epicenter'] <= dmax]

# Remove duplicate stations
df = df.drop_duplicates(
    subset=['st_latitude', 'st_longitude'],
    keep='first'
)

# Filter data based on the selected IM 
if imt == 'SA':
    selected_cols = [f"u_t{T}", f"v_t{T}", f"rotd50_t{T}", 'network_code', 'station_code', 'st_latitude', 'st_longitude']
    df = df[selected_cols]
else:
    selected_cols = ['u_pga', 'v_pga', 'rotd50_pga', 'network_code', 'station_code', 'st_latitude', 'st_longitude']
    df = df[selected_cols]

# Convert cm/s² to g

df = df.copy()

g = 9.80665  # m/s^2

if imt == 'PGA':
    df['u_pga'] = convert_accel_units(df['u_pga'],'cm/s/s', 'g')
    df['v_pga'] = convert_accel_units(df['v_pga'],'cm/s/s', 'g')
    df['rotd50_pga'] = convert_accel_units(df['rotd50_pga'],'cm/s/s', 'g')
    print("Num observations = ", len(df['u_pga'] ))
    if len(df['u_pga']) == 0:
        msg = f"Warning: no data within {dmax} km"
        print(msg)
        raise RuntimeError(msg)
else:
    df[f"u_t{T}"] = convert_accel_units(df[f"u_t{T}"],'cm/s/s', 'g')
    df[f"v_t{T}"] = convert_accel_units(df[f"v_t{T}"],'cm/s/s', 'g')
    df[f"rotd50_t{T}"] = convert_accel_units(df[f"rotd50_t{T}"],'cm/s/s', 'g')
    print("Num observations = ", len(df[f"u_t{T}"] ))
    if len(df[f"u_t{T}"] ) == 0:
        msg = f"Warning: no data within {dmax} km"
        print(msg)
        raise RuntimeError(msg)

df = df.sort_values(
    by=["network_code", "station_code"]
).reset_index(drop=True)

# Save to out files

if imt == 'PGA':
    # ROTD50 
    np.savetxt(os.path.join(out_path, f"rotD50_{imt.lower()}.txt"), df[f"rotd50_{imt.lower()}"].to_numpy())

    # GEOMETRIC MEAN 
    gm = np.sqrt(df[f"u_{imt.lower()}"].abs().to_numpy() * df[f"v_{imt.lower()}"].abs().to_numpy())
    np.savetxt(os.path.join(out_path, f"gm_{imt}.txt"), gm)

    # H MAX 
    hmax = np.maximum(df[f"u_{imt.lower()}"].abs().to_numpy(), df[f"v_{imt.lower()}"].abs().to_numpy())
    np.savetxt(os.path.join(out_path, f"hmax_{imt}.txt"), hmax)
else:
    suff = '.'.join([T.split('_')[0], T.split('_')[1][:]])

    # ROTD50 
    np.savetxt(os.path.join(out_path, f"rotd50_{imt}({suff}).txt"), df[f"rotd50_t{T}"].to_numpy())

    # GEOMETRIC MEAN 
    gm = np.sqrt(df[f"u_t{T}"].abs().to_numpy() * df[f"v_t{T}"].abs().to_numpy())
    np.savetxt(os.path.join(out_path, f"gm_{imt}({suff}).txt"), gm)

    # H MAX 
    hmax = np.maximum(df[f"u_t{T}"].abs().to_numpy(), df[f"v_t{T}"].abs().to_numpy())
    np.savetxt(os.path.join(out_path, f"hmax_{imt}({suff}).txt"), hmax)


df['station_id'] = df['network_code'].astype(str) + '_' + df['station_code'].astype(str) 

station_lats = df['st_latitude']
station_lons = df['st_longitude']

#print(df[['station_id', 'st_latitude', 'st_longitude']].head())

pois = np.column_stack((station_lats, station_lons))
np.savetxt(os.path.join(out_path, 'pois.txt'), pois, fmt='%.6f')

# Copy POIS file to INPUT_FILE

pois_file = os.path.join(out_path, 'pois.txt')
dst_folder = "INPUT_FILES"

shutil.copy(pois_file, dst_folder)

# WRITE event_stat.json FILE

eventstatjson_dir = os.path.join("INPUT_FILES/events", ev_id)
os.makedirs(eventstatjson_dir, exist_ok=True)
eventstatjson = os.path.join(eventstatjson_dir, "event_stat.json")

# Default magnitude percentiles
p50 = Mw
p16 = round((Mw - 0.3), 1)
p84 = round((Mw + 0.3), 1)

# Covariance matrix (units: km)
cov_matrix = {
    "XX": 10.0, "XY": 0.0, "XZ": 0.0,
    "YY": 10.0, "YZ": 0.0, "ZZ": 10.0
}

# Construct the .json file
data = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": coordinates
            },
            "properties": {
                "originId": 99999999,
                "version": 1,
                "time": time_str,
                "eventId": ev_id,
                "mag": p50,
                "magType": "Mw",
                "type": "earthquake",
                "place": "Italy",
                "mag_percentiles": {
                    "p16": p16,
                    "p50": p50,
                    "p84": p84
                },
                "cov_matrix": cov_matrix,
                "author": "auto(ESM)"
            }
        }
    ]
}

# Write to file
with open(eventstatjson, "w") as f:
    json.dump(data, f, indent=4)
