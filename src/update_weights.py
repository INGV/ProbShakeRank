import numpy as np
import os
import argparse
import glob
import shutil
from kagan_angle import get_kagan_angle

parser = argparse.ArgumentParser()
parser.add_argument("--ev_id", type=str, required=True, help="Event ID (e.g., 20160903_0000009)")
parser.add_argument("--beta", type=float, required=True, help="β decaying factor in the Kagan's angle formula")
args = parser.parse_args()
ev_id = args.ev_id
beta = args.beta
print("β = ", beta)

EV_DIR = os.path.join(os.getcwd(), f"INPUT_FILES/events/{ev_id}")
event_file = os.path.join(EV_DIR, f"{ev_id}.txt")
print("Event file: ", event_file)

# Read event file
with open(event_file) as f:
    line = f.readline()
    ev = line.strip().split()

Strike_ev = float(ev[4])
Dip_ev = float(ev[5])
Rake_ev = float(ev[6])

# Ensemble file
ensemble_files = glob.glob(os.path.join(EV_DIR, "list_nb0_*.txt"))

kagan_angles = []
for file in ensemble_files:
    with open(file) as f:
        for line in f:
            scen = line.strip().split()

            Strike = float(scen[4])
            Dip = float(scen[5])
            Rake = float(scen[6])

            kagan = get_kagan_angle(
                Strike, Dip, Rake,
                Strike_ev, Dip_ev, Rake_ev
            )

            kagan_angles.append(kagan)

kagan_angles = np.array(kagan_angles)

# No need to normalize here; ProbShakemap normalizes automatically
weights = np.exp(-kagan_angles / beta)

weights_file = "weights.txt"
np.savetxt(weights_file, weights, fmt="%.6f")
print("Updated scenarios weights saved in weights.txt")

# Copy the weights file to 'INPUT_FILES' folder'
source_path = os.path.join(os.getcwd(), weights_file)
destination_path = os.path.join(os.getcwd(), "INPUT_FILES/", weights_file)
shutil.copy(source_path, destination_path)
print("File weights.txt copied to 'INPUT_FILES'")