#!/bin/bash
set -e

# ---------- User Input ----------
IMT="PGA"
T="1_000" # FOR SA ONLY (will be ignored for other IMT)
EV_ID="IT-2012-0011"
FAULT_MULTIPLIER=1
STATIONS_MAX_DIST=20
BETA=5
NUM_GMPES=20
PROC=8

HOME_DIR=$(pwd)
echo "HOME_DIR: $HOME_DIR"


# ---------- Auto-set IMT_MIN and IMT_MAX based on event and IMT ----------

if [[ "$IMT" == "PGA" ]]; then
    IMT_MIN=0.01
    IMT_MAX=1
elif [[ "$IMT" == "SA" && "$T" == "1_000" ]]; then
    IMT_MIN=0.1
    IMT_MAX=1
fi


# ----------  Load Conda  ----------

source ~/miniconda3/etc/profile.d/conda.sh

conda activate probshakerank


# ---------- Set ev-related files ----------
INPUT_DIR="INPUT_FILES"

# Remove older pois file
if [ -f "$INPUT_DIR/pois.txt" ]; then
    rm "$INPUT_DIR/pois.txt"
else
    echo "No POIs file"
fi

DATA_DIR="INPUT_FILES/DATA/${EV_ID}"
cp "$DATA_DIR/pois.txt" "$INPUT_DIR/pois.txt"
echo "POIs file copied to $INPUT_DIR"

EVENT_DIR="INPUT_FILES/events/${EV_ID}"
ENSEMBLE_DIR="$HOME_DIR/INPUT_FILES/ENSEMBLE"

ENSEMBLE_FILE=$(echo "$EVENT_DIR"/*Ev_${EV_ID}_*.txt)
echo "Ensemble file: $ENSEMBLE_FILE"

# Cleaning ENSEMBLE_DIR before copying new file
rm -f "$ENSEMBLE_DIR"/* || true

cp "$ENSEMBLE_FILE" "$ENSEMBLE_DIR/"
echo "Ensemble file copied to $ENSEMBLE_DIR"


# ---------- Updating scenarios weights ----------
echo "Updating scenarios weights"

python src/update_weights.py \
  --ev_id "$EV_ID" \
  --beta "$BETA" || {
    echo "Failed: update_weights.py"
    exit 1
}



# ----------  Run ProbShakemap ----------
echo  "------------------------------"
echo  "  Running ProbShakemap with updated scenarios weights..."
echo  "------------------------------"

if [[ "$IMT" == "SA" ]]; then
    T_FLOAT="${T/_/.}"  # Convert 1_000 -> 1.000
    IMT_STR="SA(${T_FLOAT})"
else
    IMT_STR="$IMT"
fi

echo "Running with IMT: $IMT_STR"

python src/probshakemap.py \
  --ev_id "$EV_ID" \
  --imt "$IMT_STR" \
  --num_processes "$PROC" \
  --pois_file pois.txt \
  --numGMPEsRealizations "$NUM_GMPES" \
  --imt_min "$IMT_MIN" \
  --fileScenariosWeights "weights.txt" \
  --imt_max "$IMT_MAX" || {
    echo "Failed: probshakemap.py"
    exit 1
}

conda deactivate
