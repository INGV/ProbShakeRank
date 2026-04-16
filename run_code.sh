#!/bin/bash
set -e

# ---------- User Input ----------
IMT="SA"
T="1_000" # FOR SA ONLY (will be ignored for other IMT)
EV_ID="IT-2012-0011"
FAULT_MULTIPLIER=1
STATIONS_MAX_DIST=20
N_SCENS=1000
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


# ---------- Step 1: Extract ESM data ----------

# Extract data 

DATA_DIR="INPUT_FILES/DATA/${EV_ID}"
mkdir -p "$DATA_DIR"
EVENT_DIR="INPUT_FILES/events/${EV_ID}"
mkdir -p "$EVENT_DIR"

echo  "------------------------------"
echo "Step 1: Extracting ESM data..."
echo  "------------------------------"

conda activate probshakerank
python src/Get_Event_Data.py --imt "$IMT" --T "$T" --ev_id "$EV_ID" --fault_multiplier "$FAULT_MULTIPLIER" --max_dist "$STATIONS_MAX_DIST" || {
    echo "Step 1 failed: Get_Event_Data.py"
    exit 1
}
conda deactivate


# ---------- Step 2: Run SeisEnsMan  ----------
#(only if ensemble for the current event doesn't already exist)

DST_EVENT_DIR="src/SeisEnsManV2/input"

if [ ! -f "$EVENT_DIR/event_stat.json" ]; then
    echo "ERROR: Missing $EVENT_DIR/event_stat.json"
    exit 1
fi

cp "$EVENT_DIR/event_stat.json" "$DST_EVENT_DIR/event_stat.json"
echo "$EVENT_DIR/event_stat.json copied in $DST_EVENT_DIR"

# Generate ensemble for current event only if it doesn't exist, otherwise use existing file

ENSEMBLE_DIR="$HOME_DIR/INPUT_FILES/ENSEMBLE"
ENSEMBLE_EXISTS=$(find "$EVENT_DIR" -type f -name "*Ev_${EV_ID}_*.txt" | head -n 1)

# Cleaning ENSEMBLE_DIR before copying new ensemble
rm -f "$ENSEMBLE_DIR"/* || true

if [ -n "$ENSEMBLE_EXISTS" ]; then
    echo "Ensemble file already exists for event $EV_ID:"
    echo "$ENSEMBLE_EXISTS"
    
    cp "$ENSEMBLE_EXISTS" "$ENSEMBLE_DIR"
    echo "Ensemble file copied to $ENSEMBLE_DIR"
    
else
    echo  "-----------------------------"
    echo "Step 2: Running SeisEnsMan..."
    echo  "-----------------------------"
    source /usr/local/bin/SeisEnsMan/bin/activate   
    cd src/SeisEnsManV2 || exit

    mainFolder=$(pwd)
    python backup.py
    python run_ens.py --ev_id "$EV_ID" --cfg "$mainFolder/input/main.config" --event "$mainFolder/input/event_stat.json" --nb_scen "$N_SCENS" || {
        echo "Step 2 failed: run_ens.py"
        exit 1
    }
    
    cd "$HOME_DIR" || exit

    # OPTIONAL: --angles Strike Dip Rake 
    # Example: norcia 151 47 -89

    cp "$ENSEMBLE_DIR"/*Ev_"$EV_ID"_*.txt "$EVENT_DIR"/ || true
    echo "A copy of the ensemble file was saved in $EVENT_DIR"

    deactivate   
fi


# ---------- Step 3: Run ProbShakemap ----------
echo  "------------------------------"
echo  "Step 3: Running ProbShakemap..."
echo  "------------------------------"

conda activate probshakerank

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
  --imt_max "$IMT_MAX" || {
    echo "Step 3 failed: probshakemap.py"
    exit 1
}

conda deactivate
