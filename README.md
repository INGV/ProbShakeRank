# ProbShakeRank

**Probabilistic GMM ranking with explicit source uncertainty propagation**

`ProbShakeRank` is a Python workflow for ranking Ground Motion Models (GMMs) against observed seismic records (PGA, PGV, SA). Unlike traditional ranking approaches, it explicitly quantifies and propagates earthquake source uncertainty (magnitude, hypocenter location, strike/dip/rake) into the ranking process.

This enables dynamic updating of GMM performance as source characterization improves, making `ProbShakeRank` suitable for quasi-real-time applications.

The workflow is built on top of [ProbShakemap](https://github.com/INGV/ProbShakemap) (Stallone et al., 2025) and extends it with:

1. Automatic retrieval of event information from the Engineering Strong Motion (ESM) database (Luzi et al., 2016)
2. GMM performance evaluation using multiple ranking metrics:
   - Log-likelihood score (LLH), uni- and multi-variate
   - Parimutuel Gambling Score (PGS)
   - Akaike Information Criterion (AIC)
   - Bayesian Information Criterion (BIC)
3. An interactive dashboard for exploring ranking results

> Associated article: *"Dynamic GMM ranking accounting for source uncertainty"* (Stallone et al., coming soon)

---

## Table of Contents

- [Workflow](#workflow)
- [Installation](#installation)
- [Input Configuration](#input-configuration)
- [Running the Code](#running-the-code)
  - [Dynamic ranking update with scenario weights](#dynamic-ranking-update-with-scenario-weights)
- [Output Structure](#output-structure)
- [Dashboard](#dashboard)
- [Dependencies](#dependencies)
- [Citation](#citation)

---

## Workflow

```
1. Download strong-motion data from the ESM database
2. Generate an ensemble of rupture scenarios with SeisEnsMan
3. Run ProbShakemap to compute predictive ground-motion distributions at recording stations
4. Rank GMMs using LLH, PGS, multivariate LLH, AIC, and BIC
5. Launch interactive dashboard
```
See `run_code.sh` for a workflow example.

---

## Installation

##### 1. Clone the repository

```bash
git clone https://github.com/INGV/ProbShakeRank.git
cd ProbShakeRank
```

##### 2. Set up the Conda environment

```bash
conda env create -f probshakerank_environment.yml
conda activate probshakerank
```

##### 3. Set up SeisEnsMan

(skip this step if source scenarios are generated independently)

Follow the installation instructions in the [ProbShakemap repository](https://github.com/INGV/ProbShakemap). 

---

## Input Configuration

### `run_code.sh` — Workflow input parameters

| Parameter           | Description                           | Example           |
| ------------------- | ------------------------------------- | ----------------- |
| IMT                 | Intensity measure type                | PGA, PGV, SA      |
| T                   | SA period (ignored for PGA/PGV)       | 1\_000 (SA 1.0 s) |
| EV\_ID              | ESM event ID                          | IT-2012-0011      |
| FAULT\_MULTIPLIER   | Station distance scaling factor       | 1                 |
| STATIONS\_MAX\_DIST | Max station distance (km)             | 20                |
| N\_SCENS            | Number of rupture scenarios           | 1000              |
| NUM\_GMPES          | GMF realizations per scenario per GMM | 20                |
| PROC                | CPU cores for parallelization         | 8                 |


### `INPUT_FILES/gmpes.conf` — GMM configuration

Defines the sets of GMMs to rank and their associated horizontal components. Each GMM is identified by an acronym and mapped to its  `OpenQuake` GSIM class name.

### `input_file.txt` — ProbShakemap configuration

| Parameter | Description |
|---|---|
| `TectonicRegionType` | Tectonic regime (e.g., `Active Shallow Crust`) |
| `Magnitude_Scaling_Relationship` | Magnitude scaling relation (e.g., `WC1994`) |
| `Rupture_aratio` | Rupture aspect ratio |
| `Vs30file` | Path to Vs30 grid file (optional; defaults to 760 m/s) |
| `CorrelationModel` | Spatial correlation model (e.g., `JB2009`) |
| `CrosscorrModel` | Cross-correlation model (e.g., `GodaAtkinson2009`) |
| `truncation_level` | Truncation level for ground-motion sampling |
| `seed` | Random seed for reproducibility |

See [ProbShakemap repository](https://github.com/INGV/ProbShakemap) for more details.


---

## Running the Code

To run `ProbShakeRank` workflow:

```bash
bash run_code.sh
```

### Dynamic ranking update with scenario weights
 
`ProbShakeRank` supports dynamic updates of GMM rankings by re-weighting rupture scenarios based on updated source information.

Scenario weights are computed using Kagan angle similarity and a decay factor `BETA`, following Cordrie et al. (2025) and implemented in [pyPTF_data_update](https://github.com/louisecordrie/pyPTF_data_update).

To run weighted ranking:

```bash
bash run_code_weights.sh
```

This skips data retrieval and ensemble generation, and directly updates:

- scenario weights (`weights.txt`)
- ground-motion distributions via `ProbShakemap` (`--fileScenariosWeights`)
- GMM ranking

| Parameter | Description                                 | Example |
| --------- | ------------------------------------------- | ------- |
| BETA      | Decay factor controlling scenario weighting | 5       |

 
---

## Output Structure

All outputs are written to `OUTPUT/<EV_ID>/`:

```
OUTPUT/
└── <EV_ID>/
    ├── metadata.txt                         # Event Mw, time, coordinates
    ├── RANK/
    │   ├── LLH_Score_<IMT>.txt              # Per-GMM LLH scores
    │   ├── Gambling_Score_<IMT>.txt         # Per-GMM PGS scores
    │   ├── MultiIMs_Ranking.txt             # Multivariate LLH, AIC, BIC
    │   └── asll_storage.npy                 # Accumulated per-site log-likelihoods (multi-IM)
    └── OUTPUT_<GMM>_<IMT>/
        ├── <IMT>_<GMM>_mean_total_stdev_scens_OQ.npy
        ├── <IMT>_<GMM>_mean_total_stdev_scens.npy
        ├── npyFiles/
        │   └── <IMT>_<GMM>_vector_stat.npy  # Mean and std dev at all POIs
        ├── STATISTICS/
        │   ├── Stat_Mean.png                 # Spatial map of predicted mean
        │   └── Stat_ST_DEV.png              # Spatial map of predicted std dev
        ├── RANK_FIGURES/
        │   ├── Normalized_Residuals.png      # Per-site residual map
        │   ├── POI_LLH.png                   # Per-site LLH contribution map
        │   └── POI_Gambling.png              # Per-site PGS contribution map
        └── LOGS/
            └── setup_<timestamp>.log
```

---

## Dashboard

After the ranking step completes, an interactive dashboard is launched automatically, which displays:

- Predictive distribution statistics (mean and standard deviation) at each recording station, per GMM and per IMT
- Normalized residual maps
- Per-station contributions to LLH and PGS
- Summary ranking tables (LLH, PGS, AIC, BIC)

To relaunch the dashboard manually at any moment:

```bash
conda activate probshakerank
streamlit run src/dashboard.py
```

---

## Dependencies

- [OpenQuake hazardlib](https://github.com/gem/oq-engine) 
- [ProbShakemap](https://github.com/INGV/ProbShakemap) 

---

## Citation

If you use `ProbShakeRank` in your research, please cite:

```
(Stallone et al., coming soon)
Dynamic GMM ranking accounting for source uncertainty
```

and the underlying `ProbShakemap` toolbox:

```
Stallone, A., et al. (2025). ProbShakemap: A probabilistic ground-motion prediction
toolbox for urgent computing applications.
```

---

## Contact

Angela Stallone\
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Bologna, Italy\
[angela.stallone@ingv.it](mailto\:angela.stallone@ingv.it)

---

## License

This project is licensed under the terms described in `LICENSE.txt`.
