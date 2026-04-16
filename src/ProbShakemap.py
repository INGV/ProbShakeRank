import tools
import os
import numpy as np
import argparse
import time
import shutil
import logging
from datetime import datetime
import config
from configobj import ConfigObj, ConfigObjError
import sys
from openquake.hazardlib import valid
import streamlit as st
import subprocess
from numpy import save
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colormaps
from openquake.hazardlib import imt

start_time = time.time()

def setup_logging(log_dir):
    current_time = datetime.now().strftime("%Y%m%d%H%M%S")
    log_file = os.path.join(log_dir, f"setup_{current_time}.log")

    # Reset root logger by removing all handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def get_gmpes_and_hcomp(imt):

    conf_filename = os.path.join(os.getcwd(), "INPUT_FILES/gmpes.conf")
    try:
        config_gmpes = ConfigObj(conf_filename)
    except ConfigObjError as e:
        print(f"Error loading config file {conf_filename}: {e}")
        raise

    gmpe_sets = config_gmpes['gmpe_sets']
    gmpes_modules = config_gmpes['gmpe_modules']

    print("GMPEs Set(s) = ", gmpe_sets.keys())

    GMPEs_To_Rank = []
    for set in gmpe_sets.values():
        gmpes = set['gmpes']
        if isinstance(gmpes, list):
            GMPEs_To_Rank.extend(gmpes)
        else:
            GMPEs_To_Rank.append(str(gmpes))

    # Get Horiz Comp 
    hcomps_list = []
    for set_data in gmpe_sets.values():
        hcomp = set_data['hcomp']
        if isinstance(hcomp, list):
            hcomps_list.extend(hcomp)
        else:
            hcomps_list.append(str(hcomp))

    # Map acronyms to GMPE names and instantiate GSIM objects
    GMPEs_Names = {key: item[0] for key, item in gmpes_modules.items() if key in GMPEs_To_Rank}
    gmpes = {acronym: valid.gsim(name) for acronym, name in GMPEs_Names.items()} # OpenQuake equivalent of getattr

    # Map horizontal components to acronyms
    Hcomps = dict(zip(GMPEs_Names.keys(), hcomps_list))

    # Filter GMPEs with the selected IMT available
    gmpes_ok = {}
    for acronym, gmpe in gmpes.items():
        list_of_imts = ', '.join([imt.__name__ for imt in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES])   
        if imt in ['PGA', 'PGV']:
            if imt in list_of_imts:
                gmpes_ok[acronym] = gmpe
        else:
            # SA
            if imt[:2] in list_of_imts:
                # Check if chosen SA period is supported
                period = float(imt[3:-1])
                # 1. Try to get COEFFS_* from the subclass
                coeffs_attr = next(
                    (attr for attr in gmpe.__class__.__dict__ if attr.startswith("COEFFS")), None)
                # 2. If not found, look in the base class
                if not coeffs_attr:
                    coeffs_attr = next(
                        (attr for attr in dir(gmpe) if attr.startswith("COEFFS")), None)
                if coeffs_attr:
                    coeffs = getattr(gmpe, coeffs_attr)
                    if hasattr(coeffs, "sa_coeffs"):
                        supported_periods = [sa.period for sa in coeffs.sa_coeffs]
                        #print(f"Targeted COEFFS Table: {coeffs_attr}, belongs to: {gmpe.__class__.__name__}")
                        #print(f"Supported SA Periods for {gmpe.__class__.__name__}: {supported_periods}\n")
                if period in supported_periods:
                    gmpes_ok[acronym] = gmpe

    if len(gmpes) != 0:
        print("GMPEs with the requested IMT available = ", [item for _, item in gmpes.items()])     
    else:
        if imt in ['PGA', 'PGV']:
            print('No GMPEs with the requested IMT available')   
            sys.exit() 
        else:
            #SA
            if imt[:2] in list_of_imts:
                print("Chosen SA period not supported. Supported periods = ", supported_periods)
                sys.exit() 
            else:
                print('No GMPEs with the requested IMT available')   
                sys.exit() 

    # Filter gmpe_names and hcomps 
    GMPEs_Names = {k: v for k, v in GMPEs_Names.items() if k in gmpes_ok}
    Hcomps = {k: v for k, v in Hcomps.items() if k in gmpes_ok}

    return gmpes_ok, GMPEs_Names, Hcomps

def write_event_metadata(out_dir_ev, magnitude, time_str, LatEvent, LonEvent):
    metadata_path = os.path.join(out_dir_ev, "metadata.txt")
    with open(metadata_path, "w") as f:
        f.write(f"Mw: {magnitude:.1f}\n")
        f.write(f"Time: {time_str}\n")
        f.write(f"LatEvent: {LatEvent:.2f}\n")
        f.write(f"LonEvent: {LonEvent:.2f}\n")
    print(f"metadata.txt written at: {metadata_path}")

def Gambling_Score(P):
    """
    Compute gambling score
    Source: "A parimutuel gambling perspective to compare probabilistic seismicity forecasts", Zechar & Zhaung (2014)
    """
    n_models, n_pois = P.shape
    
    denom = np.sum(P, axis=0, keepdims=True)  # shape (1, n_pois)
    delta_R = -1 + n_models * (P / denom)  # shape: (n_gmpes × n_pois)

    R = np.sum(delta_R, axis=1)
    
    return R, delta_R

def Gambling_plot(pois_file, out_dir, delta_R_GMPE, Lon_Event, Lat_Event):
    
    """
    Plot point-gambling score (i.e. gambling score at each POI)
    """

    POIs_lat, POIs_lon, _, _ = tools.get_pois(pois_file)
    POIs_lat = np.array(POIs_lat)
    POIs_lon = np.array(POIs_lon)

    out_path = os.path.join(os.getcwd(), f"{out_dir}/RANK_FIGURES")
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    dim_point = 40

    min_lon, max_lon = np.min(POIs_lon), np.max(POIs_lon)
    min_lat, max_lat = np.min(POIs_lat), np.max(POIs_lat)
    lon_span = max_lon - min_lon
    lat_span = max_lat - min_lat
    buffer_lon = max(0.05 * lon_span, 0.1)  
    buffer_lat = max(0.05 * lat_span, 0.1)

    xlim_min = max(-180, min_lon - buffer_lon)
    xlim_max = min(180, max_lon + buffer_lon)
    ylim_min = max(-90, min_lat - buffer_lat)
    ylim_max = min(90, max_lat + buffer_lat)

    fig = plt.figure(figsize=(9, 6))

    latitudes = np.arange(-90, 91, 2)
    longitudes = np.arange(-180, 181, 2)
    cm = colormaps['viridis']

    m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
            llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i')
    try:
        m.drawcoastlines()
    except Exception:
        pass 

    m.drawparallels(latitudes, labels=[1,0,0,0], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])
    m.drawmeridians(longitudes, labels=[0,0,0,1], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])

    m.drawmapboundary(linewidth=2, color='black', fill_color='white')

    x, y = m(POIs_lon, POIs_lat)
    x_ev, y_ev = m(Lon_Event, Lat_Event)
    
    vmin = np.min(delta_R_GMPE)
    vmax = np.max(delta_R_GMPE)
    sc = m.scatter(x, y, c=delta_R_GMPE, s=dim_point, cmap = cm, vmin=vmin, vmax=vmax, edgecolor='k', linewidths=0.8)
    sc2 = m.scatter(x_ev, y_ev, s=100, c='k', marker='*', edgecolor='k', linewidths=0.3)

    plt.title('POI-Gambling')
    cbar = plt.colorbar(sc)
    #cbar.set_label()
    #plt.show()

    figname = out_path + '/' + 'POI_Gambling.png'
    fig.savefig(figname, bbox_inches="tight", dpi=300, pad_inches=0.1)
    plt.close(fig)


def run_main(gmpe, EnsembleSize):

    """
    Runs tools.Main()
    """

    # Set the random seed for reproducibility (must be the same as in tools.py)
    config_dict = config.load_config('input_file.txt')
    seed = config_dict['seed']
    np.random.seed(seed)

    ProbCalc = tools.Main(gmpe, args.imt, args.pois_file,
                        args.numGMPEsRealizations, args.num_processes, EnsembleSize)
    
    prob_output = ProbCalc.run_prob_analysis()

    return prob_output

if __name__ == '__main__':
    # Define command-line arguments
    parser = argparse.ArgumentParser(description='ProbShakemap Toolbox')
    input_params = parser.add_argument_group('input params')
    input_params.add_argument('--ev_id', help='ID of the event (EMSC).')
    input_params.add_argument('--imt', help='Intensity measure type (IMT). Possible choices: PGA, PGV, SA.')
    input_params.add_argument('--numGMPEsRealizations', type=int, help='Total number of GMPEs random samples')
    input_params.add_argument('--num_processes', type=int, default=1, help='Number of CPU cores for code parallelization')
    input_params.add_argument('--imt_min', type=float, help='Minimum value for the selected IMT (for plot only)')
    input_params.add_argument('--imt_max', type=float, help='Maximum value for the selected IMT (for plot only)')
    input_params.add_argument('--pois_file', help='Filename with latitude and longitude of POIs')
    input_params.add_argument('--vector_npy', action='store_true', default=False, help='Store ground motion distributions at all POIs (vector.npy)')
    input_params.add_argument('--fileScenariosWeights', default="", help='File with scenarios weights')

    # Parse command-line arguments
    args = parser.parse_args()

    # Get params
    params = tools.get_params(args.ev_id)
    Lat_Event = float(params["Lat_Event"])
    Lon_Event = float(params["Lon_Event"])
    event_dir = params["event_dir"]
    EnsembleSize = params['Ensemble_Size']
    Time = params['Time']
    Mw = params['Mw']

    # Set OUTPUT base dir
    ev_id = os.path.basename(os.path.normpath(event_dir))
    OUT_DIR = os.path.join(os.getcwd(), f"OUTPUT/{ev_id}")
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    # Set up output dir for ranking
    OUT_DIR_RANK = os.path.join(OUT_DIR, "RANK")
    if not os.path.exists(OUT_DIR_RANK):
        os.makedirs(OUT_DIR_RANK)

    write_event_metadata(OUT_DIR, Mw, Time, Lat_Event, Lon_Event)

    # Create/Load container gmpe_name: [asll arrays from different IMTs]
    ASLL_STORAGE_FILE = os.path.join(OUT_DIR_RANK, "asll_storage.npy")
    try:
        asll_storage = np.load(ASLL_STORAGE_FILE, allow_pickle=True).item()
    except FileNotFoundError:
        asll_storage = {}

    # Scenarios file
    listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
    scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]

    ###### GET GMPES #####

    GMPES, GMPEs_Names, Hcomps = get_gmpes_and_hcomp(args.imt)
    print(f"GMPEs to Rank = {GMPES}" )

    llh_file = os.path.join(OUT_DIR_RANK, f"LLH_Score_{args.imt}.txt")
    gambling_file = os.path.join(OUT_DIR_RANK, f"Gambling_Score_{args.imt}.txt")
    multi_rank_file = os.path.join(OUT_DIR_RANK, "MultiIMs_Ranking.txt")

    llh_scores = []
    n_models = len(GMPES)
    _, _, _, n_pois = tools.get_pois(args.pois_file)
    P = np.zeros((n_models, n_pois))
    for j, (key, gmpe) in enumerate(GMPES.items()):
        print(f"--------- GMPE: {GMPEs_Names[key]} ------------")

        # Get file with ground-motions observed at the POIs
        hcomp = Hcomps[key]
        data_file = os.path.join(f"./INPUT_FILES/DATA/{ev_id}", f"{hcomp}_{args.imt}.txt")    
        print(f"--------- Data file: {hcomp}_{args.imt}.txt ------------")

        # Set up output dir for current GMPE-IMT combo
        out_dir = os.path.join(OUT_DIR, f"OUTPUT_{key}_{args.imt}")
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # Set up log file
        outfile_dir = f"{out_dir}/LOGS/"
        if not os.path.exists(outfile_dir):
                    os.makedirs(outfile_dir)
        
        setup_logging(outfile_dir)
        logging.info(f"GMPE: {GMPEs_Names[key]}")    
        logging.info(f"Intensity Measure: {args.imt}")
        logging.info(f"Min {args.imt}: {args.imt_min}")
        logging.info(f"Max {args.imt}: {args.imt_max}") 
        logging.info(f"Scenarios Ensemble File: {scenarios_file[0]}") 
        logging.info(f"POIs File: {args.pois_file}")
        logging.info(f"Num GMPEsRealizations: {args.numGMPEsRealizations}")

        if args.imt is None:
            raise TypeError("Missing required argument 'imt'")
        if args.pois_file is None:
            raise TypeError("Missing required argument 'pois_file'")
        if args.numGMPEsRealizations is None:
            raise TypeError("Missing required argument 'numGMPEsRealizations'")

        prob_output = run_main(gmpe, EnsembleSize)
        SiteGmf = prob_output["SiteGmf"]
        # Mean and ST DEV for a given GMPE from OpenQuake, through all scenarios
        # (shape: (k, 4, 0, 0, N) )
        mean_total_stdev_scens_OQ = prob_output["mean_total_stdev_scens_OQ"]
        save(out_dir + '/' + f"{args.imt}_{key}_" + 'mean_total_stdev_scens_OQ.npy', mean_total_stdev_scens_OQ)

        GetStatistics = tools.GetStatistics(out_dir, SiteGmf, EnsembleSize, Lon_Event, Lat_Event, args.numGMPEsRealizations, event_dir, args.imt, 
                                            args.imt_min, args.imt_max, args.pois_file, args.num_processes, args.vector_npy, args.fileScenariosWeights)
        GetStatistics.save_statistics()
        GetStatistics.plot_statistics()
        stats = GetStatistics.calc_statistics()
        mean_and_stdev = stats['mean_and_stdev'] 

        # Mean and ST DEV from ProbShakemap predictive distibutions
        mean_total_stdev_scens = [mean_and_stdev['Mean'], mean_and_stdev['ST_DEV']]
        save(out_dir + '/' + f"{args.imt}_{key}_" + 'mean_total_stdev_scens.npy', mean_total_stdev_scens)

        # Generate plot of normalized residuals
        tools.plot_normalized_residuals(args.pois_file, out_dir, mean_total_stdev_scens, data_file, Lon_Event, Lat_Event)

        # ------ LLH Score ------
        LLH = tools.LLH_Score(data_file, mean_total_stdev_scens, args.pois_file, out_dir, Lon_Event, Lat_Event)
        asll, llh_score = LLH.calculate_llh()
        
        llh_scores.append(llh_score)

        # Generate plot of per-site contribution of POI i to GMPE j’s LLH score
        LLH.LLH_plot(asll)

        # ----- Collect asll for multivariate LLH ------
        gmpe_name = GMPEs_Names[key]
        if gmpe_name not in asll_storage:
            asll_storage[gmpe_name] = []
        asll_storage[gmpe_name].append(asll)

        # ------ Calculate probs for Gambling Score -----
        p_j = tools.compute_GMPE_probabilities(mean_total_stdev_scens, data_file)
        P[j, :] = p_j

    np.save(ASLL_STORAGE_FILE, asll_storage)

    # Save LLH Score to file
    with open(llh_file, "w") as f:
        f.write("{} {}\n".format('GMPE', 'LLH_Score'))
        for gmpe, llh_score in zip(GMPEs_Names.values(), llh_scores):
            f.write("{} {:.3f}\n".format(gmpe, llh_score))  

    # Calculate multivariate ranking
    
    llh_all_by_gmpe = {}
    aic_all_by_gmpe = {}
    bic_all_by_gmpe = {}
    asll_storage = np.load(ASLL_STORAGE_FILE, allow_pickle=True).item()
    for j, (key, gmpe) in enumerate(GMPES.items()):
        gmpe_name = GMPEs_Names[key]  
        asll_list = asll_storage.get(gmpe_name, [])

        # Stack all IMTs together (as in oq-mbtk)
        all_asll = np.hstack(asll_list)
        llh_all = -(1.0 / len(all_asll)) * np.sum(all_asll)
        llh_all_by_gmpe[gmpe_name] = llh_all

        # ------ AIC & BIC ----

        if args.imt.upper() == 'PGA':
            imt_instance = imt.PGA()
        elif args.imt.upper() == 'PGV':
            imt_instance = imt.PGV()
        elif args.imt.upper().startswith('SA'):
            period = float(args.imt.split('(')[1].strip(')')) if '(' in args.imt else 0.3
            imt_instance = imt.SA(period)
        else:
            raise ValueError(f"Unknown IMT: {args.imt}")

        # Convert from average log-likelihood in log2 to total log-likelihood in natural log
        log_likelihood = -llh_all * len(all_asll) * np.log(2)
        if args.imt in ['PGA', 'PGV']:
            coeff = gmpe.COEFFS.non_sa_coeffs[imt_instance]  
        else:
            coeff = gmpe.COEFFS.sa_coeffs[imt_instance]
        k = len(coeff)

        aic_all_by_gmpe[gmpe_name] = 2*k - 2*log_likelihood
        bic_all_by_gmpe[gmpe_name] = k*np.log(len(all_asll)) - 2*log_likelihood

    # Save multivariate rankings to file
    with open(multi_rank_file, "w") as f:
        f.write("GMPE Multivariate_LLH AIC BIC\n")
        for gmpe in llh_all_by_gmpe:
            f.write(f"{gmpe} {llh_all_by_gmpe[gmpe]:.3f} "
                    f"{aic_all_by_gmpe[gmpe]:.3f} {bic_all_by_gmpe[gmpe]:.3f}\n")

    # Calculate Gambling Score
    R, delta_R = Gambling_Score(P)
    # Generate plot of per-site contribution of POI i to GMPE j’s Gambling score
    for j, (key, gmpe) in enumerate(GMPES.items()):
        out_dir = os.path.join(OUT_DIR, f"OUTPUT_{key}_{args.imt}")
        Gambling_plot(args.pois_file, out_dir, delta_R[j,:], Lon_Event, Lat_Event)

    # Save Gambling Score to file
    with open(gambling_file, "w") as f:
        f.write("{} {}\n".format('GMPE', 'Gambling_Score'))
        for gmpe, score in zip(GMPEs_Names.values(), R):
            f.write("{} {:.3f}\n".format(gmpe, score))

    print("********* DONE! *******")
    print("--- %s seconds ---" % (time.time() - start_time))

    print("STARTING GUI...")
    
    streamlit_app = "src/dashboard.py"  
    subprocess.run(["streamlit", "run", streamlit_app])



           