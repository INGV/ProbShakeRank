from openquake.hazardlib.source import BaseRupture
from openquake.hazardlib.geo import Point, surface
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.gsim.base import ContextMaker
from openquake.hazardlib.calc.gmf import GmfComputer

from mapio.gmt import GMTGrid

import numpy as np
from scipy.stats import norm
from numpy import save
import os
import importlib

from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colormaps

from multiprocessing import Pool, Process, Lock, Manager
import config
import json
import sys

from openquake.hazardlib.geo import Point
from pyproj import Geod
geod = Geod(ellps='WGS84')
from scipy.stats import norm


# FUNCTION UTILITIES

    
def weighted_percentile(data, weights, perc):

    """
    perc : percentile in [0-1]!
    """
    ix = np.argsort(data)
    data = data[ix] # sort data
    weights = weights[ix] # sort weights
    cdf = (np.cumsum(weights) - 0.5 * weights) / np.sum(weights) # 'like' a CDF function
    return np.interp(perc, cdf, data)


def get_pois_coordinates_from_file(path, POIs_File):

    """
    Extract POIs coords from POIs file 
    """

    file = path + POIs_File

    POIs_coordinates = []
    with open(file, 'r') as f:
        for line in f:
            lat = float(line.strip().split()[0])
            lon = float(line.strip().split()[1])
            POIs_coordinates.append((lat, lon))

    return POIs_coordinates  


def get_pois(POIs_File):

    """
    Returns POIs coords and names + number of POIs in the file
    """

    path = os.path.join(os.getcwd(), "INPUT_FILES/")
    POIs_coordinates = get_pois_coordinates_from_file(path, POIs_File)
    LATs = [tup[0] for tup in POIs_coordinates]
    LONs = [tup[1] for tup in POIs_coordinates]

    POIs_lat, POIs_lon, POIs_NAMES = [], [], []

    for lat, lon in zip(LATs, LONs):

        POIs_lat.append(lat)
        POIs_lon.append(lon)
        POIs_NAMES.append(f"Site_LAT:{Point(float(lon),float(lat)).latitude}_LON:{Point(float(lon),float(lat)).longitude}")
        n_pois = len(POIs_NAMES)

    return POIs_lat, POIs_lon, POIs_NAMES, n_pois      


def get_params(ev_id):
        
    """
    Returns params needed for prob analysis
    """

    # Event dir
    event_dir = os.path.join(os.getcwd(), f"INPUT_FILES/events/{ev_id}")
    if not os.path.exists(event_dir):
        raise NotADirectoryError(f"{event_dir} is not a valid directory.")

    # Get the latitude and longitude of the event
    event_stat_file_path = os.path.join(event_dir, "event_stat.json")
    with open(event_stat_file_path, 'r', encoding='utf-8') as f:
        event_stat_file = json.load(f)
    coordinates = event_stat_file['features'][0]['geometry']['coordinates']
    Lon_Event = coordinates[0]  
    Lat_Event = coordinates[1]
    #Depth_Event = coordinates[2]
    prop = event_stat_file['features'][0]['properties']
    Time = prop['time']
    Mw = prop['mag']

    listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
    scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]

    # Get the number of scenarios in the file
    with open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r') as f:
        EnsembleSize = 0
        for _ in f:
            EnsembleSize += 1

    params = {}
    params['Lat_Event'] = Lat_Event
    params['Lon_Event'] = Lon_Event
    params['event_dir'] = event_dir
    params['Ensemble_Size'] = EnsembleSize
    params['Time'] = Time
    params['Mw'] = Mw

    return params


def load_obs(file):

    """
    Load ground-motion data from File
    """

    obs = []
    with open(file, 'r') as f:
        for line in f:
            tmp = float(line.strip().split()[0])
            obs.append((tmp))

    return np.array(obs)


def calculate_residuals(obs, mean, total_stdev):
    """
    Calculate residuals 
    """
    
    residuals = (np.log(obs) - mean) / total_stdev

    return residuals


def plot_normalized_residuals(pois_file, out_dir, mean_total_stdev_scens, data_file, Lon_Event, Lat_Event):

    POIs_lat, POIs_lon, POIs_NAMES, n_pois = get_pois(pois_file)
    POIs_lat = np.array(POIs_lat)
    POIs_lon = np.array(POIs_lon)

    out_path = os.path.join(os.getcwd(), f"{out_dir}/RANK_FIGURES")
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    obs = load_obs(data_file)
    mu = mean_total_stdev_scens[0] # shape: (n_pois,)
    total_stdev = mean_total_stdev_scens[1] # shape: (n_pois,)

    residuals = calculate_residuals(obs, mu, total_stdev) # shape: (n_pois,)

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
    cm = colormaps['seismic']

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
    
    max_abs = np.max(np.abs(residuals))
    sc = m.scatter(x, y, c=residuals, s=dim_point, cmap=cm, vmin=-max_abs, vmax=max_abs, edgecolor='k', linewidths=0.8)
    sc2 = m.scatter(x_ev, y_ev, s=100, c='k', marker='*', edgecolor='k', linewidths=0.3)

    plt.title('Normalized Residuals')
    cbar = plt.colorbar(sc)
    #cbar.set_label()
    #plt.show()

    figname = out_path + '/' + 'Normalized_Residuals.png'
    fig.savefig(figname, bbox_inches="tight", dpi=300, pad_inches=0.1)
    plt.close(fig)


def compute_GMPE_probabilities(mean_total_stdev_scens, data_file):
    """
    Compute probs for a given GMPE for gambling score 
    """

    obs = load_obs(data_file)
    obs = np.log(obs)

    mu = mean_total_stdev_scens[0] # shape: (n_pois,)
    total_stdev = mean_total_stdev_scens[1] # shape: (n_pois,)
    
    return norm.pdf(obs, loc=mu, scale=total_stdev)

##############################################################################
############################ PROBSHAKEMAP ####################################
##############################################################################

class Main:
    def __init__(self, gmpe, IMT, pois_file, NumGMPEsRealizations, num_processes, EnsembleSize):

        self.gmpe = gmpe
        self.imt = IMT
        self.pois_file = pois_file
        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.num_processes = num_processes
        self.vs30dir = os.path.join(os.getcwd(), f"INPUT_FILES/vs30")
        self.EnsembleSize = EnsembleSize
        
        self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
        self.POIs_lat = np.array(self.POIs_lat)
        self.POIs_lon = np.array(self.POIs_lon)

    def process_scenario(scen, Ensemble_Scenarios, msr, rupture_aratio, 
                        tectonicRegionType, context_maker, Site_Collection, 
                        correlation_model, crosscorr_model, gmpe, Num_Realiz):
        
        """
        For a given scenario, retrieves GMFs at all POIs
        """

        # Get scenario index 
        k = Ensemble_Scenarios.index(scen)

        Mag = float(scen[0])
        Hypocenter = Point(float(scen[1]), float((scen[2])), float(scen[3]))
        Rake = float(scen[6])
        Strike = float(scen[4])
        Dip = float(scen[5])

        planar_surface = surface.PlanarSurface.from_hypocenter(
        hypoc=Hypocenter,
        msr=msr(), 
        mag=Mag,
        aratio=rupture_aratio,
        strike=Strike,
        dip=Dip,
        rake=Rake,
        )

        source = BaseRupture(
        mag=Mag,
        rake=Rake,
        tectonic_region_type=tectonicRegionType,
        hypocenter=Hypocenter,
        surface=planar_surface
        )

        ctx = context_maker.get_ctxs([source], Site_Collection)
        gc = GmfComputer(
            source, 
            Site_Collection, 
            context_maker,
            correlation_model=correlation_model,
            cross_correl=crosscorr_model
            ) 

        # Get an array of shape (4, G, M, N) with mean and stddevs 
        # (G = number of GMPEs, M = number of IMTs, N = number of sites)
        mean_and_stdev_OQ = context_maker.get_mean_stds(ctx)  
        mean_total_stdev_OQ = [mean_and_stdev_OQ[0, 0, 0], mean_and_stdev_OQ[1, 0, 0]] # extract mean and total st dev only

        # Loop over GMPEs    
        gmf = []
        # 'gc.compute' --> Compute gmf and returns an array of shape (num_imts, num_sites, num_events) with the sampled gmf, and two arrays with shape 
        # (num_imts, num_events): sig for tau and eps for the random part
        # gmf[m] = numpy.exp(mean[g, m, :, None] + intra_res + inter_res)
        # Note that Num_Realiz[g] == num_events
        gf = gc.compute(gmpe, Num_Realiz, mean_and_stdev_OQ[:, 0, :, :])  

        # Append only first array output from gc.compute (shape: (num_imts, num_sites, num_events))  
        gmf.append(gf[0])

        return k, gmf, mean_total_stdev_OQ
    
    
    def run_prob_analysis(self):

        """
        Runs the prob analysis
        """

        print("********* STARTING PROB ANALYSIS *******")

        # Load configuration parameters
        config_dict = config.load_config('input_file.txt')

        tectonicRegionType = config_dict['tectonicRegionType']
        mag_scaling = config_dict['mag_scaling']
        rupture_aratio = config_dict['rupture_aratio']
        vs30file = config_dict['vs30file']
        CorrelationModel = config_dict['CorrelationModel']
        CrosscorrModel = config_dict['CrosscorrModel']
        vs30_clustering = config_dict['vs30_clustering']
        truncation_level = config_dict['truncation_level']
        seed = config_dict['seed']

        print("Number of source scenarios to process = ", self.EnsembleSize)
        print("Number of CPU processes: ", str(self.num_processes))

        # PRINT USER'S INPUT 
        print("TectonicRegionType: " + tectonicRegionType)
        print("Importing " + mag_scaling + " as magnitude scaling relationship")
        module = importlib.import_module('openquake.hazardlib.scalerel')
        msr = getattr(module, mag_scaling)
        print("Rupture aspect ratio: " + str(rupture_aratio))
        print("POIs file: " + self.pois_file)
        if vs30file == None:
            print("Vs30 file not provided")
        else:
            print("Vs30 file: " + vs30file)    
        print("Importing " + CorrelationModel + " as correlation model")
        module = importlib.import_module('openquake.hazardlib.correlation')
        correlation_model = getattr(module, CorrelationModel)
        print("Importing " + CrosscorrModel + " as crosscorrelation model")
        module = importlib.import_module('openquake.hazardlib.cross_correlation')
        crosscorr_model = getattr(module, CrosscorrModel)
        print("Vs30 clustering: " + str(vs30_clustering))
        print("Truncation level: " + str(truncation_level))
        print("Seed: " + str(seed))
        print("GMPE realizations per POI: " + str(self.NumGMPEsRealizations))
        print("Intensity measure type: " + str(self.imt))
       
        if vs30file is not None:
            vs30fullname = os.path.join(self.vs30dir, vs30file)
            if not os.path.isfile(vs30fullname):
                print(f"Warning: The file '{vs30fullname}' does not exist.")
                sys.exit()
            else:
                print("********* LOADING Vs30 *******")
                vs30grid = GMTGrid.load(vs30fullname)
                # Interpolate Vs30 values at POIs 
                vs30_POIs = vs30grid.getValue(self.POIs_lat, self.POIs_lon, method="nearest")

        print("********* DEFINING OpenQuake SITE COLLECTION *******")

        # Define a SiteCollection for all the POIs
        sites = []
        for i in range(len(self.POIs_NAMES)):
            site_location = Point(self.POIs_lon[i], self.POIs_lat[i])
            if vs30file == None:
                # If Vs30 file is not provided, use default value for Vs30
                site = Site(location=site_location, vs30=760., vs30measured=False, z1pt0=40., z2pt5=1.0)
            else:    
                site = Site(location=site_location, vs30=vs30_POIs[i], vs30measured=False, z1pt0=40., z2pt5=1.0)
            sites.append(site)
            
        Site_Collection = SiteCollection(sites)
        print(Site_Collection.complete)    

        print("********* BUILDING OPENQUAKE CONTEXTS *******")

        # Build OpenQuake contexts

        # # Define input parameters for ContextMaker

        imtls = {}
        imtls[self.imt] = []
        param = dict(imtls=imtls)   

        # Instantiate a ContextMaker object (Note: independent from source and sites!)
        context_maker = ContextMaker(tectonicRegionType, [self.gmpe], param)

        print("********* SAMPLING UNCERTAINTY *******")

        correlation_model = correlation_model(vs30_clustering=vs30_clustering)
        crosscorr_model = crosscorr_model(truncation_level=truncation_level)

        # Sample from current GMPE
        Num_Realiz = self.NumGMPEsRealizations  

        # Sample from the total variability of ground motion taking into account both inter- and intra-event variability (for one source scenario only)
        # gmf = exp(mu + crosscorel(tau) + spatialcorrel(phi)) --> See: https://docs.openquake.org/oq-engine/advanced/latest/event_based.html#correlation-of-ground-motion-fields

        GMPEsRealizations = [None] * self.EnsembleSize

        listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
        scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]
        f = open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r')
        print("List of scenarios: ", scenarios_file[0])

        Ensemble_Scenarios = []
        for k, line in enumerate(f):
            scen = line.strip().split(' ')
            Ensemble_Scenarios.append(scen)  

        # Set the random seed for reproducibility in OpenQuake GmfComputer
        np.random.seed(seed)

        ################################
        # SETTING MULTIPROCESSING PARAMS
        ################################

        chunk_size_default = int(self.EnsembleSize/self.num_processes) # size of each chunk of scenarios
        #print("Chunk size = ", chunk_size_default)
        last_chunk_size = chunk_size_default + self.EnsembleSize - self.num_processes * chunk_size_default # size of the last chunk
        #print("Last chunk size = ", last_chunk_size)

        # Create pool of worker processes
        with Pool(processes=self.num_processes) as pool:
            results = []

            # iterate over processes
            for i in range(self.num_processes):
                if i == self.num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == self.num_processes - 1:
                    start_idx = self.EnsembleSize - chunk_size
                    end_idx = self.EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk:
                    result = Main.process_scenario(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                            msr=msr, rupture_aratio=rupture_aratio, tectonicRegionType=tectonicRegionType,
                                            context_maker=context_maker, Site_Collection=Site_Collection, 
                                            correlation_model=correlation_model, crosscorr_model=crosscorr_model, gmpe=self.gmpe, 
                                            Num_Realiz=Num_Realiz)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        # Combine results
        mean_total_stdev_scens_OQ = [None] * self.EnsembleSize
        for result in results:
            k = result[0]
            gmf = result[1]
            mean_total_stdev_OQ = result[2]
            GMPEsRealizations[k] = gmf
            mean_total_stdev_scens_OQ[k] = mean_total_stdev_OQ

            # PRINTING INFO
            if k == 0:
                # Print this for one source scenario only
                print("IMT: ", self.imt, "-- GMPE", self.gmpe, "is sampled", Num_Realiz, "times")

        # Results AGGREGATION
        # Aggregate the generated gmf at each site for Probabilistic Shakemap

        # Structure of GMPEsRealizations
        # 1st Index: Scenario index
        # 2nd Index: GMPE index
        # 3rd Index: IMT index
        # 4th Index: Site index 

        # CHECK
        # print("SHAPE = ", len(GMPEsRealizations))
        # print("SHAPE = ", len(GMPEsRealizations[0]))
        # print("SHAPE = ", len(GMPEsRealizations[0][0]))
        # print("SHAPE = ", len(GMPEsRealizations[0][0][0]))

        # For each site, there are as many values as the number of realizations for the current GMPE 

        # PREPARE KEYS FOR SCENARIOS AND SITES
        keys_sites = [] 
        for s in range(len(sites)):
            keys_sites.append(f"Site_LAT:{Point(float(self.POIs_lon[s]), float(self.POIs_lat[s])).latitude}_LON:{Point(float(self.POIs_lon[s]), float(self.POIs_lat[s])).longitude}")

        keys_scen = []
        for k in range(self.EnsembleSize):
            keys_scen.append(f"Scenario_{k+1}") 

        # # Structure of SiteGmf 
        # 1st Index: Scenario index 
        # 2nd Index: Site index 

        # Re-structure output such that, for a given scenario, we can access GMPE realizations at each site 
        SiteGmf = np.empty((self.EnsembleSize, len(sites), self.NumGMPEsRealizations), dtype=np.float32)
        for scen in Ensemble_Scenarios:
            k = Ensemble_Scenarios.index(scen)

            SiteGmf_scen = np.empty((len(sites), self.NumGMPEsRealizations), dtype=object)
            # Loop over sites
            for s in range(len(sites)):
                SiteGmf_scen[s] = GMPEsRealizations[k][0][0][s]

            SiteGmf[k] = SiteGmf_scen

        print("********* PROB ANALYSIS DONE! *******")

        # CHECK
        # print(SiteGmf_scen.shape)
        # print(SiteGmf.shape)

        prob_output = {
            "SiteGmf": SiteGmf,
            "keys_scen": keys_scen,
            "keys_sites": keys_sites,
            "mean_total_stdev_scens_OQ": mean_total_stdev_scens_OQ
        }

        return prob_output


class GetStatistics:
    def __init__(self, out_dir, SiteGmf, EnsembleSize, Lon_Event, Lat_Event, NumGMPEsRealizations,
                  event_dir, IMT, imt_min, imt_max, pois_file, num_processes, vector_npy, fileScenariosWeights):

        self.out_dir = out_dir
        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.SiteGmf = SiteGmf
        self.EnsembleSize =  EnsembleSize
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.event_dir = event_dir
        self.imt = IMT
        self.imt_min = imt_min
        self.imt_max = imt_max
        self.pois_file = pois_file
        self.num_processes = num_processes
        self.vector_npy = bool(vector_npy)
        self.fileScenariosWeights = fileScenariosWeights

        print(f"Save vector.npy set to {self.vector_npy}")

        self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
        self.POIs_lat = np.array(self.POIs_lat)
        self.POIs_lon = np.array(self.POIs_lon)
        print("Found ", self.n_pois, "POIs")
        self.i_pois = list(range(self.n_pois))


    def process_scen_gmf(scen, Ensemble_Scenarios, site_gmf, i_pois, n_pois, num_realizations, weights):
        
        """
        Prepare for POI-level statistics
        """

        # Get scenario index 
        i_scen = Ensemble_Scenarios.index(scen)
     
        vector_scen = np.zeros([n_pois, num_realizations])
        weight_scen = np.zeros([n_pois, num_realizations]) 

        for jp in range(n_pois):

            vector_scen[jp] = site_gmf[i_scen][i_pois[jp]]
            weight_scen = weights[i_scen] / num_realizations

        return i_scen, vector_scen, weight_scen

    def calc_statistics(self):

        """
        Calculates the statistics of the ground motion distributions at all POIs
        """

        # Default weights
        weights = np.ones(self.EnsembleSize) / self.EnsembleSize

        # Load weights from file if provided
        if self.fileScenariosWeights:
            print('Weights file:', self.fileScenariosWeights)
            with open(self.fileScenariosWeights) as f:
                weights = np.array([float(line.split()[0]) for line in f.readlines()])

        weights = weights / np.sum(weights)  # Normalize weights

        listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
        scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]
        f = open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r')
        Ensemble_Scenarios = []
        for _, line in enumerate(f):
            scen = line.strip().split(' ')
            Ensemble_Scenarios.append(scen) 

        # Chunking
        chunk_size_default = int(self.EnsembleSize/self.num_processes) # size of each chunk of scenarios
        last_chunk_size = chunk_size_default + self.EnsembleSize - self.num_processes * chunk_size_default # size of the last chunk

        # Create pool of worker processes
        with Pool(processes=self.num_processes) as pool:
            results = []
        
            # iterate over processes
            for i in range(self.num_processes):
                if i == self.num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == self.num_processes - 1:
                    start_idx = self.EnsembleSize - chunk_size
                    end_idx = self.EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk: 
                    result = GetStatistics.process_scen_gmf(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                                             site_gmf=self.SiteGmf, i_pois=self.i_pois, n_pois=self.n_pois, 
                                                             num_realizations=self.NumGMPEsRealizations, weights=weights)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        # 'vector' gathers, for each POI, the ground motion distribution collecting
        # GMFs from current GMPE across all scenarios 
        vector = np.zeros((self.n_pois, self.NumGMPEsRealizations * self.EnsembleSize))
        weight = np.zeros((self.n_pois, self.NumGMPEsRealizations * self.EnsembleSize)) 
        mean_and_stdev = {name: [0] * self.n_pois for name in ['Mean', 'ST_DEV']}

        # Collect results
        for result in results:  
            i_scen = result[0]
            vector_scen = result[1]
            weight_scen = result[2]
            start = i_scen * self.NumGMPEsRealizations
            end = (i_scen + 1) * self.NumGMPEsRealizations
            vector[:, start:end] = vector_scen
            weight[:, start:end] = weight_scen

        # Calculate statistics
        for jp in range(self.n_pois):
            # Convert to log scale and calculate (weighted) stats
            mean_and_stdev['Mean'][jp] = np.sum(np.log(vector[jp, :]) * weight[jp, :]) / np.sum(weight[jp, :])
            p84 = weighted_percentile(np.log(vector[jp, :]), weight[jp, :], 0.84)
            p16 = weighted_percentile(np.log(vector[jp, :]), weight[jp, :], 0.16)
            mean_and_stdev['ST_DEV'][jp] = (p84 - p16) / 2

        stats = {'vector': vector, 'weight': weight, 'mean_and_stdev': mean_and_stdev}

        return stats
    
    def save_statistics(self):

        """
        Saves 'vector_stat.npy' with mean and st dev
        [OPTIONAL] Saves 'vector.npy' with all ground motion distributions at all POIs (can be huge!)
        """

        stats = GetStatistics.calc_statistics(self)
        vector_stat = stats['mean_and_stdev'] 
        vector = stats['vector'] 

        path = os.path.join(os.getcwd(), f"{self.out_dir}/npyFiles")
        if not os.path.exists(path):
            os.makedirs(path)
        print("********* SAVING STATISTICS *******")
        gmpe_ref = self.out_dir.split('_')[-2]
        save(path + '/' + f"{self.imt}_{gmpe_ref}" + '_vector_stat.npy', vector_stat)
        if self.vector_npy == True:
            save(path + '/' + f"{self.imt}_{gmpe_ref}" + '_vector.npy', vector)

    def plot_statistics(self):

        """
        Visualizes statistics from 'vector_stat.npy' on a map and saves it
        """

        print("********* GENERATING STATISTICS PLOTS *******")

        path = os.path.join(os.getcwd(), f"{self.out_dir}/STATISTICS")
        if not os.path.exists(path):
            os.makedirs(path)

        stats = GetStatistics.calc_statistics(self)
        vector_stat = stats['mean_and_stdev'] 

        dim_point = 40

        min_lon, max_lon = np.min(self.POIs_lon), np.max(self.POIs_lon)
        min_lat, max_lat = np.min(self.POIs_lat), np.max(self.POIs_lat)
        lon_span = max_lon - min_lon
        lat_span = max_lat - min_lat
        buffer_lon = max(0.05 * lon_span, 0.1)  
        buffer_lat = max(0.05 * lat_span, 0.1)

        xlim_min = max(-180, min_lon - buffer_lon)
        xlim_max = min(180, max_lon + buffer_lon)
        ylim_min = max(-90, min_lat - buffer_lat)
        ylim_max = min(90, max_lat + buffer_lat)

        for name in ['Mean', 'ST_DEV']:

            fig = plt.figure(figsize=(9, 6))

            latitudes = np.arange(-90, 91, 2)
            longitudes = np.arange(-180, 181, 2)

            cm = colormaps['turbo']

            m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
                    llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i')
            try:
                m.drawcoastlines()
            except Exception:
                pass 

            m.drawparallels(latitudes, labels=[1,0,0,0], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])
            m.drawmeridians(longitudes, labels=[0,0,0,1], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])

            m.drawmapboundary(linewidth=2, color='black', fill_color='white')

            x, y = m(self.POIs_lon, self.POIs_lat)
            x_ev, y_ev = m(self.Lon_Event, self.Lat_Event)
            
            tmp = vector_stat[name]
            # Note: tmp is already in log scale, imt_min and imt_max are provided in linear space by the user 
            if name == 'Mean':
                sc = m.scatter(x, y, c=tmp, vmin=np.log(self.imt_min), vmax=np.log(self.imt_max),  edgecolor='k', linewidths=0.8, s=dim_point, cmap = cm)
            else:
                sc = m.scatter(x, y, c=tmp, vmin=min(tmp), vmax=max(tmp), s=dim_point,  edgecolor='k', linewidths=0.8, cmap = cm)
            
            sc2 = m.scatter(x_ev, y_ev, s=100, c='k', marker='*')

            plt.title(f"{name} (ln scale)")
            cbar = plt.colorbar(sc)
            cbar.set_label(f"ln({self.imt})")
            #plt.show()

            figname = path + '/' + 'Stat_' + name.strip() + '.png'
            fig.savefig(figname, bbox_inches="tight", dpi=300, pad_inches=0.1)
            plt.close(fig)

        print(f"***** Figures saved in {path} *****")


class LLH_Score:
    def __init__(self, data_file, mean_total_stdev_scens, pois_file=None, out_dir=None, Lon_Event=None, Lat_Event=None):

        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.pois_file = pois_file
        self.out_dir = out_dir
        self.data_file = data_file
        self.mean_total_stdev_scens = mean_total_stdev_scens

        if pois_file is not None:
            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

    def calculate_llh(self):

        '''Inspired from 
        https://github.com/GEMScienceTools/gmpe-smtk/blob/master/smtk/residuals/gmpe_residuals.py
        '''

        obs = load_obs(self.data_file)

        # For each POI in the grid, calculate LLH score. 
        # For a given POI, residuals are calculated using mean and total st dev from ProbShakemap predictive distribution at that POI.

        mu = self.mean_total_stdev_scens[0] # shape: (n_pois,)
        total_stdev = self.mean_total_stdev_scens[1] # shape: (n_pois,)
        # print("CHECK!")
        # assert len(obs) == len(mu) == len(total_stdev), "Mismatch in POI lengths"

        residuals = calculate_residuals(obs, mu, total_stdev) # shape: (n_pois,)

        # This LLH score already incorporates source uncertainty at a given POI 
        # (mean and st dev are from ProbShakemap predictive distribution, which aggregates gmfs from all scenarios)
        asll = np.log2(norm.pdf(residuals, 0., 1.0))  # array of log-likelihoods per residual, shape: (n_pois,)
        
        # Average LLH score (averaged over POIs)
        llh = -(1.0 / float(len(asll))) * np.sum(asll)

        return asll, llh

    def LLH_plot(self, asll):

        """
        Plot point-LLH score (i.e. LLH score at each POI)
        """

        if self.POIs_lat is None or self.out_dir is None or self.Lon_Event is None or self.Lat_Event is None:
            raise ValueError("LLH_plot requires: POIs file, output directory, Lat and Lon Event.")

        out_path = os.path.join(os.getcwd(), f"{self.out_dir}/RANK_FIGURES")
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        dim_point = 40

        min_lon, max_lon = np.min(self.POIs_lon), np.max(self.POIs_lon)
        min_lat, max_lat = np.min(self.POIs_lat), np.max(self.POIs_lat)
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

        x, y = m(self.POIs_lon, self.POIs_lat)
        x_ev, y_ev = m(self.Lon_Event, self.Lat_Event)
        
        vmin = np.min(asll)
        vmax = np.max(asll)
        sc = m.scatter(x, y, c=asll, s=dim_point, cmap = cm, vmin=vmin, vmax=vmax, edgecolor='k', linewidths=0.8)
        sc2 = m.scatter(x_ev, y_ev, s=100, c='k', marker='*', edgecolor='k', linewidths=0.3)

        plt.title('POI-LLH')
        cbar = plt.colorbar(sc)
        #cbar.set_label()
        #plt.show()

        figname = out_path + '/' + 'POI_LLH.png'
        fig.savefig(figname, bbox_inches="tight", dpi=300, pad_inches=0.1)
        plt.close(fig)

        






       
