import numpy as np
from multiprocessing import Pool
import ABM_model_Active as modelV4

def run_single_seed(params):
    # initialize model
    model = modelV4.Immunology_Model(**params)
    # run model
    for _ in range(365):
        model.step()
    # return data
    return model.datacollector.model_vars

def multrun(params, num_seeds=50):
    # run the model x times using parallel processing.
    worker_args = [params for _ in range(num_seeds)]
    with Pool() as pool:
        results_list = pool.map(run_single_seed, worker_args)
    return results_list

def calculate_multrun_medians(results_list):
    # keys = different columns ["Naive", "TSCM", ...]
    keys = list(results_list[0].keys())
    # num_days = 365
    num_days = len(results_list[0][keys[0]])
    # num_seeds = 50
    num_seeds = len(results_list)
    # key_to_idx = {"Naive": 0, "TSCM": 1, ...} 
    key_to_idx = {key: i for i, key in enumerate(keys)}
    # create cubic data (first all zero) with arguments defining it's size
    data_cube = np.zeros((len(keys), num_seeds, num_days))
    # enumerate the runs
    for seed_idx, seed_data in enumerate(results_list):
        # Insert the data key by key
        for key, daily_values in seed_data.items():
            data_cube[key_to_idx[key], seed_idx, :] = daily_values
    # Now we want to calculate the medians across the runs axis giving us the complete median dataframe.
    medians = np.median(data_cube, axis=1)
    median_dict = {key: medians[key_to_idx[key]] for key in keys}
    return median_dict

import numpy as np

def measure_points(median_dict):
    # Important time points
    time_points = [0, 21, 42, 364]
    # subsets to measure.
    subsets = ["%TSCM", "%TCM", "%TEMRA"]
    
    calibration_data = {}
    
    for s in subsets:
        # Extract the 4 specific values from the 365-day history
        calibration_data[s] = [median_dict[s][i] for i in time_points]
        
    return calibration_data

def calculate_multrun_IQR(results_list):
    """
    Calculates the 25th and 75th percentiles (IQR) across seeds 
    for the uncertainty ribbons.
    """
    keys = list(results_list[0].keys())
    num_days = len(results_list[0][keys[0]])
    num_seeds = len(results_list)
    key_to_idx = {key: i for i, key in enumerate(keys)}
    
    # Create the data cube
    data_cube = np.zeros((len(keys), num_seeds, num_days))
    for seed_idx, seed_data in enumerate(results_list):
        for key, daily_values in seed_data.items():
            data_cube[key_to_idx[key], seed_idx, :] = daily_values
            
    # Calculate percentiles across the runs axis (axis=1) 
    low_bound = np.percentile(data_cube, 25, axis=1)
    high_bound = np.percentile(data_cube, 75, axis=1)
    
    # Format into dictionaries for easy plotting
    low_dict = {key: low_bound[key_to_idx[key]] for key in keys}
    high_dict = {key: high_bound[key_to_idx[key]] for key in keys}
    
    return low_dict, high_dict

def get_plot_data(params, num_seeds=50):
    """
    High-level helper to run the simulation and return all plotting 
    components: Medians and IQR bounds.
    """
    # 1. Run the parallel simulations [cite: 261]
    raw_results = multrun(params, num_seeds=num_seeds)
    
    # 2. Extract Medians
    medians = calculate_multrun_medians(raw_results)
    
    # 3. Extract IQR Bounds 
    low_iqr, high_iqr = calculate_multrun_IQR(raw_results)
    
    return medians, low_iqr, high_iqr

if __name__ == "__main__":
    
# --- Example Usage ---
    # FREE parameters
    alpha_peak = 0.28915212707842
    b_MPEC = 1
    K_mem = 243
    S_CD4 = 5738
    # FIXED parameters
    mu_N = 0.0003
    mu_TSCM = 0.0002
    mu_TCM = 0.004
    mu_TEM = 0.01
    mu_TEMRA = 0.02
    mu_MPEC = 0.02
    mu_SLEC = 0.05
    f_TSCM = 0.03
    f_TCM = 0.05
    f_TEM = 0.06
    f_TEMRA = 0.02
    t_peak = 18
    sigma = 7
    q = 0.35
    b_SLEC = b_MPEC

    params = {
        "mu_N": mu_N,
        "mu_TSCM": mu_TSCM,
        "mu_TCM": mu_TCM,
        "mu_TEM": mu_TEM,
        "mu_TEMRA": mu_TEMRA,
        "mu_MPEC": mu_MPEC,
        "mu_SLEC": mu_SLEC,
        "f_TSCM": f_TSCM,
        "f_TCM": f_TCM,
        "f_TEM": f_TEM,
        "f_TEMRA": f_TEMRA,
        "b_MPEC": b_MPEC,
        "b_SLEC": b_SLEC,
        "q": q,
        "K_mem": K_mem,
        "S_CD4": S_CD4,
        "alpha_peak": alpha_peak,
        "t_peak": t_peak,
        "sigma": sigma
    }
    multrun(num_seeds=100, params=params)