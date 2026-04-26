import numpy as np
import multiprocessing
import modelB1_model as modelB1
import random


def run_single_seed(args):
    params, seed = args
    random.seed(seed)
    model = modelB1.Immunology_Model(**params, seed=seed)
    # run model
    for _ in range(365):
        model.step()
    return model.datacollector.model_vars

def multrun(params, num_seeds=50):
    # run the model x times using parallel processing.
    worker_args = [(params, seed) for seed in range(num_seeds)]
    with multiprocessing.Pool() as pool:
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
    # 1. Run the simulation
    raw_results = multrun(params, num_seeds=num_seeds)
    
    # 2. Extract Medians
    medians = calculate_multrun_medians(raw_results)
    
    # 3. Extract IQR Bounds 
    low_iqr, high_iqr = calculate_multrun_IQR(raw_results)
    
    return medians, low_iqr, high_iqr