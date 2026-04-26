import numpy as np
import copy
import optuna

import matplotlib.pyplot as plt
from ABM_variance_Active import calculate_scale_and_RMSE
import ABM_multrun_Active as multrun 

def evaluate_rmse(params):
    """
    Runs the model with a given set of parameters and returns the RMSE.
    """
    data = multrun.multrun(params,num_seeds=250) 
    medians = multrun.calculate_multrun_medians(data)
    measured = multrun.measure_points(medians)
    
    # Assuming this function is available in your namespace
    _, RMSE = calculate_scale_and_RMSE(measured) 
    return RMSE
def temp(study):
    params = {
        "alpha_peak": 0.01, "b_MPEC": 1, "K_mem": 999_999,
        "S_CD4": 500, "q": 0.25, "mu_N": 0.0003,
        "mu_TSCM": 0.0002, "mu_TCM": 0.004, "mu_TEM": 0.01,
        "mu_TEMRA": 0.02, "mu_MPEC": 0.02, "mu_SLEC": 0.05,
        "f_TSCM": 0.03, "f_TCM": 0.05, "f_TEM": 0.06,
        "f_TEMRA": 0.02, "t_peak": 18, "sigma": 7
    }
    
    params.update(study.best_params)
    #params["b_MPEC"] = 1
    params["b_SLEC"] = params["b_MPEC"]
    #params["q"] = 0.25
    #params["K_mem"] = 250
    #params["S_CD4"] = 1000
    #params["alpha_peak"] = 0.25
#“Params: alpha=0.25, bM=1, Km=250, S=1000, q=0.25”
    print(evaluate_rmse(params))
    return




# 1. First Trial
study = optuna.load_study(
    study_name="yellowfever_250seeds_q25_ablation", 
    storage="sqlite:///yellowfever_250seeds_q25_ablation.db"
)
temp(study)
input()
# 2. Second Trial (Ablation Kmem)
study = optuna.load_study(
    study_name="ablate_kmem", # <-- Put the name you found in Step 1 here
    storage="sqlite:///50-seedtrial/ablation_kmem.db"
)
temp(study)

# 3. Third Trial (250 Seeds Normal)
study = optuna.load_study(
    study_name="typhoid_fit_250seeds", 
    storage="sqlite:///250seeds_normal.db"
)
temp(study)

# 4. Fourth Trial (250 Seeds Ablation)
study = optuna.load_study(
    study_name="typhoid_fit_250seeds_ablate", 
    storage="sqlite:///250seeds_ablation.db"
)
temp(study)