import numpy as np
import copy
import optuna

import matplotlib.pyplot as plt
from modelC1_variance import calculate_scale_and_RMSE
import modelC1_multrun as multrun 

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

def run_1d_sensitivity(base_params, free_param_names, steps=11):
    """
    Varies each explicitly 'free' parameter by +/- 50% and calculates RMSE.
    """
    sensitivity_results = {}
    
    for param_name in free_param_names:
        print(f"Sweeping parameter: {param_name}")
        
        best_val = base_params[param_name]
        lower_bound = best_val * 0.5
        upper_bound = best_val * 1.5
        
        test_values = np.linspace(lower_bound, upper_bound, steps)
        
        # Force integer values for specific parameters
        if isinstance(best_val, int) or param_name in ["b_MPEC", "S_CD4", "K_mem"]:
            test_values = np.unique(np.round(test_values).astype(int))
            
        param_rmses = []
        for val in test_values:
            test_params = copy.deepcopy(base_params)
            test_params[param_name] = val
            
            # Keep b_SLEC tied to b_MPEC if b_MPEC is the one being swept
            if param_name == "b_MPEC":
                test_params["b_SLEC"] = val
                
            rmse = evaluate_rmse(test_params)
            param_rmses.append(rmse)
            
        sensitivity_results[param_name] = {
            "values": test_values,
            "rmses": param_rmses,
            "best_val": best_val
        }
        
    return sensitivity_results

def plot_sensitivity(results, filename="sensitivity_plots.png"):
    """
    Plots RMSE vs Parameter value and saves it to a file, bypassing X11.
    """
    num_params = len(results)
    cols = 2
    rows = int(np.ceil(num_params / cols))
    
    # Handle case where there's only 1 or 2 parameters gracefully
    fig, axes = plt.subplots(rows, cols, figsize=(12, 4 * rows))
    
    # Ensure axes is always a flat array for easy iteration
    if num_params == 1:
        axes = [axes]
    else:
        axes = axes.flatten() 
    
    for idx, (param_name, data) in enumerate(results.items()):
        ax = axes[idx]
        
        ax.plot(data["values"], data["rmses"], marker='o', linestyle='-', color='b')
        ax.axvline(x=data["best_val"], color='r', linestyle='--', label=f'Best Value ({data["best_val"]:.3g})')
        
        ax.set_title(f"Sensitivity: {param_name}")
        ax.set_xlabel(f"{param_name} Value")
        ax.set_ylabel("RMSE")
        ax.grid(True, linestyle=':', alpha=0.7)
        ax.legend()
        
    # Hide any unused subplots
    for idx in range(num_params, len(axes)):
        if isinstance(axes, np.ndarray) or isinstance(axes, list):
            fig.delaxes(axes[idx])
        
    plt.tight_layout()
    
    # Save the figure instead of showing it
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    print(f"Plot saved successfully to {filename}")
    
    plt.close(fig)


if __name__ == "__main__":
    study = optuna.load_study(study_name="yellowfever_250seeds_q25_t14s5", storage="sqlite:///yellowfever_250seeds_q25_t14s5.db")

    params = {
        "alpha_peak": 0.01, "b_MPEC": 1, "K_mem": 250,
        "S_CD4": 500, "q": 0.25, "mu_N": 0.0003,
        "mu_TSCM": 0.0002, "mu_TCM": 0.004, "mu_TEM": 0.01,
        "mu_TEMRA": 0.02, "mu_MPEC": 0.02, "mu_SLEC": 0.05,
        "f_TSCM": 0.03, "f_TCM": 0.05, "f_TEM": 0.06,
        "f_TEMRA": 0.02, "t_peak": 14, "sigma": 5, "cc": 0.5
    }
    
    params.update(study.best_params)
    params["b_SLEC"] = params["b_MPEC"]

   

    free_params_to_sweep = list(study.best_params.keys())

    sensitivity_data = run_1d_sensitivity(
        base_params=params, 
        free_param_names=free_params_to_sweep, 
        steps=11
    )
    
    plot_sensitivity(sensitivity_data, filename="yellowfever_fit_sensitivity_q25_t14s5.png")