import ABM_model_Active as modelV4
import ABM_variance_Active as varianceV4
import matplotlib.pyplot as plt
import numpy as np
import ABM_scalingfactor_Active as scalingfactor


def mult_run(num_runs=20, params=None):
    subsets = ["%N", "%S", "%C", "%R", "%Total_memory"]
    titles = ["Naive (N)", "TSCM (S)", "TCM (C)", "TEMRA (R)", "Total Memory"]
    time_points = [1, 22, 43, 365]
    medians = [
        [0.349, 1.005, 0.039, 0.0], [0.0, 0.218, 0.171, 0.246],
        [0.0, 0.380, 0.338, 0.178], [0.0, 0.072, 0.091, 0.033],
        [0.0, 0.670, 0.602, 0.458]
    ]

    # Initialize storage for all runs
    all_runs_data = {subset: [] for subset in subsets}
    
    # 1. Execute Runs
    for i in range(num_runs):
        print(f"Simulation {i+1}/{num_runs}...")
        
        # Pass the dictionary to your model
        # Assuming your class is defined as: def __init__(self, parameters):
        
        model = modelV4.Immunology_Model(**params)
        for _ in range(365):
            model.step()
        
        # Apply your scaling logic
        scalingfactor.Calculate_global_factor(model)
        
        for subset in subsets:
            scaled_history = [val * model.global_factor for val in model.datacollector.model_vars[subset]]
            all_runs_data[subset].append(scaled_history)

    # 2. Plotting
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for i, subset in enumerate(subsets):
        ax = axes[i]
        data_matrix = np.array(all_runs_data[subset])
        
        # Stats for the "Mean Path" and "Stability Envelope"
        mean_path = np.mean(data_matrix, axis=0)
        std_path = np.std(data_matrix, axis=0)
        
        # Plot individual fainted runs
        for run_idx in range(num_runs):
            ax.plot(range(1, 367), data_matrix[run_idx], color='royalblue', alpha=0.1, linewidth=0.8)
        
        # Plot Mean Path
        ax.plot(range(1, 367), mean_path, color='navy', linewidth=2.5, label="Mean Model")
        
        # Plot Confidence Interval (Shaded Area)
        ax.fill_between(range(1, 367), mean_path - std_path, mean_path + std_path, 
                        color='royalblue', alpha=0.2, label="±1 Std Dev")
        
        # Plot Experimental Medians
        ax.scatter(time_points, medians[i], color='red', marker='x', s=100, zorder=10, label="Exp. Median")
        
        ax.set_title(titles[i])
        ax.grid(True, linestyle='--', alpha=0.4)
        if i == 0: ax.legend()

    axes[5].axis('off')
    plt.tight_layout()
    plt.show()





if __name__ == "__main__":
    
# --- Example Usage ---
    base_params = {'T_antigen': 22, 'initial_population': 1000, 'p_act': 0.7130575836504838, 'frac_memory': 0.81985770893272, 'S_C_Ratio': 0.13072410353235184, 'p_SC': 3.112310238392141e-05, 'p_CR': 0.0014285840051584267, 'r_N': 0.049666033337027966, 'r_S': 0.06515909784610464, 'r_C': 0.0, 'r_R': 0.0, 'd_N': 0.0035584268482424285, 'd_S': 4.731783078114052e-06, 'd_C': 0.0007002347364540852, 'd_R': 0.017535584696647337}
    mult_run(num_runs=100, params=base_params)