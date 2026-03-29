import matplotlib.pyplot as plt
import numpy as np
import ABM_multrun_Active as multrun
from ABM_variance_Active import calculate_scale_and_RMSE

def plot_model_vs_data(medians, low_iqr, high_iqr, params):
    """
    Plots the scaled model results (Median + IQR) vs experimental medians.
    """
    # 1. Extract measure points to find the optimal scaling factor
    # This ensures the plot matches the calibration logic [cite: 111, 265]
    measure_points = multrun.measure_points(medians)
    s_factor, _ = calculate_scale_and_RMSE(measure_points)
    
    # 2. Configuration for CD4 Subsets
    subsets = ["%TSCM", "%TCM", "%TEMRA"]
    titles = ["TSCM (Persistence)", "TCM (Memory)", "TEMRA (Effector)"]
    
    data_medians = {
        "%TSCM": [0.000, 0.218, 0.172, 0.246],
        "%TCM": [0.000, 0.381, 0.339, 0.178],
        "%TEMRA": [0.000, 0.072, 0.092, 0.034]
    }
    time_points = [1, 22, 43, 365]
    
    # 3. Setup the figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    days = np.arange(366) 

    for i, s in enumerate(subsets):
        ax = axes[i]
        
        # Apply scaling: (Raw Count * s_factor) = %AIM+ [cite: 110, 122]
        scaled_median = medians[s] * s_factor
        scaled_low = low_iqr[s] * s_factor
        scaled_high = high_iqr[s] * s_factor
        
        # Plot Model IQR Ribbon
        ax.fill_between(days, scaled_low, scaled_high, color='royalblue', alpha=0.2, label="Model IQR")
        
        # Plot Model Median Line
        ax.plot(days, scaled_median, color='royalblue', linewidth=2, label="Model Median")
        
        # Plot Data Medians (Red X)
        ax.scatter(time_points, data_medians[s], color='red', marker='x', s=80, 
                   zorder=5, label="Exp. Median")
        
        # Formatting
        ax.set_title(f"{titles[i]}\n(Scale S={s_factor:.1f})")
        ax.set_xlabel("Days")
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        if i == 0:
            ax.set_ylabel("Response Magnitude (%AIM+)")
            ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Example usage with your fixed/best parameters [cite: 99]
    test_params = {
        "mu_N": 0.0003, "mu_TSCM": 0.0002, "mu_TCM": 0.004, "mu_TEM": 0.01,
        "mu_TEMRA": 0.02, "mu_MPEC": 0.02, "mu_SLEC": 0.05, "f_TSCM": 0.03,
        "f_TCM": 0.05, "f_TEM": 0.06, "f_TEMRA": 0.02, "t_peak": 18,
        "sigma": 7, "q": 0.10867425382641213, "alpha_peak": 0.1060179817341049, "b_MPEC": 1, 
        "b_SLEC": 1, "K_mem": 172, "S_CD4": 1655
    }
    
    print("Gathering data for 50 seeds...")
    m, low, high = multrun.get_plot_data(test_params, num_seeds=50)
    
    plot_model_vs_data(m, low, high, test_params)