import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def replot_simulation(csv_path):
    if not os.path.exists(csv_path):
        print(f"Error: File '{csv_path}' not found.")
        return

    df = pd.read_csv(csv_path, index_col=0) 

    medians = {
        'time_points': [1, 22, 43, 365],
        'S': [0.0, 0.218, 0.171, 0.246],
        'C': [0.0, 0.380, 0.338, 0.178],
        'R': [0.0, 0.072, 0.091, 0.033],
        'Total-memory': [0.0, 0.670, 0.602, 0.458]
    }

    # --- MSE CALCULATION ---
    total_mse = 0
    count = 0
    
    for cell in ['S', 'C', 'R', 'Total-memory']:
        if cell in df.columns:
            for i, t in enumerate(medians['time_points']):
                # Find the simulation value at that day (t)
                if t in df.index:
                    sim_val = df.loc[t, cell]
                    exp_val = medians[cell][i]
                    total_mse += (sim_val - exp_val)**2
                    count += 1
    
    avg_mse = total_mse / count if count > 0 else 0
    # -----------------------

    max_day = df["runtime_days"].iloc[-1] if "runtime_days" in df.columns else df.index[-1]
    p_act = df["param_p_act"].iloc[0] if "param_p_act" in df.columns else "N/A"
    
    cell_types = ['S', 'C', 'R', 'Total-memory']
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, cell_type in enumerate(cell_types):
        ax = axes[i]
        if cell_type in df.columns:
            ax.plot(df.index, df[cell_type], label='Simulation', color='tab:blue', linewidth=2.5)
        
        if cell_type in medians:
            ax.scatter(medians['time_points'], medians[cell_type], color='red', s=60, zorder=5)

        ax.set_title(cell_type, fontweight='bold')
        ax.set_xlim(0, max_day)
        ax.set_ylim(0, 1.05)
        ax.grid(True, linestyle=':', alpha=0.7)

    # Adding the Error Score to the title
    plt.suptitle(f"Replay: {os.path.basename(csv_path)}\nParams: p_act={p_act} | Duration={max_day} Days | MSE Score: {avg_mse:.6f}", 
                 fontsize=14, y=0.98)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    plt.show()

# Usage:
# replot_simulation("immuno_results_test.csv")
replot_simulation("immuno_results_224554.csv")