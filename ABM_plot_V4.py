import matplotlib.pyplot as plt

def plot_model_vs_data(run):
    time_points = [1, 22, 43, 365]
    indices = [0, 21, 42, 364]
    
    # Data labels and medians
    subsets = ["%N", "%S", "%C", "%R", "%Total_memory"]
    titles = ["Naive (N)", "TSCM (S)", "TCM (C)", "TEMRA (R)", "Total Memory"]
    medians = [
        [0.349, 1.005, 0.039, 0.0],    # N
        [0.0, 0.218, 0.171, 0.246],    # S
        [0.0, 0.380, 0.338, 0.178],    # C
        [0.0, 0.072, 0.091, 0.033],    # R
        [0.0, 0.670, 0.602, 0.458]     # Tot Mem
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, subset in enumerate(subsets):
        ax = axes[i]
        
        # 1. Plot the full Model trajectory
        full_history = run.datacollector.model_vars[subset]
        ax.plot(range(1, 367), [val * run.global_factor for val in full_history], 
                label="Model (Full)", color='royalblue', alpha=0.6)
        
        # 2. Highlight the specific Model time-points used for fitting
        model_points = [full_history[idx] * run.global_factor for idx in indices]
        ax.scatter(time_points, model_points, color='blue', zorder=5, label="Model Points")
        
        # 3. Plot the Experimental Medians
        ax.scatter(time_points, medians[i], color='red', marker='x', s=100, 
                   linewidth=2, zorder=6, label="Exp. Median")
        
        ax.set_title(titles[i])
        ax.set_xlabel("Days")
        ax.set_ylabel("Response Magnitude (%)")
        ax.grid(True, linestyle='--', alpha=0.5)
        if i == 0:
            ax.legend()

    # Hide the empty 6th subplot
    axes[5].axis('off')
    
    plt.tight_layout()
    plt.show()