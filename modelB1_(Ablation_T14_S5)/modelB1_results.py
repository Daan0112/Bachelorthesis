import optuna
import modelB1_plot as plotter
import modelB1_multrun as multrun

study = optuna.load_study(study_name="yellowfever_250seeds_q25_t14s5_ablation", storage="sqlite:///yellowfever_250seeds_q25_t14s5_ablation.db")

# Base params
params = {
    "alpha_peak": 0.01, "b_MPEC": 1, "K_mem": 999_999,
    "S_CD4": 500, "q": 0.25, "mu_N": 0.0003,
    "mu_TSCM": 0.0002, "mu_TCM": 0.004, "mu_TEM": 0.01,
    "mu_TEMRA": 0.02, "mu_MPEC": 0.02, "mu_SLEC": 0.05,
    "f_TSCM": 0.03, "f_TCM": 0.05, "f_TEM": 0.06,
    "f_TEMRA": 0.02, "t_peak": 14, "sigma": 5
}
# Update to optimized params
params.update(study.best_params)
params["b_SLEC"] = params["b_MPEC"]

m, low, high = multrun.get_plot_data(params, num_seeds=250)
plotter.plot_model_vs_data(m, low, high, params, saving="yellowfever_250seeds_q25_t14s5_ablation.png")
print(params)

