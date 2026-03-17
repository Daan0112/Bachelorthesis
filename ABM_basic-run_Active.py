import ABM_model_Active as modelV4
import ABM_variance_Active as varianceV4
import ABM_plot_Active as plotV4


def run_simulation(params_dict):
    """Helper to run the model and return the run object."""
    model = modelV4.Immunology_Model(**params_dict)
    for _ in range(365):
        model.step()
    variance = varianceV4.Calculate_variance2median(model)
    print(variance)
    plotV4.plot_model_vs_data(model)
    return model


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
run_simulation(params)
