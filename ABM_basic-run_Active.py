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

base_params = {'T_antigen': 22, 'initial_population': 1000, 'p_act': 0.7130575836504838, 'frac_memory': 0.81985770893272, 'S_C_Ratio': 0.13072410353235184, 'p_SC': 3.112310238392141e-05, 'p_CR': 0.0014285840051584267, 'r_N': 0.049666033337027966, 'r_S': 0.06515909784610464, 'r_C': 0.0, 'r_R': 0.0, 'd_N': 0.0035584268482424285, 'd_S': 4.731783078114052e-06, 'd_C': 0.0007002347364540852, 'd_R': 0.017535584696647337}
run_simulation(base_params)
