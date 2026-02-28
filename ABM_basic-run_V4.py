import ABM_model_V4 as modelV4
import ABM_variance_V4 as varianceV4
import ABM_plot_V4 as plotV4


def run_simulation(params_dict):
    """Helper to run the model and return the run object."""
    model = modelV4.Immunology_Model(**params_dict)
    for _ in range(365):
        model.step()
    variance = varianceV4.Calculate_variance2median(model)
    print(variance)
    plotV4.plot_model_vs_data(model)
    return model

base_params =  {'T_antigen': 22, 'initial_population': 1000, 'p_act': 0.01713332054406515, 'frac_memory': 0.2672530433147862, 'S_C_Ratio': 0.7791262745678639, 'p_SC': 8.130671443154906e-06, 'p_CR': 0.011627423834149758, 'r_N': 0.09530478896550659, 'r_S': 0.062412966322281306, 'r_C': 0.0, 'r_R': 0.0, 'd_N': 0.05735070697698946, 'd_S': 1.4479126825448941e-05, 'd_C': 5.578563936899805e-05, 'd_R': 0.010173207088739552}
run_simulation(base_params)
