import ABM_scalingfactor_Active as scalingfactor

def Calculate_variance2median(run):
    time_points = [0, 21, 42, 364]
    N_modelsubset = [run.datacollector.model_vars["%N"][i] for i in time_points]
    S_modelsubset = [run.datacollector.model_vars["%S"][i] for i in time_points]
    C_modelsubset = [run.datacollector.model_vars["%C"][i] for i in time_points]
    R_modelsubset = [run.datacollector.model_vars["%R"][i] for i in time_points]

    N_datasubset_median = [0.349, 1.005, 0.039, 0.0]
    S_datasubset_median = [0.0, 0.218, 0.171, 0.246]
    C_datasubset_median = [0.0, 0.380, 0.338, 0.178]
    R_datasubset_median = [0.0, 0.072, 0.091, 0.033]

    scalingfactor.Calculate_global_factor(run)

    loss_N = 0
    loss_S = 0
    loss_C = 0
    loss_R = 0
    for i in range(4):
        add_loss_N = (N_modelsubset[i] * run.global_factor - N_datasubset_median[i])**2
        add_loss_S = (S_modelsubset[i] * run.global_factor - S_datasubset_median[i])**2
        add_loss_C = (C_modelsubset[i] * run.global_factor - C_datasubset_median[i])**2
        add_loss_R = (R_modelsubset[i] * run.global_factor - R_datasubset_median[i])**2
        loss_N += add_loss_N
        loss_S += add_loss_S
        loss_C += add_loss_C
        loss_R += add_loss_R
    return [loss_N/sum(N_datasubset_median), loss_S/sum(S_datasubset_median), loss_C/sum(C_datasubset_median), loss_R/sum(R_datasubset_median)] 