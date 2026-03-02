def Calculate_global_factor(run):
    time_points = [0, 21, 42, 364]
    N_modelsubset = [run.datacollector.model_vars["%N"][i] for i in time_points]
    S_modelsubset = [run.datacollector.model_vars["%S"][i] for i in time_points]
    C_modelsubset = [run.datacollector.model_vars["%C"][i] for i in time_points]
    R_modelsubset = [run.datacollector.model_vars["%R"][i] for i in time_points]

    N_datasubset_median = [0.349, 1.005, 0.039, 0.0]
    S_datasubset_median = [0.0, 0.218, 0.171, 0.246]
    C_datasubset_median = [0.0, 0.380, 0.338, 0.178]
    R_datasubset_median = [0.0, 0.072, 0.091, 0.033]
    numerator=0
    denominator=0
    for i in range(4):
        numerator+=(N_modelsubset[i]*N_datasubset_median[i] +
                    S_modelsubset[i]*S_datasubset_median[i]
                    + C_modelsubset[i]*C_datasubset_median[i]
                    + R_modelsubset[i]*R_datasubset_median[i])
        denominator+=(N_modelsubset[i]**2 +
                      S_modelsubset[i]**2
                      + C_modelsubset[i]**2
                      + R_modelsubset[i]**2)
    try:
        run.global_factor=numerator/denominator
    except ZeroDivisionError:
        run.global_factor=1.0
        print("Warning: Denominator is zero. Global factor set to 1.0 by default.")