def calculate_scale_and_RMSE(model_medians):
    #Experimental Medians
    data_TSCM = [0.0, 0.218, 0.171, 0.246]
    data_TCM = [0.0, 0.380, 0.338, 0.178]
    data_TEMRA = [0.0, 0.072, 0.091, 0.034]
    # Peak values for normalization
    norm_TSCM = max(data_TSCM)  # 0.246
    norm_TCM = max(data_TCM)    # 0.380
    norm_TEMRA = max(data_TEMRA) # 0.091
    # Calculate Scaling Factor 
    numerator = 0
    denominator = 0
    for i in range(4):
    # Each term is divided by its respective max(y)^2
        numerator += (
            (model_medians['%TSCM'][i] * data_TSCM[i]) / (norm_TSCM**2) + 
            (model_medians['%TCM'][i] * data_TCM[i]) / (norm_TCM**2) + 
            (model_medians['%TEMRA'][i] * data_TEMRA[i]) / (norm_TEMRA**2)
        )

        denominator += (
            (model_medians['%TSCM'][i]**2) / (norm_TSCM**2) + 
            (model_medians['%TCM'][i]**2) / (norm_TCM**2) + 
            (model_medians['%TEMRA'][i]**2) / (norm_TEMRA**2)
        )

    s_factor = numerator / denominator if denominator != 0 else 1
    # Calculate RMSE
    total_error = 0

    for i in range(4):
        err_S = ((model_medians["%TSCM"][i] * s_factor - data_TSCM[i]) / norm_TSCM)**2
        err_C = ((model_medians["%TCM"][i] * s_factor - data_TCM[i]) / norm_TCM)**2
        err_R = ((model_medians["%TEMRA"][i] * s_factor - data_TEMRA[i]) / norm_TEMRA)**2
        total_error += (err_S + err_C + err_R)

    rmse = (total_error / 12)**0.5 # 12 points (3 subsets * 4 timepoints)

    return s_factor, rmse