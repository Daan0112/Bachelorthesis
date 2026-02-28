import ABM_model_V4 as modelV4
import ABM_variance_V4 as varianceV4
import ABM_plot_V4 as plotV4
import optuna

study = optuna.create_study(
    study_name="typhoid_fit", 
    storage="sqlite:///my_study.db",
    load_if_exists=True
)


model = modelV4.Immunology_Model(**study.best_params)
for _ in range(365):
    model.step()
variance = varianceV4.Calculate_variance2median(model)
print(variance)
print(f'total loss = {sum(variance)}')
plotV4.plot_model_vs_data(model)