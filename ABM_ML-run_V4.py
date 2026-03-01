import ABM_model_V4 as modelV4
import ABM_variance_V4 as varianceV4
import optuna

def objective(trial):
    # --- 1. TIME & POPULATION ---
    # suggest_int for whole numbers
    T_antigen = trial.suggest_int("T_antigen", 10, 30)
    initial_population = trial.suggest_int("initial_population", 1000, 1000)
    # global_factor = trial.suggest_float("global_factor", 1e-7, 0.1, log=True)
    # --- 2. ACTIVATION & DIFFERENTIATION ---
    # Log scale is better for p_act because 0.0001 and 0.001 are both valid
    p_act = trial.suggest_float("p_act", 1e-7, 1, log=True)
    frac_memory = trial.suggest_float("frac_memory", 0, 1)
    S_C_Ratio = trial.suggest_float("S_C_Ratio", 0, 1)
    # --- 3. TRANSITION PROBABILITIES ---
    p_SC = trial.suggest_float("p_SC", 1e-7, 1, log=True)
    p_CR = trial.suggest_float("p_CR", 1e-7, 1, log=True)
    # --- 4. PROLIFERATION RATES (r) ---
    r_N = trial.suggest_float("r_N", 0, 0.1)
    r_S = trial.suggest_float("r_S", 0, 0.1)
    r_C = trial.suggest_float("r_C", 0, 0.0)
    r_R = trial.suggest_float("r_R", 0, 0.0)
    # --- 5. DEATH RATES (d) ---
    d_N = trial.suggest_float("d_N", 1e-7, 0.1, log=True)
    d_S = trial.suggest_float("d_S", 1e-7, 0.1, log=True)
    d_C = trial.suggest_float("d_C", 1e-7, 0.1, log=True)
    d_R = trial.suggest_float("d_R", 1e-7, 0.1, log=True)

    # --- 6. PACKING INTO DICTIONARY ---
    params = {
        "T_antigen": T_antigen,
        "initial_population": initial_population,
        # "global_factor": global_factor,
        "p_act": p_act,
        "frac_memory": frac_memory,
        "S_C_Ratio": S_C_Ratio,
        "p_SC": p_SC,
        "p_CR": p_CR,
        "r_N": r_N,
        "r_S": r_S,
        "r_C": r_C,
        "r_R": r_R,
        "d_N": d_N,
        "d_S": d_S,
        "d_C": d_C,
        "d_R": d_R
    }

    # --- 7. RUN & EVALUATE ---
    model = modelV4.Immunology_Model(**params)
    for _ in range(365):
        model.step()
    
    # Use your existing loss function
    losses = varianceV4.Calculate_variance2median(model)
    
    # Return the sum of all losses to the optimizer
    return sum(losses)

# 4. Start the "Machine Learning" Search
study = optuna.create_study(
    study_name="typhoid_fit", 
    storage="sqlite:///my_study.db", 
    load_if_exists=True,
    direction="minimize"
)
study.optimize(objective, n_trials=1000) 

print("Best parameters found:", study.best_params)
