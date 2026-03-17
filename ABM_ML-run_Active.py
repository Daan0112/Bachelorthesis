import ABM_multrun_Active as multrun
import optuna

def objective(trial):
    # FREE parameters
    alpha_peak = trial.suggest_float("alpha_peak",0.01,0.3)
    b_MPEC = trial.suggest_int("b_MPEC",1,4) # 5
    K_mem = trial.suggest_int("K_mem",50,500)
    S_CD4 = trial.suggest_int("S_CD4",500,10000)
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

    
    return loss

# 4. Start the "Machine Learning" Search
study = optuna.create_study(
    study_name="typhoid_fit", 
    storage="sqlite:///my_study.db", 
    load_if_exists=True,
    direction="minimize"
)
study.optimize(objective, n_trials=100) 

print("Best parameters found:", study.best_params)
