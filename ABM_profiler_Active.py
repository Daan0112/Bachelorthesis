import cProfile
import pstats
import ABM_model_Active as modelV4

def profile_model():
    # Use your fixed parameters
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
    model = modelV4.Immunology_Model(**params)
    for _ in range(365):
        model.step()

# Run the profiler
profiler = cProfile.Profile()
profiler.enable()
profile_model()
profiler.disable()

# Print the results sorted by 'tottime' (time spent in the function itself)
stats = pstats.Stats(profiler).sort_stats('tottime')
stats.print_stats(20) # Show the top 20 bottlenecks