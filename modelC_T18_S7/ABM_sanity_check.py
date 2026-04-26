import ABM_multrun_Active as multrun
import numpy as np

def run_check(custom_params):
    # Base parameters for comparison
    base_params = {
        "mu_N": 0.0003, "mu_TSCM": 0.0002, "mu_TCM": 0.004, "mu_TEM": 0.01,
        "mu_TEMRA": 0.02, "mu_MPEC": 0.02, "mu_SLEC": 0.05, "f_TSCM": 0.03,
        "f_TCM": 0.05, "f_TEM": 0.06, "f_TEMRA": 0.02, "t_peak": 14,
        "sigma": 5, "q": 0.35, "alpha_peak": 0.2, "b_MPEC": 2, 
        "b_SLEC": 2, "K_mem": 250, "S_CD4": 1000
    }
    # Update with specific check parameters
    base_params.update(custom_params)
    results = multrun.multrun(base_params, num_seeds=10)
    medians = multrun.calculate_multrun_medians(results)

    return medians

# 1. Set alpha_peak = 0: All pools stay at zero
m1 = run_check({"alpha_peak": 0})
print("--- Check 1 ('alpha_peak': 0) ---")
print("Max Values:", {k: float(round(max(v), 4)) for k, v in m1.items()})
if max(m1['%TSCM']) == 0 and max(m1['%TCM']) == 0 and max(m1['%TEMRA']) == 0 and max(m1['%TEM']) == 0.0:
    print('SUCCESS',"\n")
else:
    print('FAILURE',"\n")

# 2. Set q = 0 (All MPEC): TEM/TEMRA get NO signal
m2 = run_check({"q": 0})
print("--- Check 2 ('q': 0) ---")
print("Max Values:", {k: float(round(max(v), 4)) for k, v in m2.items()})
if max(m2['%TEMRA']) == 0 and max(m2['%TEM']) == 0.0:
    print('SUCCESS',"\n")
else:
    print('FAILURE',"\n")

# 3. Set q = 1 (All SLEC): TSCM/TCM get NO signal
m3 = run_check({"q": 1})
print("--- Check 3 ('q': 1) ---")
print("Max Values:", {k: float(round(max(v), 4)) for k, v in m3.items()})
if max(m3['%TSCM']) == 0 and max(m3['%TCM']) == 0:
    print('SUCCESS',"\n")
else:
    print('FAILURE',"\n")

# 4. Set all death rates = 0: Pools grow monotonically
m4 = run_check({
    "mu_N": 0, "mu_TSCM": 0, "mu_TCM": 0, "mu_TEM": 0, 
    "mu_TEMRA": 0, "mu_MPEC": 0, "mu_SLEC": 0
})
print("--- Check 4 ('mu_@': 0) ---")
total_end = m4['Total_Live'][-1]
total_start = m4['Total_Live'][0]
total_max = max(m4['Total_Live'])
print(f"Start Population: {total_start} | End Population: {total_end} | Max Population: {total_max}")
if total_end > total_start and total_max == total_end:
    print('SUCCESS',"\n")
else:
    print('FAILURE',"\n")


# 5. 365 days
m5 = run_check({})
print("--- Check 5 (Check all data for lowest value) ---")
min_val = min([min(v) for v in m5.values()])
print(f"Minimum value found: {min_val} (Should be >= 0)")
if min_val >= 0:
    print('SUCCESS',"\n")
else:
    print('FAILURE',"\n")

# 6. Plot normal
import matplotlib.pyplot as plt
import numpy as np

def plot_sanity_check(medians):
    """
    Plots the median results from run_check for a visual sanity confirmation.
    """
    plt.figure(figsize=(10, 6))
    
    # Use the days range from your results (typically 366 points)
    days = np.arange(len(next(iter(medians.values()))))
    
    # Plotting relevant subsets for CD4 dynamics
    plt.plot(days, medians['%Naive'], label='Naive', color='gray', linestyle='--')
    plt.plot(days, medians['%TSCM'], label='TSCM', color='blue')
    plt.plot(days, medians['%TCM'], label='TCM', color='green')
    plt.plot(days, medians['%TEMRA'], label='TEMRA', color='red')
    plt.plot(days, medians['%TEM'], label='TEM', color='yellow')
    plt.plot(days, medians['%SLEC'], label='SLEC', color='orange')
    plt.plot(days, medians['%MPEC'], label='MPEC', color='cyan')
    plt.plot(days, medians['%Total_memory'], label='Total Memory', color='black', linewidth=2)

    plt.title("Sanity Check 6: Raw Median Dynamics")
    plt.xlabel("Days")
    plt.ylabel("Percentage (%)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

plot_sanity_check(m5)
