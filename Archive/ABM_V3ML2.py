import mesa
# from mesa.visualization import SolaraViz, make_plot_component
# import solara
# import math
# import numpy as np
import matplotlib.pyplot as plt
# import itertools
import random
# import pandas as pd
# import altair as alt
# from mesa.visualization.utils import update_counter
# import datetime
import copy
import optuna

def antigen(model):
    return model.steps < model.T_antigen

class CD4Cell(mesa.Agent):
    def __init__(self, model, cell_type = 'N', age = 0, path = "Undefined"):
        super().__init__(model)
        self.cell_type = cell_type  # N, S, C, R <-- Naive, TSCM, TCM, TEMRA
        self.age = age

        #### Select path on initialization. ####
        if path == "Undefined":
            self.path = self.select_path()
        else:
            self.path = path

    ####################################
    ## ## ## # Helper methods # ## ## ##
    ####################################

    def select_path(self):  # Function to define path
        if random.random() < self.model.frac_memory:
            if random.random() < self.model.S_C_Ratio:
                return "S"
            else:
                return "C"
        else:
            return "R"

    ####################################
    ## ## ## ## Main methods ## ## ## ##
    ####################################

    # Activate T-cell is detected
    def activate(self):
        # N to S/R
        if self.cell_type == "N" and antigen(self.model) and random.random() < self.model.p_act:
            self.cell_type = self.path
            self.activation_time = self.model.steps
        elif self.cell_type == "S" and random.random() < self.model.p_SC:
            self.cell_type = "C"
        elif self.cell_type == "C" and random.random() < self.model.p_CR:
            self.cell_type = "R"

    def proliferate(self):
        if antigen(self.model):    
            if self.cell_type == "N" and random.random() < self.model.r_N:
                CD4Cell(self.model, cell_type = 'N')
            elif self.cell_type == "S" and random.random() < self.model.r_S:
                CD4Cell(self.model, cell_type = 'S')
            elif self.cell_type == "C" and random.random() < self.model.r_C:
                CD4Cell(self.model, cell_type = 'C')
            elif self.cell_type == "R" and random.random() < self.model.r_R:
                CD4Cell(self.model, cell_type = 'R')

    # Age the T-cell and handle transitions or death
    def death_age(self):
        self.age += 1
        # Death based on cell type
        if self.cell_type == "N" and random.random() < self.model.d_N:
            self.remove()
            self.death_time = self.model.steps
        elif self.cell_type == "S" and random.random() < self.model.d_S:
            self.remove()
            self.death_time = self.model.steps
        elif self.cell_type == "C" and random.random() < self.model.d_C:
            self.remove()
            self.death_time = self.model.steps
        elif self.cell_type == "R" and random.random() < self.model.d_R:
            self.remove()
            self.death_time = self.model.steps

    ###################################
    ## ## ## ## Step method ## ## ## ##
    ###################################

    def step(self):
        self.activate()
        self.proliferate()
        self.death_age()

class Immunology_Model(mesa.Model):
    """
    Manager class to run Model
    """
    def __init__(self,
                 T_antigen = 23,
                 p_act=0.00022,
                 frac_memory=0.75,
                 S_C_Ratio=0.6,
                 p_SC=0.0015,
                 p_CR=0.0015,
                 r_N=0.02,
                 r_S=0.075,
                 r_C=0.0,
                 r_R=0.0,
                 d_N=0.0005,
                 d_S=0.0001,
                 d_C=0.0016,
                 d_R=0.011,
                 initial_population=5_000,
                 global_factor=1):
        super().__init__()
        # Time antigen present
        self.T_antigen = int(T_antigen)
        # Activation rates
        self.p_act = float(p_act)
        # Path fractions
        self.frac_memory = float(frac_memory)
        self.S_C_Ratio = float(S_C_Ratio)
        # Pathway rates
        self.p_SC = float(p_SC)
        self.p_CR = float(p_CR)
        # Proliferate rates
        self.r_N = float(r_N)
        self.r_S = float(r_S)
        self.r_C = float(r_C)
        self.r_R = float(r_R)
        # Death rates
        self.d_N = float(d_N)
        self.d_S = float(d_S)
        self.d_C = float(d_C)
        self.d_R = float(d_R)
        # Initial population
        self.initial_population = int(initial_population)
        self.global_factor = global_factor

        # create CD4Cells
        for i in range(self.initial_population):
           CD4Cell(self, cell_type = 'N')


        # initiate data collector
        self.datacollector = mesa.DataCollector(
            model_reporters={
                "%N": lambda m: (sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'N') / self.initial_population)*100,
                "%S": lambda m: (sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'S') / self.initial_population)*100,
                "%C": lambda m: (sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'C') / self.initial_population)*100,
                "%R": lambda m: (sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'R') / self.initial_population)*100,
                "%Total_memory": lambda m: (sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type in ['S', 'C', 'R']) / self.initial_population)*100,
                "Total_Live": lambda m: len(m.agents)
            }
)
        # Advance the model by one step
        self.datacollector.collect(self)

    def step(self):
        for agent_class in self.agent_types:
            self.agents_by_type[agent_class].shuffle_do("step")
        # collect model level data
        self.datacollector.collect(self)

def Calculate_variance2median(run):
    time_points = [0, 21, 42, 364]
    N_modelsubset = [run.datacollector.model_vars["%N"][i] for i in time_points]
    S_modelsubset = [run.datacollector.model_vars["%S"][i] for i in time_points]
    C_modelsubset = [run.datacollector.model_vars["%C"][i] for i in time_points]
    R_modelsubset = [run.datacollector.model_vars["%R"][i] for i in time_points]
    Tot_mem_modelsubset = [run.datacollector.model_vars["%Total_memory"][i] for i in time_points]

    N_datasubset_median = [0.349, 1.005, 0.039, 0.0]
    S_datasubset_median = [0.0, 0.218, 0.171, 0.246]
    C_datasubset_median = [0.0, 0.380, 0.338, 0.178]
    R_datasubset_median = [0.0, 0.072, 0.091, 0.033]
    Tot_mem_datasubset_median = [0.0, 0.670, 0.602, 0.458]
    loss_N = 0
    loss_S = 0
    loss_C = 0
    loss_R = 0
    loss_Tot_mem = 0
    for i in range(4):
        add_loss_N = (N_modelsubset[i] * run.global_factor - N_datasubset_median[i])**2
        add_loss_S = (S_modelsubset[i] * run.global_factor - S_datasubset_median[i])**2
        add_loss_C = (C_modelsubset[i] * run.global_factor - C_datasubset_median[i])**2
        add_loss_R = (R_modelsubset[i] * run.global_factor - R_datasubset_median[i])**2
        add_loss_Tot_mem = (Tot_mem_modelsubset[i] * run.global_factor - Tot_mem_datasubset_median[i])**2
        loss_N += add_loss_N
        loss_S += add_loss_S
        loss_C += add_loss_C
        loss_R += add_loss_R
        loss_Tot_mem += add_loss_Tot_mem
    return [loss_N, loss_S, loss_C, loss_R, loss_Tot_mem]

def plot_model_vs_data(run):
    time_points = [1, 22, 43, 365]
    indices = [0, 21, 42, 364]
    
    # Data labels and medians
    subsets = ["%N", "%S", "%C", "%R", "%Total_memory"]
    titles = ["Naive (N)", "TSCM (S)", "TCM (C)", "TEMRA (R)", "Total Memory"]
    medians = [
        [0.349, 1.005, 0.039, 0.0],    # N
        [0.0, 0.218, 0.171, 0.246],    # S
        [0.0, 0.380, 0.338, 0.178],    # C
        [0.0, 0.072, 0.091, 0.033],    # R
        [0.0, 0.670, 0.602, 0.458]     # Tot Mem
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, subset in enumerate(subsets):
        ax = axes[i]
        
        # 1. Plot the full Model trajectory
        full_history = run.datacollector.model_vars[subset]
        ax.plot(range(1, 367), [val * run.global_factor for val in full_history], 
                label="Model (Full)", color='royalblue', alpha=0.6)
        
        # 2. Highlight the specific Model time-points used for fitting
        model_points = [full_history[idx] * run.global_factor for idx in indices]
        ax.scatter(time_points, model_points, color='blue', zorder=5, label="Model Points")
        
        # 3. Plot the Experimental Medians
        ax.scatter(time_points, medians[i], color='red', marker='x', s=100, 
                   linewidth=2, zorder=6, label="Exp. Median")
        
        ax.set_title(titles[i])
        ax.set_xlabel("Days")
        ax.set_ylabel("Response Magnitude")
        ax.grid(True, linestyle='--', alpha=0.5)
        if i == 0:
            ax.legend()

    # Hide the empty 6th subplot
    axes[5].axis('off')
    
    plt.tight_layout()
    plt.show()

#def run_simulation(params_dict):
#    """Helper to run the model and return the run object."""
#    model = Immunology_Model(**params_dict)
#    for _ in range(365):
#        model.step()
#    return model

#base_params =  {'T_antigen': 22, 'initial_population': 1000, 'global_factor': 0.0024843846383350164, 'p_act': 0.030317289499438103, 'frac_memory': 0.8250177555065171, 'S_C_Ratio': 0.317240200260811, 'p_SC': 1.2065068520807962e-06, 'p_CR': 0.006112218789151671, 'r_N': 0.08341126901474147, 'r_S': 0.09443847460848312, 'r_C': 0.0, 'r_R': 0.0, 'd_N': 0.009108535626386806, 'd_S': 1.920910467625338e-07, 'd_C': 0.010761146990473012, 'd_R': 6.462961373287432e-05}
#run = run_simulation(base_params)
#variance = Calculate_variance2median(run)
#print(variance)
#plot_model_vs_data(run)



def objective(trial):
    # --- 1. TIME & POPULATION ---
    # suggest_int for whole numbers
    T_antigen = trial.suggest_int("T_antigen", 10, 30)
    initial_population = trial.suggest_int("initial_population", 1000, 1000)
    global_factor = trial.suggest_float("global_factor", 1e-7, 0.1, log=True)
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
        "global_factor": global_factor,
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
    model = Immunology_Model(**params)
    for _ in range(365):
        model.step()
    
    # Use your existing loss function
    losses = Calculate_variance2median(model)
    
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

