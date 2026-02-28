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
                 initial_population=5_000):
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

def Calculate_variance2median(run, global_factor = 1.0):
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
        add_loss_N = (N_modelsubset[i] * global_factor - N_datasubset_median[i])**2
        add_loss_S = (S_modelsubset[i] * global_factor - S_datasubset_median[i])**2
        add_loss_C = (C_modelsubset[i] * global_factor - C_datasubset_median[i])**2
        add_loss_R = (R_modelsubset[i] * global_factor - R_datasubset_median[i])**2
        add_loss_Tot_mem = (Tot_mem_modelsubset[i] * global_factor - Tot_mem_datasubset_median[i])**2
        loss_N += add_loss_N
        loss_S += add_loss_S
        loss_C += add_loss_C
        loss_R += add_loss_R
        loss_Tot_mem += add_loss_Tot_mem
    return [loss_N, loss_S, loss_C, loss_R, loss_Tot_mem]

def plot_model_vs_data(run, global_factor=1.0):
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
        ax.plot(range(1, 367), [val * global_factor for val in full_history], 
                label="Model (Full)", color='royalblue', alpha=0.6)
        
        # 2. Highlight the specific Model time-points used for fitting
        model_points = [full_history[idx] * global_factor for idx in indices]
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

def run_simulation(params_dict):
    """Helper to run the model and return the run object."""
    model = Immunology_Model(**params_dict)
    for _ in range(365):
        model.step()
    return model

def optimize_naive_phase(base_params, global_factor, iterations=30):
    current_params = copy.deepcopy(base_params)
    
    # Get initial loss using your index 0
    current_loss = get_optimization_score(current_params, global_factor)
    
    print(f"Starting Phase 1 (Naive). Initial Loss_N: {current_loss:.6f}")

    for i in range(iterations):
        key = random.choice(NAIVE_KEYS)
        original_val = current_params[key]
        
        # Nudge the parameter
        nudge = random.uniform(0.8, 1.2)
        current_params[key] = original_val * nudge
        
        # Calculate new loss using your index 0
        new_loss = get_optimization_score(current_params, global_factor)
        
        if new_loss < current_loss:
            current_loss = new_loss
            print(f"Iter {i}: {key} improved to {current_params[key]:.6f} | Loss_N: {current_loss:.6f}")
        else:
            # Revert
            current_params[key] = original_val
            
    return current_params

def get_optimization_score(params, global_factor):
    """
    Uses YOUR function to calculate loss and returns only the 
    Naive loss (index 0) for this phase.
    """
    # 1. Run the model with current parameters
    run = Immunology_Model(**params)
    for _ in range(365):
        run.step()
    
    # 2. Call YOUR function
    # It returns [loss_N, loss_S, loss_C, loss_R, loss_Tot_mem]
    all_losses = Calculate_variance2median(run, global_factor)
    
    # 3. Return only the Naive loss (index 0)
    return all_losses[0]

base_params = {
    "T_antigen": 23,
    "p_act": 0.037253,
    "frac_memory": 0.75,
    "S_C_Ratio": 0.6,
    "p_SC": 0.0015,
    "p_CR": 0.0015,
    "r_N": 0.039907,
    "r_S": 0.075,
    "r_C": 0.0,
    "r_R": 0.0,
    "d_N": 0.000542,
    "d_S": 0.0001,
    "d_C": 0.0016,
    "d_R": 0.011,
    "initial_population": 5000,
}

run = run_simulation(base_params)
variance = Calculate_variance2median(run,0.004)
print(variance)
# To use it:
plot_model_vs_data(run, global_factor=0.004)
# Parameters for Naive fitting
NAIVE_KEYS = ["p_act", "r_N", "d_N"]
N_TARGETS = [0.349, 1.005, 0.039, 0.0]
TIME_INDICES = [0, 21, 42, 364]
DAYS = [1, 22, 43, 365]
s = 0.004
best_params_after_naive = optimize_naive_phase(base_params, s, iterations=50)
