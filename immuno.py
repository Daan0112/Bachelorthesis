
# # Set-up environment

# %%
import mesa
from mesa.visualization import SolaraViz, make_plot_component
import solara
import math
import numpy as np
import matplotlib.pyplot as plt
import itertools
import random
import pandas as pd
import altair as alt
from mesa.visualization.utils import update_counter
import datetime
from mesa.discrete_space import OrthogonalMooreGrid

# %% [markdown]
# # Help Functions

# %%
def antigen(model):
    return model.steps < model.T_antigen

# %% [markdown]
# # CD4Cell
# ######(formerly called: T-cell agent)

# %%
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
        self.death_age()

# %% [markdown]
# # Model

# %%
class Immunology_Model(mesa.Model):
    """
    Manager class to run Model
    """
    def __init__(self,
                 T_antigen = 23,
                 p_act=0.04,
                 frac_memory=0.75,
                 S_C_Ratio=0.6,
                 d_N=0.0005,
                 d_S=0.0001,
                 d_C=0.0016,
                 d_R=0.011,
                 initial_population=1_000):
        super().__init__()
        self.T_antigen = int(T_antigen)
        self.p_act = float(p_act)
        self.frac_memory = float(frac_memory)
        self.S_C_Ratio = float(S_C_Ratio)
        self.d_N = float(d_N)
        self.d_S = float(d_S)
        self.d_C = float(d_C)
        self.d_R = float(d_R)
        self.initial_population = int(initial_population)

        # create CD4Cells
        for i in range(self.initial_population):
           CD4Cell(self, cell_type = 'N')


        # initiate data collector
        self.datacollector = mesa.DataCollector(
            model_reporters={
                "N": lambda m: sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'N') / max(1, len(m.agents)),
                "S": lambda m: sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'S') / max(1, len(m.agents)),
                "C": lambda m: sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'C') / max(1, len(m.agents)),
                "R": lambda m: sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type == 'R') / max(1, len(m.agents)),
                "Total-memory": lambda m: sum(1 for agent in m.agents_by_type[CD4Cell] if agent.cell_type in ['S', 'C', 'R']) / max(1, len(m.agents)),
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

# %% [markdown]
# # Set parameters

# %%
T_antigen   = 23        # duration of antigen presence (days)

p_act       = 0.04      # fraction of N activated per day when antigen present
frac_memory = 0.75      # fraction of activated cells that follow the memory path (S/C)
                        # (1 - frac_memory) follow the TEMRA path directly (R)

S_C_Ratio = 0.6         # fraction of memory frac going to S


# Death rates (CD4-specific; memory is long-lived)
d_N         = 0.0005    # naive death
d_S         = 0.0001    # TSCM: very slow decay
d_C         = 0.0016    # TCM: slow decay
d_R         = 0.011     # TEMRA: faster decay

# %% [markdown]
# # Solara Set-up

# %%
model_params = {
    "initial_population": {
        "type": "InputText",       # Use InputInt for whole numbers
        "value": 1000,
        "label": "Initial Agents",
    },
    "T_antigen": {
        "type": "InputText",
        "value": 23,
        "label": "Antigen Duration (Days)",
    },
    "p_act": {
        "type": "InputText",     # Use InputFloat for decimal rates
        "value": 0.04,
        "label": "Activation Prob (p_act)",
        "step": 0.01,             # Optional: adds +/- buttons next to the box
    },
    "frac_memory": {
        "type": "InputText",
        "value": 0.75,
        "label": "Memory Fraction",
    },
    "S_C_Ratio": {
        "type": "InputText",
        "value": 0.6,
        "label": "S/C Ratio",
    },
    # Death Rates (Daily Decay)
    "d_N": {
        "type": "InputText",
        "value": 0.0005,
        "label": "Death: Naive",
        "step": 0.0001
    },
    "d_S": {
        "type": "InputText",
        "value": 0.0001,
        "label": "Death: TSCM",
        "step": 0.0001
    },
    "d_C": {
        "type": "InputText",
        "value": 0.0016,
        "label": "Death: TCM",
        "step": 0.0001
    },
    "d_R": {
        "type": "InputText",
        "value": 0.011,
        "label": "Death: TEMRA",
        "step": 0.001
    },
}

# %%




@solara.component
def SingleCellChart(model, cell_type, median_data=None):
    # This ensures the chart updates when you click 'Step'
    update_counter.get() 
    
    # Get the data and prepare for Altair
    df = model.datacollector.get_model_vars_dataframe().reset_index()
    df = df.rename(columns={"index": "Step"})

    if df.empty or cell_type not in df.columns:
        return solara.Markdown(f"Waiting for **{cell_type}** data...")

    # 1. Simulation Line
    line = alt.Chart(df).mark_line(strokeWidth=3).encode(
        x=alt.X("Step:Q", title="Time (Days)"),
        y=alt.Y(f"{cell_type}:Q", title="Fraction", scale=alt.Scale(domain=[0, 1])),
        color=alt.value("steelblue")
    )

    chart = line

    # 2. Add Experimental Medians if provided
    if median_data:
        exp_df = pd.DataFrame({
            "Step": median_data['time_points'],
            "Value": median_data['values']
        })
        points = alt.Chart(exp_df).mark_point(size=80, filled=True, color='red').encode(
            x="Step:Q",
            y="Value:Q"
        )
        chart = line + points

    return solara.FigureAltair(chart.properties(width='container', height=250).interactive())






@solara.component
def ImmunologyDashboard(model):

    if isinstance(model, type):
        return solara.Markdown("Initializing model instance...")
    # Your real medians from the dataset
    medians = {
        'time_points': [1, 22, 43, 365],
        'S': [0.0, 0.218, 0.171, 0.246],
        'C': [0.0, 0.380, 0.338, 0.178],
        'R': [0.0, 0.072, 0.091, 0.033],
        'Total-memory': [0.0, 0.670, 0.602, 0.458]
    }

    with solara.Column():
        solara.Markdown("## CD4+ T-Cell Differentiation Model")
        
        # Row 1: TSCM and TCM
        with solara.Row():
            with solara.Card("Stem Cell Memory (TSCM)"):
                SingleCellChart(model, "S", {"time_points": medians['time_points'], "values": medians['S']})
            with solara.Card("Central Memory (TCM)"):
                SingleCellChart(model, "C", {"time_points": medians['time_points'], "values": medians['C']})
        
        # Row 2: TEMRA and Total Memory
        with solara.Row():
            with solara.Card("Effector Memory (TEMRA)"):
                SingleCellChart(model, "R", {"time_points": medians['time_points'], "values": medians['R']})
            with solara.Card("Total Memory Pool"):
                SingleCellChart(model, "Total-memory", {"time_points": medians['time_points'], "values": medians['Total-memory']})

        # Data Export Section with Metadata
        with solara.Card("Data Export"):
            def get_csv():
                # 1. Grab the current data as a dictionary
                data_dict = model.datacollector.model_vars
                
                # 2. Convert to DataFrame FIRST (creates a separate copy)
                df_export = pd.DataFrame(data_dict)
                
                # 3. Add parameters to the EXPORT copy only
                params = ["p_act", "frac_memory", "S_C_Ratio", "d_N", "d_S", "d_C", "d_R", "initial_population", "T_antigen"]
                for p in params:
                    df_export[f"param_{p}"] = getattr(model, p, "N/A")
                
                df_export["runtime_days"] = int(df_export.index[-1]) if not df_export.empty else 0
                
                return df_export.to_csv(index=True)
            
            import datetime
            ts = datetime.datetime.now().strftime("%H%M%S")
            solara.FileDownload(get_csv, filename=f"immuno_results_{ts}.csv", label="Download Results + Params")


# %% [markdown]
# # Solara (Main)

# %%
# Create initial model instance
model = Immunology_Model()
page = SolaraViz(
    model,
    model_params=model_params,
    components=[ImmunologyDashboard],
    name="CD4 Cell Model",
)
page


