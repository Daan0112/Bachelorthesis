import mesa
import random

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
                 # global_factor=1
                 ):
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
        self.global_factor = 'Undefined'

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