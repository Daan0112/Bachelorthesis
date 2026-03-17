import mesa
import random
import math

def gaussian_pulse(model):
    return model.alpha_peak * math.exp( -0.5 * (( model.steps - model.t_peak ) / model.sigma)**2)

class CD4Cell(mesa.Agent):
    def __init__(self, model, cell_type = 'Naive'):
        super().__init__(model)
        self.cell_type = cell_type  # Naive, TSCM, TCM, TEMRA, TEM SLEC, MPEC 

    def cap(self):
        memory_count = self.model.counts["TSCM"] + self.model.counts["TCM"]
        occupancy = memory_count / self.model.K_mem
        cap = max(0, 1 - occupancy)
        return cap

    # Activate T-cell is detected
    def activate(self):
        random_value = random.random()
        # N to S/R
        if self.cell_type == "Naive" and random_value < gaussian_pulse(self.model):
            self.model.counts["Naive"] -= 1
            if random.random() < self.model.q:
                #SLEC Route
                self.cell_type = "SLEC"
                for _ in range(2**self.model.b_SLEC - 1):
                    CD4Cell(self.model, cell_type = "SLEC")
                self.model.counts["SLEC"] += 2**self.model.b_SLEC
            else:
                #MPEC Route
                self.cell_type = "MPEC"
                for _ in range(2**self.model.b_MPEC - 1):
                    CD4Cell(self.model, cell_type = "MPEC")
                self.model.counts["MPEC"] += 2**self.model.b_MPEC
        
        # Add mpec route
        elif self.cell_type == "MPEC" and random_value < self.model.f_TSCM*self.cap():
            self.model.counts["MPEC"] -= 1
            self.cell_type = "TSCM"
            self.model.counts["TSCM"] += 1
        # The random value - f_TSCM is to not influence chances due to cells already changed to TSCM.
        elif self.cell_type == "MPEC" and random_value - self.model.f_TSCM*self.cap() < self.model.f_TCM*self.cap():
            self.model.counts["MPEC"] -= 1
            self.cell_type = "TCM"
            self.model.counts["TCM"] += 1
        
        # Add Slec route
        elif self.cell_type == "SLEC" and random_value < self.model.f_TEMRA:
            self.model.counts["SLEC"] -= 1
            self.cell_type = "TEMRA"
            self.model.counts["TEMRA"] += 1
        # The random value - f_TSCM is to not influence chances due to cells already changed to TSCM.
        elif self.cell_type == "SLEC" and random_value - self.model.f_TSCM < self.model.f_TEM:
            self.model.counts["SLEC"] -= 1
            self.cell_type = "TEM"
            self.model.counts["TEM"] += 1

    # Age the T-cell and handle transitions or death
    def death_age(self):
        # Death based on cell type
        if self.cell_type == "Naive" and random.random() < self.model.mu_N:
            self.model.counts["Naive"] -= 1
            self.remove()
        elif self.cell_type == "TSCM" and random.random() < self.model.mu_TSCM:
            self.model.counts["TSCM"] -= 1
            self.remove()
        elif self.cell_type == "TCM" and random.random() < self.model.mu_TCM:
            self.model.counts["TCM"] -= 1
            self.remove()
        elif self.cell_type == "TEMRA" and random.random() < self.model.mu_TEMRA:
            self.model.counts["TEMRA"] -= 1
            self.remove()
        elif self.cell_type == "TEM" and random.random() < self.model.mu_TEM:
            self.model.counts["TEM"] -= 1
            self.remove()
        elif self.cell_type == "SLEC" and random.random() < self.model.mu_SLEC:
            self.model.counts["SLEC"] -= 1
            self.remove()
        elif self.cell_type == "MPEC" and random.random() < self.model.mu_MPEC:
            self.model.counts["MPEC"] -= 1
            self.remove()

    ###################################
    ## ## ## ## Step method ## ## ## ##
    ###################################

    def step(self):
        self.activate()
        self.death_age()

class Immunology_Model(mesa.Model):
    """
    Manager class to run Model
    """
    def __init__(self,
                 mu_N = 0.0003,     # FIXED
                 mu_TSCM = 0.0002,  # FIXED
                 mu_TCM = 0.004,    # FIXED
                 mu_TEM = 0.01,     # FIXED
                 mu_TEMRA = 0.02,   # FIXED
                 mu_MPEC = 0.02,    # FIXED
                 mu_SLEC = 0.05,    # FIXED
                 f_TSCM = 0.03,     # FIXED
                 f_TCM = 0.05,      # FIXED
                 f_TEM = 0.06,      # FIXED
                 f_TEMRA = 0.02,    # FIXED
                 t_peak = 18,       # FIXED
                 sigma = 7,         # FIXED
                 q = 0.35,          # FIXED
                 alpha_peak = 0.01, # FREE [0.01-0.3]
                 b_MPEC = 1,        # FREE [1-5]
                 b_SLEC = 1,        # FIXED to b_MPEC 
                 K_mem = 50,        # FREE [50-500]
                 S_CD4=500,         # FREE [500-10000]
                 ):
        super().__init__()
        self.mu_N = mu_N
        self.mu_TSCM = mu_TSCM
        self.mu_TCM = mu_TCM
        self.mu_TEM = mu_TEM
        self.mu_TEMRA = mu_TEMRA
        self.mu_MPEC = mu_MPEC
        self.mu_SLEC = mu_SLEC
        self.f_TSCM = f_TSCM
        self.f_TCM = f_TCM
        self.f_TEM = f_TEM
        self.f_TEMRA = f_TEMRA
        self.t_peak = t_peak
        self.sigma = sigma
        self.q = q
        self.alpha_peak = alpha_peak
        self.b_MPEC = b_MPEC
        self.b_SLEC = b_SLEC
        self.K_mem = K_mem
        self.S_CD4 = S_CD4

        self.counts = {
            "Naive": S_CD4, "TSCM": 0, "TCM": 0, "TEM": 0, 
            "TEMRA": 0, "SLEC": 0, "MPEC": 0
        }

        # create CD4Cells
        for _ in range(self.S_CD4):
           CD4Cell(self, cell_type = 'Naive')

        
        # initiate data collector
        self.datacollector = mesa.DataCollector(
            model_reporters={
                "%Naive": lambda m: (m.counts["Naive"] / self.S_CD4)*100,
                "%TSCM": lambda m: (m.counts["TSCM"] / self.S_CD4)*100,
                "%TCM": lambda m: (m.counts["TCM"] / self.S_CD4)*100,
                "%TEMRA": lambda m: (m.counts["TEMRA"] / self.S_CD4)*100,
                "%TEM": lambda m: (m.counts["TEM"] / self.S_CD4)*100,
                "%SLEC": lambda m: (m.counts["SLEC"] / self.S_CD4)*100,
                "%MPEC": lambda m: (m.counts["MPEC"] / self.S_CD4)*100,
                "%Total_memory": lambda m: ((m.counts["TSCM"]+m.counts["TCM"]+m.counts["TEMRA"]+m.counts["TEM"]) / self.S_CD4)*100,
                "Total_Live": lambda m: sum(m.counts.values())
            }
        )
        # Advance the model by one step
        self.datacollector.collect(self)

    def step(self):
        for agent_class in self.agent_types:
            self.agents_by_type[agent_class].shuffle_do("step")
        # collect model level data
        self.datacollector.collect(self)

if __name__=="__main__":
    run = Immunology_Model()
    for _ in range(100):
        run.step()