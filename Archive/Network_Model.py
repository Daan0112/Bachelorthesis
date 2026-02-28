import mesa
import mesa.agent
import matplotlib.pyplot as plt
import networkx as nx
import random

###### Helper functions ######

###### Infection Agents ######
class Person_Agent(mesa.Agent):
    def __init__(self, model, stage = 'Susceptible'):
        super().__init__(model)
        self.stage = stage
    
    ######################################
    ########    Helper methods    ########
    ######################################
    
    ####################################
    ########    Main methods    ########
    ####################################
                                     
    def infect(self):
        self.stage = 'Infected'
    # Change cell state to healthy
    def recover(self):
        self.stage = 'Recovered'
    
    def infect_other(self):
        # 1. identify all neighbours
        Neighbour_agents = self.model.grid.get_neighbors(self.pos,include_center = False, radius = 1)
        # 2. determine healthy cells
        Susceptible_Neighbour_agents = [agent for agent in Neighbour_agents if agent.stage == 'Susceptible']
        for agent in Susceptible_Neighbour_agents:
            if random.random() < 0.1:
                agent.infect()
            

    def step(self):
        if self.stage == 'Infected':
            self.infect_other()
            if random.random() < 0.05:
                self.recover()

###### Model ######


class Immunology_Model(mesa.Model):
    """
    Manager class to run Model
    """
    def __init__(self, initial_population=500, network = None):
        super().__init__()
        #initiate mesa network class
        self.G = nx.connected_watts_strogatz_graph(initial_population, k=4, p=0.6)
        self.grid = mesa.space.NetworkGrid(self.G)
        #nx.draw(self.G)
        #plt.show()
        # Create agents for every person
        for node in self.G.nodes():
            a = Person_Agent(self,stage='Susceptible')
            self.grid.place_agent(a, node)

        # initiate data collector
        self.datacollector = mesa.DataCollector(
            model_reporters={
                "Susceptible_People": lambda m: sum(1 for agent in m.agents_by_type[Person_Agent] if agent.stage == 'Susceptible'),
                "Infected_People": lambda m: sum(1 for agent in m.agents_by_type[Person_Agent] if agent.stage == 'Infected'),
                "Recovered_People": lambda m: sum(1 for agent in m.agents_by_type[Person_Agent] if agent.stage == 'Recovered'),
            }
        )
        # Advance the model by one step
        self.datacollector.collect(self)
    
    def infect_patient_zero(self):
        patient_zero = random.choice(self.grid.get_all_cell_contents())
        patient_zero.infect()

    def step(self):

        for agent_class in self.agent_types:
            self.agents_by_type[agent_class].shuffle_do("step")
        # collect model level data
        self.datacollector.collect(self)

model = Immunology_Model()
model.infect_patient_zero()
# for _ in range(200):
#     model.step()

# results = model.datacollector.get_model_vars_dataframe()

# print(results)
# # Suppose your model is `model`
# G = model.G

# # Define color map based on agent state
# color_map = []
# for node in G.nodes():
#     agent = model.grid.get_cell_list_contents([node])[0]
#     if agent.stage == "Susceptible":
#         color_map.append("blue")
#     elif agent.stage == "Infected":
#         color_map.append("red")
#     elif agent.stage == "Recovered":
#         color_map.append("green")

# # Draw the graph
# pos = nx.spring_layout(G)  # or layout from your model
# nx.draw(G, pos, node_color=color_map, with_labels=False)
# plt.show()
pos = nx.spring_layout(model.G)
for i in range(200):
    model.step()
    plt.clf()
    color_map = []
    for node in model.G.nodes():
        agent = model.grid.get_cell_list_contents([node])[0]
        color_map.append(
            {"Susceptible": "blue", "Infected": "red", "Recovered": "green"}[agent.stage]
        )
    nx.draw(model.G, pos, node_color=color_map, with_labels=False)
    plt.title(f"Step {i}")
    plt.pause(0.3)