
import mesa
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

#Helper functions

def get_distance(pos_1,pos_2):
  '''
  Calculate the Euclidean distance between two points
  used in trade.move()
  '''
  x1,y1 = pos_1
  x2,y2 = pos_2
  dx = x1 - x2
  dy = y1 - y2
  distance = math.sqrt(dx**2 + dy**2)
  return distance

def flatten(list_of_lists):
  '''
  helper function for model datacollector for trade price
  collapses agent prices lists into one list
  '''
  return [item for sublist in list_of_lists for item in sublist]

def geometric_mean(pricelist):
  '''
  helper function for model datacollector for trade price
  returns the geometric mean of a list
  '''
  return np.exp(np.log(pricelist).mean())

def get_trade(agent):
  '''
  For agent reporters in data collector

  filters agent and return agent trade partners 
  '''
  if type(agent) == Trader:
    return agent.trade_partners
  else: return None

# Resource classes

class Sugar(mesa.Agent):
  """
  Sugar:
  - contains an amount of sugar
  - grows 1 amount of sugar at each turn
  """

  def __init__(self, unique_id, model, pos, max_sugar):
    super().__init__(model)
    #self.pos = pos
    self.amount = max_sugar
    self.max_sugar = max_sugar
  def step(self):
    '''
    Sugar Growth function, adds one unit of sugar each step until max amount.

    '''
    self.amount = min([self.max_sugar, self.amount+1])

class Spice(mesa.Agent):
  """
  Spice:
  - contains an amount of spice
  - grows 1 amount of spice at each turn
  """

  def __init__(self, unique_id, model, pos, max_spice):
    super().__init__(model)
    #self.pos = pos
    self.amount = max_spice
    self.max_spice = max_spice

  def step(self):
    '''
    Spice Growth function, adds one unit of spice each step until max amount.

    '''
    self.amount = min([self.max_spice, self.amount+1])

#Trader class

class Trader(mesa.Agent):
  """
  Trader:
  - has a metabolism of sugar and spice
  - harvest and trade sugar and spice to survive
  """


  def __init__(self, unique_id, model, pos, moore=False, sugar=0,
              spice=0, metabolism_sugar=0, metabolism_spice=0,
              vision=0):
    super().__init__(model)
    #self.pos = pos
    self.moore = moore
    self.sugar = sugar
    self.spice = spice
    self.metabolism_sugar = metabolism_sugar
    self.metabolism_spice = metabolism_spice
    self.vision = vision
    self.prices = []
    self.trade_partners = []

  def get_sugar(self,pos):
    '''
    used in self.get_sugar_amount()
    used in self.eat()
    '''
    this_cell = self.model.grid.get_cell_list_contents(pos)
    for agent in this_cell:
      if type(agent) is Sugar:
        return agent
    return None

  def get_spice(self,pos):
    '''
    used in self.get_spice_amount()
    used in self.eat()
    '''
    this_cell = self.model.grid.get_cell_list_contents(pos)
    for agent in this_cell:
      if type(agent) is Spice:
        return agent
    return None

  def get_sugar_amount(self,pos):
    '''
    used in self.move() as part of self.calculate_welfare()
    '''

    sugar_patch = self.get_sugar(pos)
    if sugar_patch:
      return sugar_patch.amount
    return 0

  def get_spice_amount(self,pos):
    '''
    used in self.move() as part of self.calculate_welfare()
    '''

    spice_patch = self.get_spice(pos)
    if spice_patch:
      return spice_patch.amount
    return 0

  def get_trader(self,pos):
    '''
    helper function used in self.trade_with_neighbors()
    '''

    this_cell = self.model.grid.get_cell_list_contents(pos)
    for agent in this_cell:
      if isinstance(agent,Trader):
        return(agent)


  def is_occupied_by_other(self,pos):
    '''
    Helper function part 1 of self.move()
    '''
    if pos == self.pos:
      #agent can stay on it's own spot
      return False
    this_cell = self.model.grid.get_cell_list_contents(pos)
    for a in this_cell:
      if isinstance(a,Trader):
        return True
    return False

  def calculate_welfare(self, sugar, spice):
    '''
    helper function part 2 of self.move()
    '''

    #calculate total resources
    m_total = self.metabolism_sugar + self.metabolism_spice
    #Cobb-Douglas functional form starting on p.97 in Growing Artificial Societies
    return sugar**(self.metabolism_sugar/m_total) * spice**(self.metabolism_spice/m_total)

  def is_starved(self):
    '''
    Helper function for self.maybe_die()
    '''

    return (self.sugar <= 0) or (self.spice <= 0)

  def calculate_MRS(self):
    '''
    helper function for self.trade()

    determine what the trader agent needs and can give up
    '''

    return (self.spice/self.metabolism_spice) / (self.sugar/self.metabolism_sugar)

  def calculate_sell_spice_amount(self,price):
    '''
    helper function for self.maybe_sell_spice()
    '''

    if price >= 1:
      sugar = 1
      spice = int(price)
    else:
      sugar = int(1/price)
      spice = 1
    return sugar, spice

  def sell_spice(self,other,sugar,spice):
    '''
    used in self.maybe_sell_spice()

    exchanges sugar and spice between traders
    '''
    self.sugar += sugar
    other.sugar -= sugar
    self.spice -= spice
    other.spice += spice
  
  def maybe_sell_spice(self,other,price,welfare_self,welfare_other):
    '''
    helper function for self.trade()
    '''

    sugar_exchanged, spice_exchanged = self.calculate_sell_spice_amount(price)

    # Assess new sugar and spice amount - what if change did occur
    self_sugar = self.sugar + sugar_exchanged
    other_sugar = other.sugar - sugar_exchanged
    self_spice = self.spice - spice_exchanged
    other_spice = other.spice + spice_exchanged

    # double check to ensure agents have resources

    if((self_sugar<=0) or
       (other_sugar<=0) or
       (self_spice<=0) or
       (other_spice<=0)):
       return False

    # trade criteria #1 - are both agents better off?
    both_agents_better_off = (
        (welfare_self < self.calculate_welfare(self_sugar,self_spice)) and
        (welfare_other < other.calculate_welfare(other_sugar,other_spice))
        )

    # trade criteria #2 is their mrs crossing
    mrs_not_crossing = self.calculate_MRS() > other.calculate_MRS()

    if not (both_agents_better_off and mrs_not_crossing):
      return False

    # Criteria met, execute trade
    self.sell_spice(other,sugar_exchanged,spice_exchanged)
    return True


  def trade(self,other):
    '''
    helper function used in trade_with_neighbors()

    other is a trader agent object
    '''
    # sanity check to verify code is working as expected
    assert self.sugar > 0
    assert self.spice > 0
    assert other.sugar > 0
    assert other.spice > 0

    # calculate parginal rate of subsitution in Growing Artificial Societies
    mrs_self = self.calculate_MRS()
    mrs_other = other.calculate_MRS()

    #calculate each agents welfare
    welfare_self = self.calculate_welfare(self.sugar,self.spice)
    welfare_other = other.calculate_welfare(other.sugar, other.spice)

    if math.isclose(mrs_self, mrs_other, rel_tol=1e-02):
      return

    # calculate price
    price = math.sqrt(mrs_self*mrs_other)

    if mrs_self > mrs_other:
      # self is a sugar buyer, spice seller
      sold = self.maybe_sell_spice(other,price,welfare_self,welfare_other)
      # no trade - criteria not met
      if not sold:
        return
    else:
      # self is a spice buyer, sugar seller
      sold = other.maybe_sell_spice(self,price,welfare_other,welfare_self)
      # no trade - criteria not met
      if not sold:
        return

    # Capture data
    self.prices.append(price)
    self.trade_partners.append(other.unique_id)

    # Continue trading
    self.trade(other)

  def move(self):
    '''
    Function for trader agent to identify optimal move for each step in 4 parts
    1 - identify all possible moves
    2 - determine which move maximize welfare
    3 - find closest best option
    4 - move
    '''

    # 1. identify all possible moves
    neighbors = [i for i in self.model.grid.get_neighborhood(self.pos, self.moore, True, self.vision) if not self.is_occupied_by_other(i)]

    # 2. determine which move maximizes welfare

    welfares = [self.calculate_welfare(
        self.sugar + self.get_sugar_amount(pos),
        self.spice + self.get_spice_amount(pos))
     for pos in neighbors]

    # 3. Find closest best option
    # find the highest welfare in welfares
    max_welfares = max(welfares)

    # get the index of the max welfare cells
    candidate_indices = [i for i in range(len(welfares)) if math.isclose(welfares[i],max_welfares,rel_tol = 1e-02)]

    # convert index to position of those cells
    candidates = [neighbors[i] for i in candidate_indices]

    # find the closest
    min_dist = min(get_distance(self.pos,pos) for pos in candidates)

    final_candidates = [pos for pos in candidates if math.isclose(get_distance(self.pos,pos),min_dist,rel_tol=1e-02)]

    self.random.shuffle(final_candidates)

    # 4. Move Agent
    self.model.grid.move_agent(self, final_candidates[0])


  def eat(self):
    '''
    Function for agents to get local resources and consume sugar and spice
    '''
    #get sugar
    sugar_patch = self.get_sugar(self.pos)

    if sugar_patch:
      self.sugar += sugar_patch.amount
      sugar_patch.amount = 0
    self.sugar -= self.metabolism_sugar

    #get spice
    spice_patch = self.get_spice(self.pos)

    if spice_patch:
      self.spice += spice_patch.amount
      spice_patch.amount = 0
    self.spice -= self.metabolism_spice

  def maybe_die(self):
    '''
    Function to remove Traders who have consumed all their sugar or spice
    '''

    if self.is_starved():
      self.model.grid.remove_agent(self)
      self.remove()

  def trade_with_neighbors(self):
    '''
    Function for traders agents to decide who to trade with in three parts

    1- identify neighbors who can trade
    2- trade
    3- collect data
    '''

    neighbour_agents = [self.get_trader(pos) for pos in self.model.grid.get_neighborhood(self.pos,self.moore,False,self.vision) if self.is_occupied_by_other(pos)]

    if len(neighbour_agents) == 0:
      return
    # iterate through traders in neighboring cells and trade
    for agent in neighbour_agents:
      if agent:
        self.trade(agent)
    return

#Model Class

class SugarscapeG1mt(mesa.Model):
  """
  Manager class to run Sugarscape with Traders
  """


  def __init__(self, width=50,height=50, initial_population=200, endownment_min=25,
               endownment_max=50, metabolism_min=1, metabolism_max=5, vision_min=1, vision_max=5):
    super().__init__()
    #Initiate width and heigh of sugarscape
    self.width = width
    self.height = height
    #initiate population attributes
    self.initial_population = initial_population
    self.endownment_min = endownment_min
    self.endownment_max = endownment_max
    self.metabolism_min = metabolism_min
    self.metabolism_max = metabolism_max
    self.vision_min = vision_min
    self.vision_max = vision_max

    #initiate mesa grid class
    self.grid = mesa.space.MultiGrid(self.width, self.height, torus=False)

    # initiate data collector
    self.datacollector = mesa.DataCollector(
        model_reporters={
        "Trader": lambda m: len(m.agents_by_type[Trader]),
         "Trade Volume": lambda m: sum(len(a.trade_partners) for a in self.agents_by_type[Trader]),
         "Price": lambda m: geometric_mean(flatten(a.prices for a in self.agents_by_type[Trader]))},
        agent_reporters = {
        "Trade Network": lambda a: get_trade(a)
        }
    )

    for agent_class in self.agent_types:
      self.agents_by_type[agent_class].shuffle_do("step")

    #read in landscape file from supplmentary material
    sugar_distribution = np.genfromtxt("sugar-map.txt")
    spice_distribution = np.flip(sugar_distribution, 1)

    agent_id = 0
    for _,(x,y) in self.grid.coord_iter():

      max_sugar = sugar_distribution[x,y]
      if max_sugar > 0:
        sugar = Sugar(agent_id, self, (x,y), max_sugar)
        self.grid.place_agent(sugar,(x,y))
        agent_id += 1

      max_spice = spice_distribution[x,y]
      if max_spice > 0:
        spice = Spice(agent_id, self, (x,y), max_spice)
        self.grid.place_agent(spice,(x,y))
        agent_id += 1
    for i in range(self.initial_population):
      #get agent position
      x = self.random.randrange(self.width)
      y = self.random.randrange(self.height)
      #see Growing Artificial Societies p.108 for initialization
      #give agents initial endownment
      sugar = int(self.random.uniform(self.endownment_min,self.endownment_max+1))
      spice = int(self.random.uniform(self.endownment_min,self.endownment_max+1))
      #give agents initial metabolism
      metabolism_sugar = int(self.random.uniform(self.metabolism_min,self.metabolism_max+1))
      metabolism_spice = int(self.random.uniform(self.metabolism_min,self.metabolism_max+1))
      #give agents vision
      vision = int(self.random.uniform(self.vision_min,self.vision_max))
      #create Trader object
      trader = Trader(agent_id, self, (x,y), moore=False, sugar=sugar, spice=spice,
                      metabolism_sugar=metabolism_sugar, metabolism_spice=metabolism_spice, vision=vision)
      #place agent
      self.grid.place_agent(trader,(x,y))
      agent_id += 1

  def randomize_traders(self):
    '''
    helper function for self.step()

    put traders in random list
    '''
    trader_shuffle = list(self.agents_by_type[Trader])
    self.random.shuffle(trader_shuffle)

    return trader_shuffle

  def step(self):
    '''
    Step function for model
    '''
    for sugar in self.agents_by_type[Sugar]:
      sugar.step()
    for spice in self.agents_by_type[Spice]:
      spice.step()

    # Randomize agent
    trader_shuffle = self.randomize_traders()

    for agent in trader_shuffle:
      agent.prices = []
      agent.trade_partners = []
      agent.move()
      agent.eat()
      agent.maybe_die()

    trader_shuffle = self.randomize_traders()

    for agent in trader_shuffle:
      agent.trade_with_neighbors()
  
    # collect model level data
    self.datacollector.collect(self)

  def run_model(self,step_count=100):
    for i in range(step_count):
      self.step()


#Run Sugarscape


model = SugarscapeG1mt()
model.run_model()

#Analyze data

results = model.datacollector.get_model_vars_dataframe()

results

# plot number of agents per time step
results.plot(y = "Trader", use_index=True)
plt.show()
# plot trade price per step
y = list(results["Price"])
x = range(100)
plt.scatter(x,y, s=1)
plt.show()
# Plot trade volume

plt.bar(results.index[10:], results["Trade Volume"][10:])
plt.show()
# Plot trade volume improved

for i in range(100):
  plt.vlines(i,0,results["Trade Volume"][i])

agent_results = model.datacollector.get_agent_vars_dataframe()
agent_results = agent_results[agent_results["Trade Network"].notnull()]
agent_results

#create graph object
G = nx.Graph()

#add agent keys to make initial mode
G.add_nodes_from([agent.unique_id for agent in list(model.agents_by_type[Trader])])

#create edge list
for idx, row in agent_results.iterrows():
 if len(row["Trade Network"]) > 0:
   for agent in row["Trade Network"]:
     G.add_edge(idx[1], agent)   

print(nx.node_connectivity(G))
print(nx.average_clustering(G))
try:
 nx.diameter(G)
except:
 print("failed")
nx.global_efficiency(G)

degree = [d for n, d in G.degree()]
plt.hist(degree)
plt.show()
nx.draw_shell(G)
plt.show()
#Batch Run and Analysis

# params = {"width": 50, "height": 50, "vision_min": range(1,3), "metabolism_max": [3,5]}

# results_batch = mesa.batch_run(
#     SugarscapeG1mt,
#     parameters = params,
#     iterations = 1,
#     number_processes = 1,
#     data_collection_period = 1,
#     display_progress = True)

# import pandas as pd

# results_df = pd.DataFrame(results_batch)
# results_df

# plt.scatter(results_df["Step"], results_df["Price"], s=0.75)

# results_explore = results_df[results_df["metabolism_max"]==3]
# results_explore

# plt.scatter(results_explore["Step"], results_explore["Price"], s=0.75)

# for i in range(4):
#   results_explore = results_df[results_df["RunId"] == i]
#   plt.plot(results_explore["Step"],results_explore["Trader"])