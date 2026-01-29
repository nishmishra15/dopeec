# HyperParameter Tuning .py file

'''Import libraries'''
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sp

"""
Load input file with required data
"""
days20 = pd.read_excel('file location')
home = pd.read_excel('file location')

'''
Define initial conditions and parameters
1. Read variables from the input file
2. Define initial params, control volume parameters
3. Functions to create mass and energy balance models for indoor environment
'''
'''
# Function to estimate pollutant concentration without filter
# Define symbolic variables
'''
t = sp.symbols('t')
Cin = sp.Function('Cin')(t)
initial_c, ach, emi_pm, dr = sp.symbols('initial_c ach emi_pm dr')
# Define the differential equation
# ach(/h), ir(ug/s), V(m^3), emi_pm(ug/s), gamma(/s), Cin(ug/m3), t(s)
eqn = sp.Eq(sp.Derivative(Cin, t), (ir * ach) / (V * 0.5) - (Cin * ach / 3600) + (emi_pm / V) - (dr * Cin))
# Define the initial condition
initial_condition = {Cin.subs(t, 0): initial_c}
# Solve the differential equation with the initial condition
solution = sp.dsolve(eqn, ics=initial_condition)
C_indoor = solution.rhs                                      # Indoor conc in ug/m^3

'''
# Function to estimate air parameters
# T(Celsius), RH(%)
'''
def air_parameters(T, RH):
  sat_press = (np.exp(34.494 - (4924.99 / (T+ 237.1)))) / (T + 105) ** 1.57   # Saturation vapour pressure of indoor air at any time (in Pascal)
  vap_press = RH * sat_press/ 100    # Vapour pressure of indoor air at any time (in Pascal)
  w = 621.9907 * vap_press/ P_atm   # Fraction of water in indoor air at any time t (in kg H2O/kg of air)
  mixratio = 621.9907 * vap_press / (P_atm - vap_press) # Mixratio of indoor air at any time (gm H2O/kg dry air)
  h= T* (1.01 + 0.00189 * mixratio) + 2.5 * mixratio    # Enthalpy of indoor air at any time (kJ/kg)
  return sat_press, vap_press, w, mixratio, h

'''
# Estimate power consumption in Watts
# ach_dvs(/h), ach_hvac(/h), h_net(kJ/kg), h_hv(kJ/kg), power(Watts)
# Input: ACH of DVS and HVAC, Enthaply of mixed air at t-1, Enthalpy of supply air at t
# Output: HVAC, DVS, and Total Power Consumption
'''
def power_consumption(achdvs,achhvac,h_net_prev,h_hv_current):
  q_hvac = achhvac * V / 3600              # m^3/s
  q_dvs = achdvs * V / 3600                # m^3/s
  power_hv = (1000 * (((q_hvac + q_dvs) * rho_air * (h_net_prev - h_hv_current)) / (EER)) \
                            + dP_dvs * (q_hvac + q_dvs) / eff_dvs)                      
  power_dvs = (dP_dvs * achdvs * V / (3600 * eff_dvs))                                 
  power_tot = (power_hv + power_dvs)                                                    
  return power_hv, power_dvs, power_tot

'''
# Function for estimating mixed air (dvs air and indoor air being sent to HVAC)
# Returns air properties of mixed air (dvs+indoor)
'''
def mixed_air(achdvs,achhvac,hout,hin,win,wout):
  q_hvac = achhvac * V / 3600               # m^3/s
  q_dvs = achdvs * V / 3600                 # m^3/s
  h_mixed = (q_dvs * tstep * 60 * hout + q_hvac * tstep * 60 * hin)/ ((q_dvs + q_hvac) * tstep * 60)     # (in kJ/kg)
  w_mixed = (tstep * 60 * q_hvac * win + tstep * 60 * q_dvs * wout) / ((q_dvs + q_hvac) * tstep * 60)    # (in kg H2O/kg of air)
  vap_press_mixed = w_mixed * P_atm / 621.9907                                                           # (in Pascal)
  mixratio_mixed = 621.9907 * vap_press_mixed / (P_atm - vap_press_mixed)                                # (gm H2O/kg dry air)
  temp_mixed = (h_mixed - 2.5 * mixratio_mixed)/ (1.01 + 0.00189 * mixratio_mixed)                       # in Celsius
  sat_press_mixed = (np.exp(34.494 - (4924.99 / (temp_mixed + 237.1)))) / (temp_mixed + 105) ** 1.57     # (in Pascal)
  rh_mixed = 100 * vap_press_mixed / sat_press_mixed                                                     # in %
  return h_mixed, w_mixed, vap_press_mixed, mixratio_mixed, temp_mixed, sat_press_mixed, rh_mixed

'''
# HVAC Operation based on indoor air temperature and mixed air cooling only (no consideration of Power)
# Input: Mixed air temp and RH at any time, Indoor temp at any time
# Output: Supply air temp properties at next time
'''
def hvac_operation(temp_net,temp_in,rh_net):
  if (temp_net - 10 < 15):        # mixed air temp at t
        temp_hv = 15              # Temp of HVAC air at t+1
  else:
        temp_hv = temp_net - 10    # Temp of HVAC air at t+1
  if (temp_in >= T_set + 0.5):    # Indoor temp at t >=Set Temp+0.5
        rh_hv = RH_hvac
        sat_press_hv,vap_press_hv,w_hv,mixratio_hv,h_hv=air_parameters(temp_hv,rh_hv)
  else:
        temp_hv = temp_net
        rh_hv = rh_net
        sat_press_hv,vap_press_hv,w_hv,mixratio_hv,h_hv=air_parameters(temp_hv, rh_hv)
  return temp_hv, rh_hv, sat_press_hv, vap_press_hv, w_hv, mixratio_hv, h_hv

'''
Air properties at next time (t+1)
'''
def air_properties_next_time(achdvs,achhvac,w_hv,win_prev,h_hv,hin_prev,Tout_prev,Tin_prev):
  # Calculate w_bm for next time step
  q_hvac = achhvac * V / 3600
  q_dvs = achdvs * V / 3600
  win = win_prev + ((tstep * 60 * (q_hvac + q_dvs) * rho_air * (w_hv - win_prev)) / (V * rho_air))     # (in kg H2O/kg of air)
  vap_pressin = win * P_atm / 621.9907                                                                 # Pascals    
  hin = ((q_hvac + q_dvs) * tstep * 60 * h_hv/ V) + \
            ((V - (q_hvac + q_dvs) * tstep * 60) * hin_prev/ V) + \
            (alpha * tstep * 60 * (Tout_prev - Tin_prev) / (V * rho_air))                              # kJ/kg
  mixratioin = 621.9907 * vap_pressin / (P_atm - vap_pressin)                                          # (gm H2O/kg dry air)
  Tin = (hin - 2.5 * mixratioin) / (1.01 + 0.00189 * mixratioin)                                       # Celsius
  sat_pressin = (np.exp(34.494 - (4924.99 / (Tin + 237.1)))) / (Tin + 105) ** 1.57                    # Pascals
  RHin = 100 * vap_pressin / sat_pressin                                                              # %

  return win,vap_pressin,hin,mixratioin,Tin,sat_pressin,RHin
    
'''
Update supply air parameters if HVAC Power is more than Max HVAC Power
'''
def update_parameters(achdvs,achhvac,Power_hv,h_mixed):
  q_hvac = achhvac * V / 3600
  q_dvs = achdvs * V / 3600
  Power_hv = Max_power_hvac
  h_hv = h_mixed - (Power_hv - dP_dvs * (q_hvac + q_dvs) / eff_dvs) * EER / (1000 * (q_hvac + q_dvs) * rho_air)
  # Define the symbolic variable T
  T = sp.symbols('T')
  # Define the expressions
  var = (621.9907 * RH_hvac * (sp.exp(34.494 - (4924.99 / (T + 237.1)))) / (100 * (T + 105) ** 1.57))
  var1 = ((RH_hvac * sp.exp(34.494 - 4924.99 / (T + 237.1))) / (100 * (T + 105) ** 1.57))
  h = T * (1.01 + 0.00189 * var / (P_atm - var1)) + 2.5 * var / (P_atm - var1)
  # Define the equation to solve
  eqn = h_hv - h
  # Convert the equation to a callable function
  eqn_func = sp.lambdify(T, eqn)
  # Solve the equation using fsolve
  initial_guess = 25  # Initial guess for T
  T_solution = fsolve(eqn_func, initial_guess)[0]
  # Update T_hvac for the next time step
  temp_hv = T_solution
  sat_press_hv,vap_press_hv,w_hv,mixratio_hv,h_hv=air_parameters(temp_hv,RH_hvac)
  return temp_hv, sat_press_hv, vap_press_hv, w_hv, mixratio_hv, h_hv, Power_hv

"""
Create a separate Python environment and install all dependencies such as 
pytorch, gymnasium, and other relevant libraries
"""

import gymnasium as gym
from gymnasium import spaces
import random
from collections import deque
import torch
import torch.nn as nn
import torch.optim as optim
from itertools import product
import sys
import copy

class IAQEnv(gym.Env):
    def __init__(self,emi_pm_matrix,pm_dr_matrix,temp_outdoor,rh_outdoor,pm_outdoor):
        self.max_ventilation_rate = 10    # Maximum ventilation rate
        self.min_ventilation_rate = 0.5     # Minimum ventilation rate
        self.action_space_intervals = int((self.max_ventilation_rate - self.min_ventilation_rate) / 0.1)

        # Discretized action space: select ventilation rate from 0 to 10 at intervals of 0.1
        self.action_space = spaces.Discrete(self.action_space_intervals + 1)

        # Observation space: 
        # Outdoor: PM (0), Temp (1), RH (2)
        # Indoor: PM (3), Current Ventilation Rate (4), Indoor Temp (5), Indoor RH (6), Total Power (7)
        self.observation_space = spaces.Box(low=np.array([0, 0, 0, 0, self.min_ventilation_rate,0,0,0]),
                                         high=np.array([np.inf, np.inf, 100, np.inf, self.max_ventilation_rate, np.inf, 100, np.inf]),
                                            dtype=np.float32)

        # Store the PM, outdoor temp and rh matrix
        self.current_time_step = 0
        self.emi_pm_matrix = emi_pm_matrix
        self.pm_dr_matrix = pm_dr_matrix
        self.temp_outdoor = temp_outdoor
        self.rh_outdoor = rh_outdoor
        self.pm_outdoor = pm_outdoor
        
        self.inital_temp_out = self.temp_outdoor[0]
        self.inital_rh_out = self.rh_outdoor[0]
        self.initial_pm_out = self.pm_outdoor[0]
        
        self.initial_pm_in = 1.28            # Initial IAQ (ug/m3)
        self.initial_ventilation_rate = 0.5  # Initial ventilation rate (per hour)
        self.initial_temp_in = 25            # Intial indoor temperature (C)
        self.initial_rh_in = 60              # Intial relative humidity (%)
        self.intitial_power_total = 45                # in Watts
        
        self.state = np.array([self.initial_pm_out, self.inital_temp_out, self.inital_rh_out, self.initial_pm_in, \
                     self.initial_ventilation_rate, self.initial_temp_in, self.initial_rh_in, self.intitial_power_total],dtype=np.float32)

    def step(self, action):
        # Map the action to the actual ventilation rate
        new_ventilation_rate = (action / self.action_space_intervals) * \
                               (self.max_ventilation_rate - self.min_ventilation_rate) + self.min_ventilation_rate

        # Update the ventilation rate based on the action
        new_ventilation_rate = np.clip(new_ventilation_rate, 0, self.max_ventilation_rate)

        # Get the step-wise PM emission, temperature and RH from the matrix based on the current ventilation rate
        if self.current_time_step < len(self.emi_pm_matrix):
            step_pm_emit = self.emi_pm_matrix[self.current_time_step]
            step_pm_dr = self.pm_dr_matrix[self.current_time_step]
            step_temp_out = self.temp_outdoor[self.current_time_step]
            step_rh_out = self.rh_outdoor[self.current_time_step]
            step_pm_out = self.pm_outdoor[self.current_time_step]
        else:
            # If we reach the end of the matrix, use the last value
            step_pm_emit = self.emi_pm_matrix[-1]
            step_pm_dr = self.pm_dr_matrix[-1]
            step_temp_out = self.temp_outdoor[-1]
            step_rh_out = self.rh_outdoor[-1]
            step_pm_out = self.pm_outdoor[-1]

        # Upadte the indoor pollutant concentration using defined fundction for mass balance
        
        # Estimate the property of outdoor air and indoor air property, and update the state
        # Outdoor: PM (0), Temp (1), RH (2)
        # Indoor: PM (3), Current Ventilation Rate (4), Indoor Temp (5), Indoor RH (6), Total Power (7)
        

        # Increment the current time step
        self.current_time_step += 1

        # Update the state
        if self.current_time_step < len(self.emi_pm_matrix):
            new_temp_out = self.temp_outdoor[self.current_time_step]
            new_rh_out = self.rh_outdoor[self.current_time_step]
            new_pm_out = self.pm_outdoor[self.current_time_step]
        else:
            # If we reach the end of the matrix, use the last value
            new_temp_out = self.temp_outdoor[-1]
            new_rh_out = self.rh_outdoor[-1]
            new_pm_out = self.pm_outdoor[-1]

        self.state = np.array([new_pm_out, new_temp_out, new_rh_out, new_pm_in, \
                             new_ventilation_rate, new_temp_in, new_rh_in, power_total],dtype=np.float32)
        
        # Define a reward function
        # Outdoor: PM (0), Temp (1), RH (2)
        # Indoor: PM (3), Current Ventilation Rate (4), Indoor Temp (5), Indoor RH (6), Total Power (7)
        # Maximize reward = -W1*EN-W2*PM_Conc-W3*dT
        reward = -(W1*self.state[7] + W2*np.maximum(0, self.state[3] - C_max) ** 2 + W3*abs(self.state[5]-T_set)/T_buffer)

        # Check if the episode is done
        done = self.current_time_step >= len(self.emi_pm_matrix)
        return self.state, reward, done, {} 
       
    def reset(self):
        # Reset the state to the initial IAQ and ventilation rate
        self.state = np.array([self.initial_pm_out, self.inital_temp_out, self.inital_rh_out, self.initial_pm_in, \
                     self.initial_ventilation_rate, self.initial_temp_in, self.initial_rh_in, self.intitial_power_total],dtype=np.float32)   
        self.current_time_step = 0
        return self.state
    
    def render(self, mode=None):
        print(f"Current IAQ: {self.state[0]}, Current Ventilation Rate: {self.state[1]}")

# Input to the Environment: emi_pm_matrix, pm_dr_matrix, temp_outdoor, rh_outdoor, pm_outdoor
# env(emi_pm_matrix,pm_dr_matrix,temp_outdoor,rh_outdoor,pm_outdoor)
env = IAQEnv(pm_emi_to_env,gamma_to_env,T_to_env,RH_to_env,pm_out_to_env)

# Define hyperparameters to tune
hyperparams = {
    'learning_rate': [0.001],
    'gamma': [0.99],
    'batch_size': [512],
    'epsilon_start': [1.0],
    'epsilon_min': [0.01],
    'epsilon_decay': [0.99],
    'target_update_freq': [200],
    'hidden_layers': [(128, 256, 128)]#, (128, 256, 128, 64)]
}

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")

class ReplayBuffer:
    def __init__(self, capacity=20000):
        self.buffer = deque(maxlen=capacity)    
    def put(self, state, action, reward, next_state, done):
        self.buffer.append([state, action, reward, next_state, done])    
    def sample(self, batch_size):
        sample = random.sample(self.buffer, batch_size)
        states, actions, rewards, next_states, done = map(np.asarray, zip(*sample))
        return np.stack(states), actions, rewards, np.stack(next_states), done    
    def size(self):
        return len(self.buffer)

class ActionStateModel(nn.Module):
    def __init__(self, state_dim, action_dim, hidden_layers, learning_rate, epsilon_start, epsilon_decay, epsilon_min):
        super(ActionStateModel, self).__init__()
        self.state_dim  = state_dim
        self.action_dim = action_dim
        self.epsilon = epsilon_start  # Starting epsilon for exploration
        
        # Create hidden layers dynamically based on input configuration
        layers = []
        input_size = state_dim
        for hidden_size in hidden_layers:
            layers.append(nn.Linear(input_size, hidden_size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.2))
            input_size = hidden_size
        layers.append(nn.Linear(input_size, action_dim))
        self.model = nn.Sequential(*layers)
        self.model.to(device)        
        self.optimizer = optim.Adam(self.parameters(), lr=learning_rate)
        self.criterion = nn.MSELoss()
    
    def forward(self, x):
        return self.model(x)
    
    def predict(self, state, num_samples=50):
        with torch.no_grad():
            state_tensor = torch.tensor(state, dtype=torch.float32).to(device)
            q_values_samples = []
            self.model.train()

            for _ in range(num_samples):
                q_values_samples.append(self.forward(state_tensor))

            q_values_samples = torch.stack(q_values_samples)
            mean_q_values = q_values_samples.mean(dim=0)
            var_q_values = q_values_samples.var(dim=0)
            return mean_q_values.cpu().numpy(), var_q_values.cpu().numpy()
    
    def get_action(self, state):
        mean_q_values, var_q_values = self.predict(state)
        epsilon = self.epsilon
        self.epsilon *= epsilon_decay
        self.epsilon = max(self.epsilon, epsilon_min)        
        if random.random() < epsilon:
            return random.randint(0, self.action_dim - 1), var_q_values.max(), mean_q_values
        else:
            return np.argmax(mean_q_values), var_q_values.max(), mean_q_values

    def train_step(self, states, targets):
        states_tensor = torch.tensor(states, dtype=torch.float32).to(device)
        targets_tensor = torch.tensor(targets, dtype=torch.float32).to(device)        
        self.optimizer.zero_grad()
        outputs = self.forward(states_tensor)
        loss = self.criterion(outputs, targets_tensor)
        loss.backward()
        self.optimizer.step()
        return loss.item()

class Agent:
    def __init__(self, env, state_dim, action_dim, learning_rate, gamma, batch_size, epsilon_start, epsilon_min, \
                 epsilon_decay, hidden_layers, target_update_freq):
        self.env = env
        self.state_dim = state_dim
        self.action_dim = action_dim
        self.model = ActionStateModel(self.state_dim, self.action_dim, hidden_layers, learning_rate, \
                                      epsilon_start, epsilon_decay, epsilon_min).to(device)
        self.target_model = ActionStateModel(self.state_dim, self.action_dim, hidden_layers, learning_rate, \
                                             epsilon_start, epsilon_decay, epsilon_min).to(device)
        self.target_update()
        self.buffer = ReplayBuffer()
        self.losses = []
        self.total_rewards = [] 
        self.uncertainties = []
        self.mean_q_values_list = []               
        self.gamma = gamma
        self.batch_size = batch_size
        self.target_update_freq = target_update_freq
        
        # CUDA Streams for asynchronous operations
        self.stream = torch.cuda.Stream()

    def target_update(self):
        self.target_model.load_state_dict(self.model.state_dict())
    
    def replay(self):
        if self.buffer.size() < self.batch_size:
            return        
        losses = []
        # for _ in range(10):
        # Create CUDA events and stream
        with torch.cuda.stream(self.stream):
                states, actions, rewards, next_states, done = self.buffer.sample(self.batch_size)
                states = torch.tensor(states, dtype=torch.float32).to(device, non_blocking=True)
                next_states = torch.tensor(next_states, dtype=torch.float32).to(device, non_blocking=True)
                rewards = torch.tensor(rewards, dtype=torch.float32).to(device, non_blocking=True)
                done = torch.tensor(done, dtype=torch.float32).to(device, non_blocking=True)
                actions = torch.tensor(actions, dtype=torch.long).to(device, non_blocking=True) 
                targets = self.target_model.forward(states)
                next_q_values = self.target_model.forward(next_states)[range(self.batch_size), \
                                torch.argmax(self.model.forward(next_states), dim=1)]
                targets[range(self.batch_size), actions] = rewards + (1 - done) * next_q_values * self.gamma                
                self.stream.synchronize()                
                loss = self.model.train_step(states.cpu().detach().numpy(), targets.cpu().detach().numpy())
                losses.append(loss)
        # self.losses.append(np.mean(losses))
    
    def train(self, max_episodes):
        global best_reward, best_target_model, best_model
        timestep_count = 0
        for ep in range(max_episodes):
            done, total_reward = False, 0
            state = self.env.reset()
            ep_mean_q_values = []
            ep_uncertainty = []
            while not done:
                action, uncertainty, mean_q_value = self.model.get_action(state)
                ep_uncertainty.append(uncertainty)
                ep_mean_q_values.append(mean_q_value)
                
                next_state, reward, done, _ = self.env.step(action)
                self.buffer.put(state, action, reward, next_state, done)
                total_reward += reward
                state = next_state
                timestep_count += 1
                if self.buffer.size() >= self.batch_size:
                    self.replay()
                if timestep_count % self.target_update_freq == 0:
                    self.target_update()
            self.uncertainties.append(np.mean(ep_uncertainty))
            self.mean_q_values_list.append(np.mean(ep_mean_q_values))
            if total_reward > best_reward:
                best_reward = total_reward
                best_target_model = copy.deepcopy(self.target_model)
                best_model = copy.deepcopy(self.model)
                torch.save(best_target_model.state_dict(), 'W1_W2_1by10_global_best_target_model_unceratinty_meanQaction.pth')
                torch.save(best_model.state_dict(), 'W1_W2_1by10__global_best_model_unceratinty_meanQaction.pth')            
            # if self.buffer.size() >= self.batch_size:
            #     self.replay()              
            self.total_rewards.append(total_reward)
            print(f'EP{ep} EpReward={total_reward} EpLosses={self.losses[-1] if self.losses else 0}')
        return self.losses, self.total_rewards, self.uncertainties, self.mean_q_values_list, best_target_model, best_model

print('Env Created, ActionState Model Defined, and Agent Activated')

best_reward = float('-inf')
best_target_model = None
best_model = None

# Define number of episodes
num_episodes = 1500

# Initialize a list to store results
results_list = []

# Generate all hyperparameter combinations
combinations = list(product(*hyperparams.values()))

for i, combination in enumerate(combinations):
    learning_rate, gamma, batch_size, epsilon_start, epsilon_min, epsilon_decay, target_update_freq, hidden_layers = combination
    print(f"Running combination {i+1}/{len(combinations)}: {combination}")
    
    # Initialize the agent with the current combination of hyperparameters
    agent = Agent(env, env.observation_space.shape[0], env.action_space.n, learning_rate, gamma, batch_size, epsilon_start, \
                  epsilon_min, epsilon_decay, hidden_layers, target_update_freq)
    
    # Train the agent
    ep_losses, ep_rewards, ep_uncertainties, ep_mean_q, _, _ = agent.train(max_episodes=num_episodes)
    
    # Append episode-wise data to results_list
    for episode in range(num_episodes):
        results_list.append({
            'learning_rate': learning_rate,
            'gamma': gamma,
            'batch_size': batch_size,
            'epsilon_start': epsilon_start,
            'epsilon_min': epsilon_min,
            'epsilon_decay': epsilon_decay,
            'target_update_freq': target_update_freq,
            'hidden_layers': hidden_layers,
            'episode': episode + 1,
            'reward': ep_rewards[episode],
            'loss': ep_losses[episode],
            'mean_q_value': ep_mean_q[episode],
            'uncertainty': ep_uncertainties[episode]
        })

    # # Convert results_list to a DataFrame
    results_df = pd.DataFrame(results_list)
     
    # # Save the results to a CSV file
    results_df.to_csv(f"W1_W2_1by10_DDQN_Tuning_Combination_Uncertainty_meanQaction_{i+1}.csv", index=False)
    print(f'Results saved')

