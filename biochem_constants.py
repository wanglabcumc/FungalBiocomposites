import math
import numpy as np
import globals4simulation as tng

th_quorum = 0.035  # threshold for quorum sensors
r_AiiA=10  # the decay rate boost over agar background by 1 unit of AiiA
propagator_strength = 1  # for tuning the expression of cinI in individual strains
sender_strength = 1

def colony_activity(t):  # define metabolic activity of a colony
    activity = 1 / (1 + np.exp((8 - t))) / (1 + np.exp((t - 34) / 4))
    return activity

tau_cell = 0.5 / math.log(2)  # double time 0.5 hr
tau_enzyme = 8 / math.log(2)  # enzyme degradation/inactivation time 4 hr
tau_TF = 0.5 / math.log(2)  # repressor dilution time 0.5 hr

h_mc = 2  # hill coefficient
h_quorum = 2
h_repressor = 2

th_repressor = 0.05  # threshold repressor concentration to repress a promoter

inv_fold = 1/100  # fold change of promoters

DD_qS = 2 / tng.dx ** 2  # diffusion coefficient in terms of grid
rr_qS = 0.2  # AHL degradation rate hr-1

g = 1  # for globally tuning production strength

# intermediates to speed up simulation
decay_TF = (1 - tng.dt / tau_TF)
decay_cell = (1 - tng.dt / tau_cell)
decay_enzyme = (1 - tng.dt / tau_enzyme)
