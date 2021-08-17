import logging
logging.basicConfig(level=logging.INFO)  # logging.INFO DEBUG

condition_name = None  # choose name of the simulation, use None to automatically generate
T = 24  # simulation duration in hour  # 1200
flag_coarse = False  # choose low / high spatial resolution setting
flag_record = False  # choose to record the time course or not
flag_snapshot = True  # choose to take snapshot every simulated hour

if flag_coarse:
    dt = 0.02  # simulation temporal resolution in hour [other value 0.02 #
    dx = 0.5  # simulation spatial resolution in mm [other value: 0.5 #
    colony_radius = 0.5  # in mm
else:  # high spatial resolution setting
    dt = 0.001  # simulation temporal resolution in hour [other value 0.02 #
    dx = 0.1  # simulation spatial resolution in mm [other value: 0.5 #
    colony_radius = 0.5 # in mm

x_distance = 2.5  # in mm, center to center distance between colonies
y_distance = 2.5  # in mm
margin_distance = 10  # in mm
nr_colonies, nc_colonies = 10, 10  # rows and columns of colonies
rank0 = (0,0)  # reference colony at upper left corner of the colony array

n_repressors = 1
n_quorum_channels = 1

# file name setting
loc4condition = None
pickle_fname = '_run_conditions.pickle'
parameter_sets = None
printed_map_design_fname = 'z_positions.txt'
result_locs = None
simulation_name = None

# constants to be derived / imported
nr, nc = None, None  # nr_colonies*2 - 1 #
strain2id_dict = None
# temporal constants
nstep_start = 0
nstep_end = None
n_hr = None
# spatial constants
c_margin = None
c_step = None
r_step = None
r0, c0 = None, None


def dprt(msg):
    from logging import info
    info(f'>>> {msg}')
    pass