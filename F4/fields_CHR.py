# global spatial variables

# constants
coordinate_map = None

# to be initialized before running simulation
strain_id_map = None
basic_strain_id_map = None
AHLs = None  # [AHL], may be several orthogonal ones
qIs = None  # [quorum producer], may be several orthogonal ones
mC = None  # [mCherry]
aiiA = None
decay_AHL = None  # AHL decary rate at each grid point
repressors = None  # [repressor], may be several orthogonal ones


all_colonie_mask = None
n_colonies = None
r_probe, c_probe = None, None

colony_positions = None