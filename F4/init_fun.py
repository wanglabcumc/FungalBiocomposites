import numpy as np

import globals4simulation as tng
import basic_fun as bf
import fields_CHR as mfd
import biochem_constants as mbc
import traces_CHR as mtc
import os
import pickle

# ----
def set_temporal_points():
    mtc.t = np.linspace(0, tng.dt * tng.nstep_end, tng.nstep_end, endpoint=False)
    mtc.activity = mbc.colony_activity(mtc.t)


def wipe_fields():
    mfd.AHLs = np.zeros((tng.n_quorum_channels, tng.nr, tng.nc))  # np.zeros((6, nr, nc))
    mfd.qIs = np.zeros((tng.n_quorum_channels, tng.nr, tng.nc))
    mfd.mC = np.zeros((tng.nr, tng.nc))
    mfd.aiiA = np.ones((tng.nr, tng.nc))  # np.zeros((nr, nc))  #
    mfd.repressors = 1e-3 * np.ones((tng.n_repressors, tng.nr, tng.nc))


def place_colony(canvas, colony_rank, strain_id):
    colony_position = bf.rank2pos(colony_rank)
    bf.gen_circle(canvas, colony_position, radius=tng.colony_radius/tng.dx, set_value=strain_id)


def mod_colony_map(parameter_set):
    mfd.strain_id_map = mfd.basic_strain_id_map.copy()
    print(parameter_set)
    for i in range(1):
        if parameter_set[i] == '0':
            place_colony(mfd.strain_id_map, tng.rank0, -1)
        else:
            place_colony(mfd.strain_id_map, tng.rank0, tng.strain2id_dict['sender'])


def gen_colony_map():
    tng.nr = 2 * tng.c_margin + 1 + (tng.nr_colonies - 1) * tng.r_step
    tng.nc = 2 * tng.c_margin + 1 + (tng.nc_colonies - 1) * tng.c_step  # +   # 401, 401  #  201, 251  #
    tng.r0, tng.c0 = tng.c_margin, tng.c_margin

    tng.strain2id_dict = {'sender': 1, 'propagator': 2}
    mfd.basic_strain_id_map = -1 * np.ones((tng.nr, tng.nc), dtype=np.int8)
    for n in range(tng.nr_colonies):
        for c in range(tng.nc_colonies):
            place_colony(mfd.basic_strain_id_map, (n,c), strain_id=2)

    mfd.coordinate_map = np.zeros((tng.nr, tng.nc, 2), dtype=np.uint16)
    for i in range(tng.nr):
        for j in range(tng.nc):
            mfd.coordinate_map[i, j] = i, j


def gen_colony_map_from_design(design_filename):
    # placing colonies
    lines = []
    with open(design_filename, 'r') as f:
        for line in f:
            lines.append(line)

    coordinates = []
    cell_types = []
    tng.strain2id_dict = {}
    current_id = 0
    for i in range(len(lines)):  # read coordinate and cell type assignment
        words = lines[i].strip().split(', ')
        # print(words)
        num_string = words[0]
        for letter in num_string.strip():
            j = int(letter)
            r = i + j
            c = j - i
            coordinates.append(bf.rank2pos((r, c)))

        for word in words[1:]:
            cell_types.append(word)
            if word not in tng.strain2id_dict:
                tng.strain2id_dict[word] = current_id
                current_id += 1

    coordinates = np.array(coordinates, dtype=np.int16)
    coordinates[:, 0] = coordinates[:, 0] - coordinates[:, 0].min()
    coordinates[:, 1] = coordinates[:, 1] - coordinates[:, 1].min()

    tng.nr_colonies = coordinates[:, 0].max() + 1
    tng.nc_colonies = coordinates[:, 1].max() + 1

    # r0, c0 = c_margin, c_margin  # int(c_margin + c_step * (n_row-1)/2)

    tng.nr = 2 * tng.c_margin + 1 + (tng.nr_colonies - 1) * tng.r_step
    tng.nc = 2 * tng.c_margin + 1 + (tng.nc_colonies - 1) * tng.c_step  # +   # 401, 401  #  201, 251  #

    # placing the colonies according to design
    mfd.basic_strain_id_map = -1 * np.ones((tng.nr, tng.nc), dtype=np.int8)
    for rank, strain_name in zip(coordinates[:, ], cell_types):
        r, c = rank
        mfd.basic_strain_id_map[bf.rank2pos((r, c))] = tng.strain2id_dict[strain_name]

    mfd.coordinate_map = np.zeros((tng.nr, tng.nc, 2), dtype=np.uint16)
    for i in range(tng.nr):
        for j in range(tng.nc):
            mfd.coordinate_map[i, j] = i, j


def update_colony_info():
    mfd.all_colonie_mask = (mfd.strain_id_map >= 0)
    mfd.n_colonies = mfd.all_colonie_mask.sum()
    mfd.colony_positions = mfd.coordinate_map[mfd.all_colonie_mask[:]]


def wipe_traces():
    set_temporal_points()
    if tng.flag_record:
        mtc.repressors = np.zeros((tng.n_repressors, tng.nstep_end, mfd.n_colonies))
        mtc.qIs = np.zeros((tng.n_quorum_channels, tng.nstep_end, mfd.n_colonies))
        mtc.AHLs = np.zeros((tng.n_quorum_channels, tng.nstep_end, mfd.n_colonies))
        mtc.mC = np.zeros((tng.nstep_end, mfd.n_colonies))
        mtc.AiiA = np.zeros((tng.nstep_end, mfd.n_colonies))

def set_derived_const():
    tng.c_margin = int(tng.margin_distance / tng.dx)
    tng.c_step = int(tng.x_distance / tng.dx)  # 2.5
    tng.r_step = int(tng.y_distance/ tng.dx)
    tng.r0, tng.c0 = tng.c_margin, tng.c_margin  # int(c_margin + c_step * (n_row-1)/2)
    tng.n_hr = int(1 / tng.dt)
    tng.nstep_end = int(tng.T / tng.dt)

    set_temporal_points()


def load_or_init(parameter_set, data_loc = None):
    # if pickle exists: load into init condition
    # else create a new condition with parameter_set
    pickle_path = os.path.join(tng.loc4condition, tng.pickle_fname)

    if not os.path.exists(pickle_path):  # start anew
        wipe_fields()  # place colonies, wipe to 0 or 1
        mod_colony_map(parameter_set)
        update_colony_info()
        wipe_traces()  # assign memory space to store traces
        tng.nstep_start = 0

    else:  # load pickle
        print('found pickle of previous simulation')
        with open(pickle_path, 'rb') as f:
            fields, traces, old_nstep_end = pickle.load(f)

        tng.nstep_start = old_nstep_end
        if tng.nstep_start >= tng.nstep_end:
            return 1

        mfd.strain_id_map, mfd.repressors, mfd.qIs, mfd.aiiA, mfd.AHLs, mfd.mC = fields
        update_colony_info()

        wipe_traces()  # assign memory space to store traces
        old_repressors, old_qIs, old_aiiA, old_qSs, old_mC = traces
        if tng.flag_record:
            mtc.repressors[:, 0:old_nstep_end, :] = old_repressors
            mtc.qIs[:, 0:old_nstep_end, :] = old_qIs
            mtc.AiiA[0:old_nstep_end, :] = old_aiiA
            mtc.AHLs[:, 0:old_nstep_end, :] = old_qSs
            mtc.mC[0:old_nstep_end, :] = old_mC

        return 0




