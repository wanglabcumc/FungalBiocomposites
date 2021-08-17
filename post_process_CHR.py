import numpy as np
import matplotlib.pyplot as plt
plt.interactive(False)
import traces_CHR as mtc

import basic_fun as bf
import fields_CHR as mfd
import globals4simulation as tng
import pickle
import os

# ----

def take_snap_shot(k, stack_loc):
    # k is time_step_index
    base_name = str(int((k + 1) / tng.n_hr)).zfill(4)
    png_path = os.path.join(stack_loc, 'AHL_' + base_name + '.png')

    qS = mfd.AHLs.sum(axis=0)
     # save raw nparray, takes more space, get the option to select range
    np_path = os.path.join(stack_loc, f'raw_AHL_{base_name}.npy')
    np.save(np_path, qS)


def gen_summary_fig(save_impath=None):
    # generate summary after all conditions are simulated
    print('generating report')
    fig, axs = plt.subplots(2,3)
    for i, result_loc in enumerate(tng.result_locs):
        # get the pickle, and read AHL field
        pickle_path = os.path.join(result_loc, tng.pickle_fname)
        with open(pickle_path, 'rb') as f:
            # graph 1
            fields, traces, old_nstep_end = pickle.load(f)
            AHL_plot = fields[4]
            bf.show_field(AHL_plot.sum(axis=0), fig=fig,
                          ax=axs[i, 0], log = True)  # , vmax=1, vmin=0

            # graph 2
            if tng.flag_record:
                repressors, qIs, AiiA, qSs, mC = traces
                t = mtc.t
                bf.show_plot(t, qSs[0],ax = axs[i,1], log = True)

            # graph 3
            all_AHL = AHL_plot.sum(axis=0)
            ahl_line = all_AHL[tng.r0,:]
            x = range(len(ahl_line))
            bf.show_plot(x, ahl_line, ax = axs[i,2], vmin = 0, log = False)
            # bf.show_plot(t, qIs[0],ax = axs[i,2], log = True)

    if save_impath != None:
        plt.savefig(save_impath, dpi=300)
        os.startfile(save_impath)


def save_fields(pickle_path):
    # save simulation raw data
    fields = [mfd.strain_id_map, mfd.repressors, mfd.qIs, mfd.aiiA, mfd.AHLs, mfd.mC]
    traces = [mtc.repressors, mtc.qIs, mtc.AiiA, mtc.AHLs, mtc.mC]
    number = tng.nstep_end
    with open(pickle_path, 'wb') as f:
        pickle.dump([fields, traces, number], f)