import os
import basic_fun as bf
import biochem_constants as mbc
import fields_CHR as mfd
import init_fun
import globals4simulation as tng
import post_process_CHR as mpp
import traces_CHR as mtc

# ----


def update_senders(cell_type, activity, out_channel):
    cell_type_id = tng.strain2id_dict[cell_type]
    mask = (mfd.strain_id_map == cell_type_id)
    d_pConst = tng.dt * activity
    mfd.qIs[out_channel][mask[:]] = \
        mbc.decay_enzyme * mfd.qIs[out_channel][mask[:]] \
        + d_pConst * mbc.sender_strength

    mfd.mC[mask[:]] = \
        mfd.mC[mask[:]] + d_pConst


def update_propagators(cell_type, activity, in_channel, out_channel=None):
    cell_type_id = tng.strain2id_dict[cell_type]
    mask = (mfd.strain_id_map == cell_type_id)

    d_pConst = tng.dt * activity
    d_pQ = d_pConst * \
                  bf.aHill(mbc.inv_fold, 1,
                           mfd.AHLs[in_channel][mask[:]], mbc.h_quorum, mbc.th_quorum)

    mfd.aiiA[mask[:]] = mbc.decay_enzyme * mfd.aiiA[mask[:]] \
                        + d_pConst
    mfd.qIs[out_channel][mask[:]] = mbc.decay_enzyme * mfd.qIs[out_channel][mask[:]] \
                                        + d_pQ * mbc.propagator_strength

    mfd.mC[mask[:]] = mfd.mC[mask[:]] + d_pQ


def update_colonies(activity):
    update_senders('sender', activity, out_channel= 0)
    update_propagators('propagator', activity, in_channel=0, out_channel=0)


def update_AHLs():
    mfd.decay_AHL = 1 - tng.dt*mbc.rr_qS*(1 + mbc.r_AiiA * mfd.aiiA)
    for i in range(tng.n_quorum_channels):
        mfd.AHLs[i] = mfd.AHLs[i] * mfd.decay_AHL \
                      + tng.dt * mbc.g * mfd.qIs[i] \
                      + tng.dt * mbc.DD_qS * bf.laplacian(mfd.AHLs[i]) \

        mfd.AHLs[i][0, :] = mfd.AHLs[i][1, :]
        mfd.AHLs[i][-1, :] = mfd.AHLs[i][-2, :]
        mfd.AHLs[i][:, 0] = mfd.AHLs[i][:, 1]
        mfd.AHLs[i][:, -1] = mfd.AHLs[i][:, -2]


def simu_1condition(parameter_set, result_loc=tng.loc4condition):
    if not os.path.isdir(result_loc):
        os.makedirs(result_loc)

    flag_skip = init_fun.load_or_init(parameter_set, result_loc)

    if flag_skip:
        return

    for k in range(tng.nstep_start, tng.nstep_end):
        # update
        update_colonies(mtc.activity[k])
        update_AHLs()
        # record time trace
        if tng.flag_record:
            for i in range(tng.n_repressors):
                mtc.repressors[i,k,:] = mfd.repressors[i][mfd.all_colonie_mask]
            for i in range(tng.n_quorum_channels):
                mtc.qIs[i,k,:] = mfd.qIs[i][mfd.all_colonie_mask]
                mtc.AHLs[i, k, :] = mfd.AHLs[i][mfd.all_colonie_mask]

            mtc.mC[k] = mfd.mC[mfd.all_colonie_mask]
            mtc.AiiA[k] = mfd.aiiA[mfd.all_colonie_mask]

        # take snapshot every simulated hour
        if (k + 1) % tng.n_hr == 0:
            mpp.take_snap_shot(k, tng.loc4condition)

    mpp.save_fields(os.path.join(tng.loc4condition, tng.pickle_fname))  # save numerical results
