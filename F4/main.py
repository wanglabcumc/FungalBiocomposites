#! python3.8
import time
import os
import basic_fun as mbf
import simulation_fun as msf
import init_fun as mif
import globals4simulation as tng
import post_process_CHR as mpp

# ----


def main():
    mif.set_derived_const()
    # mif.gen_colony_map_from_design(tng.printed_map_design_fname)
    mif.gen_colony_map()

    tng.parameter_sets = ['0', '1']    #
    simu_conditions(tng.parameter_sets, loc4condition_set= tng.condition_name)  # '200806_1953'


def simu_conditions(parameter_sets, loc4condition_set=None):
    # simulate a set of conditions
    tng.result_locs = []
    if loc4condition_set == None:
        loc4condition_set = mbf.datetime_stamp()
        if not os.path.isdir(loc4condition_set):
            os.makedirs(loc4condition_set)

    tng.simulation_name = loc4condition_set

    for parameter_set in parameter_sets:  # [0]:  # [1] #len(parameter_sets)
        condition_str = parameter_set
        tng.loc4condition = os.path.join(loc4condition_set, condition_str)

        print(f'simulating: {condition_str}')
        start_time = time.time()
        msf.simu_1condition(parameter_set, result_loc=tng.loc4condition)  # old_dir_name)  #

        print("  run time: {:.2f} s".format((time.time() - start_time)))
        print(f'    results saved in {tng.loc4condition}')
        tng.result_locs.append((tng.loc4condition))

    mpp.gen_summary_fig(f'{tng.simulation_name}.png')


if __name__ == '__main__':
    main()