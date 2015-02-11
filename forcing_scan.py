#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import shutil
import multiprocessing

from common import KuramotoSimulation, write_preset_to_file
from preset import get_preset, plot_preset
from crit_k import calc_crit_k, plot_crit_k

def run_simulation(params):
    preset_config = params['preset_config']
    forcing_frequency = params['forcing_frequency']
    forcing_strength = params['forcing_strength']
    simul = params['simul']

    preset_name = 'forcing_{}_{}'.format(forcing_frequency, forcing_strength)

    # 1) genrate new preset

    preset = get_preset(preset_config)
    write_preset_to_file(preset_name, preset)

    plot_preset(preset_name, preset, save=True)

    # 2) run simulation

    simul.forcing_strength = forcing_strength # 20.5
    simul.forcing_freq = forcing_frequency # 30.0
    simul.run(preset_name, r_file=True, mean_file=True, mean_vel_file=True)


if __name__ == '__main__':
    simul = KuramotoSimulation()
    simul.N_steps = 2000
    simul.dt = 0.01
    simul.noise = 0.0
    simul.coupling_type = 'kuramoto'
    simul.dump_interval = 10
    simul.freq_modulation_enabled = False
    simul.k_modulation_enabled = False
    simul.forcing_strength = 0.0
    simul.forcing_freq = 0.0

    preset_config = {
        "N": 500,
        "connectivity": "all_to_all",

        # used when connectivity == 'random_remove_bidirectional' or
        #           connectivity == 'random_remove_unidirectional'
        #"conn_remove_coeff": 0, 

        # used only when connectivity == 'layered_1'
        #"N_layers": 2, 
        
        "K": 0.5,

        # can be "lorentz", "normal" or "uniform"
        "freq_dist": "uniform",
        "freq_dist_scale": 10,
        "freq_dist_location": 0.0,

        "freq_dist_no_unary": False,
    }

    crit_k_config = {
        'K_min': 0.0,
        'K_max': 30.0,
        'K_steps': 20,
        'K_runs': 10,
        'K_skip_init_steps': 1000
    }

    all_connectivities = [
        "all_to_all",
        "linear_unidirectional",
        "linear_bidirectional",
        "circular_unidirectional",
        "circular_bidirectional",
        "layered_1",
        "layered_2",
        "random_remove_unidirectional",
        "random_remove_bidirectional"
    ]

    conn_remove_coeffs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # preview presets
    if 0:
        for _ in range(10):
            preset = get_preset(preset_config)
            plot_preset('preset preview', preset, save=False)
        sys.exit()

    # calculate critical K
    if 0:
        preset_name = 'forcing_scan_crit_k'
        K_values, r_mean = calc_crit_k(preset_name, preset_config, crit_k_config, simul)
        plot_crit_k(preset_name, K_values, r_mean, save=True)
        plot_crit_k(preset_name, K_values, r_mean, save=False)
        sys.exit()
    else:
        preset_config['K'] = 15.5

    forcing_presets = []
    for forcing_strength in range(0, 20, 4):
        for forcing_frequency in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
            forcing_presets.append((forcing_frequency, forcing_strength))

    # run simulation
    if 1:
        ctx = multiprocessing.get_context('spawn')
        pool = ctx.Pool(multiprocessing.cpu_count())

        args = []
        for forcing_frequency, forcing_strength in forcing_presets:
            args.append({
                'preset_config': preset_config.copy(),
                'forcing_frequency': forcing_frequency,
                'forcing_strength': forcing_strength,
                'simul': simul
            })
        #print(args)
        pool.map(run_simulation, args)
        pool.close()

    # visualize
    if 1:
        for forcing_frequency, forcing_strength in forcing_presets:
            preset_name = 'forcing_{}_{}'.format(forcing_frequency, forcing_strength)
            
            os.system('nice python3 gen_plots.py ' + preset_name)

    # copy files in one place
    if 1:
        for forcing_frequency, forcing_strength in forcing_presets:
            input_dir = 'dump_forcing_{}_{}'.format(forcing_frequency, forcing_strength)
            
            f_name = '{}_{}'.format(forcing_frequency, forcing_strength)
            
            dest_dir = 'forcing_out'
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)

            shutil.copy2(os.path.join(input_dir, 'hist.avi'), os.path.join(dest_dir, f_name + '_hist.avi'))
            shutil.copy2(os.path.join(input_dir, 'phase.avi'), os.path.join(dest_dir,  f_name + '_phase.avi'))
            shutil.copy2(os.path.join(input_dir, 'r.png'), os.path.join(dest_dir,  f_name + '_r.png'))
            shutil.copy2(os.path.join(input_dir, 'mean.png'), os.path.join(dest_dir,  f_name + '_mean.png'))
            shutil.copy2(os.path.join(input_dir, 'mean_vel.png'), os.path.join(dest_dir,  f_name + '_mean_vel.png'))
