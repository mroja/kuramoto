#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# calculate critical K for given preset

import os
import sys
import json
import shutil
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

from common import save_plot, write_preset_to_file
from preset import get_preset, plot_preset


def calc_r_for_k(params):
    run_id = params['run_id']
    preset_name = params['preset_name']
    preset_config = params['preset_config']
    K_runs = params['K_runs']
    K_skip_init_steps = params['K_skip_init_steps']
    simulation = params['simulation']

    print('Current K: {} (run_id: {})'.format(preset_config['K'], run_id))
    
    r_sum = 0.0
    for run in range(K_runs):
        print('Run %i of %i' % (run, K_runs))

        # create new initial data
        current_preset_name = preset_name + '_' + str(run_id)
        
        preset = get_preset(preset_config)
        write_preset_to_file(current_preset_name, preset)

        r_file_name = os.path.join('dump_' + preset_name, 'crit_k', 'r_{}_{}.txt'.format(run_id, run))

        simulation.run(current_preset_name, r_file=r_file_name)
        
        try:
            os.remove(current_preset_name + '.preset')
        except Exception:
            pass

        with open(r_file_name, 'r') as f:
            lines = f.readlines()
        
        r_hist = [float(x) for x in lines[0].split()]
        r = np.mean(r_hist[K_skip_init_steps:])
        r_sum += r

        #shutil.move(r_file_name, os.path.join('dump_' + preset_name, 'k_crit', 'r_{}_{}.txt'.format(run_id, run)))

    return r_sum / float(K_runs)


def calc_crit_k(preset_name, preset_config, crit_k_config, simulation):
    K_min   = float(crit_k_config['K_min'])
    K_max   = float(crit_k_config['K_max'])
    K_steps = int(crit_k_config['K_steps'])
    K_runs  = int(crit_k_config['K_runs'])
    K_skip_init_steps = int(crit_k_config['K_skip_init_steps'])

    assert(simulation.N_steps -K_skip_init_steps > 100)

    K_values = np.linspace(K_min, K_max, num=K_steps, endpoint=True)

    # create output directory if not exists
    out_dir = os.path.join('dump_' + preset_name, 'crit_k')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, 'k_val.txt'), 'w') as f:
        for i, k_val in enumerate(K_values):
            f.write('%i %.6f\n' % (i, k_val))

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    
    args = []
    for i, k_val in enumerate(K_values):
        cfg = preset_config.copy()
        cfg['K'] = k_val
        args.append({
            'run_id': i,
            'preset_name': preset_name,
            'preset_config': cfg,
            'K_runs': K_runs,
            'K_skip_init_steps': K_skip_init_steps,
            'simulation': simulation
        })

    r_mean = pool.map(calc_r_for_k, args)

    pool.close()

    r_mean = np.array(r_mean)

    with open(os.path.join('dump_' + preset_name, 'crit_k.txt'), 'w') as f:
        for kv, rm in zip(K_values, r_mean):
            f.write('%f %f\n' % (kv, rm))

    return K_values, r_mean


def plot_crit_k(preset_name, K_values, r_mean, save=True, close=True):
    plt.figure()
    plt.xlabel('k')
    plt.ylabel('r')
    
    plt.xlim(np.min(K_values), np.max(K_values))
    plt.ylim(0, 1)
    
    plt.plot(K_values, r_mean)
    
    # TODO: plot theoretical value
    
    if save:
        save_plot(os.path.join('dump_' + preset_name, 'crit_k'))
    else:
        plt.show()

    if close:
        plt.close()
