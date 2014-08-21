#!/usr/bin/env python
# -*- coding: utf-8 -*-

# calculate critical K for given preset

import os
import sys
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt

from common import save_plot

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "Usage: calc_crit_k.py preset_name [--no-plot]"
        sys.exit()

    preset_name = sys.argv[1]
    with open(preset_name + '.preset.json', 'r') as f:
        config = json.load(f)

    K_min   = float(config['K_min'])
    K_max   = float(config['K_max'])
    K_steps = int(config['K_steps'])
    K_runs  = int(config['K_runs'])
    K_skip_init_steps = int(config['K_skip_init_steps'])

    K_values = np.linspace(K_min, K_max, num=K_steps, endpoint=True)

    out_dir = os.path.join('dump_' + preset_name, 'k_crit')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, 'k_val.txt'), 'w') as f:
        for i, k_val in enumerate(K_values):
            f.write('%i %.6f\n' % (i, k_val))

    r_mean = np.zeros(len(K_values))
    for i, k_val in enumerate(K_values):
        print 'Current K: %f (index %i)' % (k_val, i)
        r_sum = 0
        for run in range(K_runs):
            print 'Run %i of %i' % (run, K_runs)

            # create new initial data
            os.system('python gen_preset.py %s %.6f' % (preset_name, k_val))
            os.system('kuramoto_simulation %s --only-r' % (preset_name,))

            f_name = os.path.join('dump_' + preset_name, 'steps', 'r.txt')
            with open(f_name, 'r') as f:
                lines = f.readlines()
            r_hist = [float(x) for x in lines[0].split()]
            r = np.mean(r_hist[K_skip_init_steps:])
            r_sum += r

            shutil.move(f_name, os.path.join(out_dir, str(i) + '.' + str(run) + '.txt'))

        r_mean[i] = r_sum / float(K_runs)


    if len(sys.argv) >= 3 and sys.argv[2] != '--no-plot':
        plt.xlabel('k')
        plt.ylabel('r')
        plt.xlim(K_min, K_max)
        plt.ylim(0, 1)
        plt.plot(K_values, r_mean)
        # TODO: plot theoretical value
        save_plot(os.path.join('dump_' + preset_name, 'crit_k'))

    with open(os.path.join('dump_' + preset_name, 'crit_k.txt'), 'w') as f:
        for z in zip(K_values, r_mean):
            f.write('%f %f\n' % (z[0], z[1]))
