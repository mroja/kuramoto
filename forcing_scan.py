#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json

from common import KuramotoSimulation, write_preset_to_file
from preset import get_preset, plot_preset
from crit_k import calc_crit_k, plot_crit_k


simul = KuramotoSimulation()
simul.N_steps = 1500
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
    'K_max': 3.0,
    'K_steps': 20,
    'K_runs': 10,
    'K_skip_init_steps': 2000
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

# calculate critical K
if 1:
    K_values, r_mean = calc_crit_k('forcing_scan_crit_k', preset_config, crit_k_config, simul)
    plot_crit_k(preset_name, K_values, r_mean, save=True)
    plot_crit_k(preset_name, K_values, r_mean, save=False)
    sys.exit()


conn_remove_coeffs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

name_i = 10
for fs in [10]: # xrange(0, 20, 4):
    for ff in [10]:#, 0.5, 1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
        
        preset_name = 'osc_' + str(name_i)
        
        # 1) genrate new preset

        preset = get_preset(preset_config)
        write_preset_to_file(preset_name, preset)

        plot_preset(preset_name, preset, save=True)

        # 2) run simulation

        simul.forcing_strength = fs # 20.5
        simul.forcing_freq = ff # 30.0
        simul.run(preset_name, r_only=True, r_file='r.txt')

        # 3) generate plots from dumps
        
        #os.system('python gen_plots.py ' + name)
