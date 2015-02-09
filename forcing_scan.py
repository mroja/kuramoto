#!/usr/bin/env python

import os
import sys
import json
from common import SimulationConfig

config = SimulationConfig()
config.N_steps = 1500
config.dt = 0.01
config.noise = 0.0
config.coupling_type = 'kuramoto'
config.dump_interval = 10

config.freq_modulation_enabled = False
config.freq_modulation_ampl = 0.0
config.freq_modulation_freq = 1.0
config.freq_modulation_offset = 0.0
config.k_modulation_enabled = False
config.k_modulation_ampl = 0.0
config.k_modulation_freq = 1.0
config.k_modulation_offset = 0.0

preset_config = {
    "N": 200,
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

name_i = 10
for fs in [10]: # xrange(0, 20, 4):
    for ff in [10]:#, 0.5, 1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
        
        name = 'osc_' + str(name_i)
        
        # 1) genrate preset file

        with open(name + '.preset.json', 'w') as f:
            f.write(json.dumps(preset_config))
        
        os.system('python gen_preset.py ' + name)

        # 2) run simulation

        config.forcing_strength = fs # 20.5
        config.forcing_freq = ff # 30.0
        config.write_to_file(name)
        
        cmd = 'kuramoto_simulation ' + name + ' --only-r -r r.txt'
        print cmd
        os.system(cmd)

        # 3) run simulation
        
        #os.system('python gen_plots.py ' + name)
