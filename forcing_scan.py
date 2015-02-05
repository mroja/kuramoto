#!/usr/bin/env python

import os
import sys
from common import SimulationConfig

config = SimulationConfig()
config.N_steps = 1500
config.dt = 0.01
config.noise = 0.0
config.coupling_type = 'kuramoto'
config.dump_interval = 10

config.freq_modulation_enabled = True
config.freq_modulation_ampl = 1.0
config.freq_modulation_freq = 1.0
config.freq_modulation_offset = 0.0
config.k_modulation_enabled = True
config.k_modulation_ampl = 1.0
config.k_modulation_freq = 1.0
config.k_modulation_offset = 0.0

get_preset_config = {
    "N": 200,
    "connectivity": "all_to_all",
    "K": 0.5,
    #"freq_dist": "lorentz",
    #"freq_dist": "normal",
    "freq_dist": "uniform",
    "freq_dist_no_unary": false,
    "freq_dist_scale": 10,
    "freq_dist_location": 0.0
}

connectivities = [
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
for fs in xrange(0, 20, 4):
    for ff in [0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
        name = 'osc_' + str(name_i)
        config.forcing_strength = fs # 20.5
        config.forcing_freq = ff #30.0
        config.write_to_file(name)
        os.system('python gen_preset.py ' + name)
        os.system('kuramoto_simulation.exe ' + name + ' --only-r -r r.txt')
        os.system('python gen_plots.py ' + name)
