#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

# N     - number of oscillators
# freq  - oscillator frequency (N elements)
# phase - oscillator phases (N elements)
# k     - coupling coefficients (NxN matrix)
DataPreset = namedtuple('DataPreset', ['N', 'k', 'freq', 'phase'])

class SimulationConfig(object):
    def __init__(self):
        self.N_steps = 0
        self.dt = 0.01
        self.noise = 0.0
        self.dump_interval = 0
        self.coupling_type = 'kuramoto'
        self.forcing_strength = 0.0
        self.forcing_freq = 0.0
        self.freq_modulation_enabled = False
        self.k_modulation_enabled = False

    def write_to_file(self, preset_name):
        conf_str = \
            '{N_steps:d}\n' + \
            '{dt:f}\n' +  \
            '{noise:f}\n' + \
            '{dump_interval:d}\n' + \
            '{coupling_type}\n' + \
            '{forcing_strength:f}\n' + \
            '{forcing_freq:f}\n' + \
            '{freq_modulation_enabled}\n' + \
            '{k_modulation_enabled}'
        args = {
           'N_steps':                 self.N_steps,
           'dt':                      self.dt,
           'noise':                   self.noise,
           'dump_interval':           self.dump_interval,
           'coupling_type':           self.coupling_type,
           'forcing_strength':        self.forcing_strength,
           'forcing_freq':            self.forcing_freq,
           'freq_modulation_enabled': 'freq_modulation' if self.freq_modulation_enabled else 'no_freq_modulation',
           'k_modulation_enabled':    'k_modulation' if self.k_modulation_enabled else 'no_k_modulation',
        }
        conf = conf_str.format(**args)
        with open(preset_name + '.conf.txt', 'w') as f:
            f.write(conf)

    def read_from_file(self, preset_name):
        with open(preset_name + '.conf.txt', 'r') as f:
            lines = f.readlines()
            self.N_steps = int(lines[0])
            delf.dt = float(lines[1])
            self.noise = float(lines[2])
            self.dump_interval = int(lines[3])
            self.coupling_type = lines[4]
            self.forcing_strength = float(lines[5])
            self.forcing_freq = float(lines[6])
            self.freq_modulation_enabled = True if lines[7] == 'freq_modulation' else False
            self.k_modulation_enabled = True if lines[11] == 'k_modulation' else False

def load_preset_from_file(name):
    with open(name + '.preset', 'rb') as f:
        N = struct.unpack('I', f.read(4))[0]
        freq = np.fromfile(f, dtype=np.float64, count=N)
        phase = np.fromfile(f, dtype=np.float64, count=N)
        k = np.fromfile(f, dtype=np.float64, count=N*N)
    return DataPreset(N, k, freq, phase)

def write_preset_to_file(file_name, preset):
    with open(file_name + '.preset', 'wb') as f:
        f.write(struct.pack('I', preset.N))
        preset.freq.astype(np.float64).tofile(f)
        preset.phase.astype(np.float64).tofile(f)
        preset.k.astype(np.float64).tofile(f)

def write_freq_modul_data_to_file(file_name, freq_ampl, freq_freq, freq_offset):
    with open(file_name + '.fm.preset', 'wb') as f:
        freq_ampl.astype(np.float64).tofile(f)
        freq_freq.astype(np.float64).tofile(f)
        freq_offset.astype(np.float64).tofile(f)

def write_k_modul_data_to_file(file_name, k_ampl, k_freq, k_offset):
    with open(file_name + '.km.preset', 'wb') as f:
        k_ampl.astype(np.float64).tofile(f)
        k_freq.astype(np.float64).tofile(f)
        k_offset.astype(np.float64).tofile(f)

def save_plot(path, ext='png', close=True):
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
    if not os.path.exists(directory):
        os.makedirs(directory)
    savepath = os.path.join(directory, filename)
    # print("Saving figure to '%s'" % savepath)
    plt.savefig(savepath)
    if close:
        plt.close()

