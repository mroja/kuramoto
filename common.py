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
        self.freq_modulation_ampl = 0.0
        self.freq_modulation_freq = 0.0
        self.freq_modulation_offset = 0.0
        self.k_modulation_enabled = False
        self.k_modulation_ampl = 0.0
        self.k_modulation_freq = 0.0
        self.k_modulation_offset = 0.0

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
            '{freq_modulation_ampl:f}\n' + \
            '{freq_modulation_freq:f}\n' + \
            '{freq_modulation_offset:f}\n' + \
            '{k_modulation_enabled}\n' + \
            '{k_modulation_ampl:f}\n' + \
            '{k_modulation_freq:f}\n' + \
            '{k_modulation_offset:f}'
        args = {
           'N_steps':                 self.N_steps,
           'dt':                      self.dt,
           'noise':                   self.noise,
           'dump_interval':           self.dump_interval,
           'coupling_type':           self.coupling_type,
           'forcing_strength':        self.forcing_strength,
           'forcing_freq':            self.forcing_freq,
           'freq_modulation_enabled': 'freq_modulation' if self.freq_modulation_enabled else 'no_freq_modulation',
           'freq_modulation_ampl':    self.freq_modulation_ampl,
           'freq_modulation_freq':    self.freq_modulation_freq,
           'freq_modulation_offset':  self.freq_modulation_offset,
           'k_modulation_enabled':    'k_modulation' if self.k_modulation_enabled else 'no_k_modulation',
           'k_modulation_ampl':       self.k_modulation_ampl,
           'k_modulation_freq':       self.k_modulation_freq,
           'k_modulation_offset':     self.k_modulation_offset
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
            self.freq_modulation_ampl = float(lines[8])
            self.freq_modulation_freq = float(lines[9])
            self.freq_modulation_offset = float(lines[10])
            self.k_modulation_enabled = True if lines[11] == 'k_modulation' else False
            self.k_modulation_ampl = float(lines[12])
            self.k_modulation_freq = float(lines[13])
            self.k_modulation_offset = float(lines[14])

def load_preset_from_file(name):
    with open(name + '.bin', 'rb') as f:
        N = struct.unpack('I', f.read(4))[0]
        freq = np.fromfile(f, dtype=np.float64, count=N)
        phase = np.fromfile(f, dtype=np.float64, count=N)
        k = np.fromfile(f, dtype=np.float64, count=N*N)
    return DataPreset(N, k, freq, phase)

def write_preset_to_file(file_name, preset):
    with open(file_name, 'wb') as f:
        f.write(struct.pack('I', preset.N))
        preset.freq.astype(np.float64).tofile(f)
        preset.phase.astype(np.float64).tofile(f)
        preset.k.astype(np.float64).tofile(f)

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
