#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO:
# - scaling of freq distributions

import os
import sys
import math
import json
import struct
import numpy as np
from common import DataPreset, write_preset_to_file

# Available presets:
# * simple_freq
# * simple_phase
# * simple_cluster
# * simple_anti_phase
# * [custom]

def get_simple_preset(preset, g_dtype=np.float64):
    if preset == 'simple_freq': # setup for freq sync
        N = 8
        k = np.zeros((N,N), dtype=g_dtype)
        omega = np.ones(N, dtype=g_dtype)
        phase = np.zeros(N, dtype=g_dtype)
        i = 0
        while i < (N-1):
            if i%2:
                a = 4.0
            else:
                a = 2.0
            k[i+1,i] = k[i,i+1] = a
            i += 1
        #k[N-1,0] = k[0,N-1] = 4.0 # driving oscillator ???
        k *= 1.0/np.sqrt(3.0)
        omega[0] = omega[2] = omega[4] = omega[6] = 3.0 * np.sqrt(3.0)
        omega[1] = omega[3] = omega[5] = omega[7] = 1.0 * np.sqrt(3.0)
    elif preset == 'simple_phase': # setup for phase sync
        N = 8
        k = np.zeros((N,N), dtype=g_dtype)
        omega = np.ones(N, dtype=g_dtype)
        phase = np.zeros(N, dtype=g_dtype)
        i = 0
        while i < (N-1):
            if i%2:
                a = 2
            else:
                a = 1
            k[i+1,i] = k[i,i+1] = a
            i += 1
        k[N-1,0] = k[0,N-1] = 2
        # all omegas == 1
        phase = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6], dtype=g_dtype)
    elif preset == 'simple_cluster': # setup for cluster sync
        N = 8
        k = np.zeros((N,N), dtype=g_dtype)
        omega = np.ones(N, dtype=g_dtype)
        phase = np.zeros(N, dtype=g_dtype)
        e1 = 1.0/32.0
        e2 = 1.0/34.0
        e3 = 1.0/35.0
        e4 = np.sqrt(3.0)
        k = np.array([[0., 3., 0., e1, e1, e1, e3, e3],
                      [3., 0., 2., e1, e1, e1, e3, e3],
                      [0., 2., 0., e1, e1, e1, e3, e3],
                      [e1, e1, e1, 0., e4, 0., e2, e2],
                      [e1, e1, e1, e4, 0., 1., e2, e2],
                      [e1, e1, e1, 0., 1., 0., e2, e2],
                      [e3, e3, e3, e2, e2, e2, 0., 1.],
                      [e3, e3, e3, e2, e2, e2, 1., 0.]], dtype=g_dtype)
        omega[0] = omega[1] = omega[2] = 3.0/32.0 + (3.0*np.sqrt(3.0))/35.0
        omega[3] = omega[4] = omega[5] = 1.0/34.0 + (2.0*np.sqrt(3.0))/35.0
        omega[6] = omega[7] = np.sqrt(3.0)/70.0
        phase = np.array([0.5, 1.0, 3.0, 4.0, 2.0, 1.5, 3.0, 3.5], dtype=g_dtype)
    elif preset == 'simple_anti_phase': # setup for anti-phase sync
        N = 7
        k = np.zeros((N,N), dtype=g_dtype)
        omega = np.ones(N, dtype=g_dtype)
        phase = np.zeros(N, dtype=g_dtype)
        e1 = np.sqrt(3.0)
        k = np.array([[ 0.,     e1-1.0, 0., -2.,    0.,    0.,    0.],
                      [ e1-1.0, 0.,     0., -2.,    0.,    0.,    0.],
                      [ 0.,     0.,     0., -5.,    2.,    0.,   -3.],
                      [-2.,    -2.,    -5.,  0.,    0.,    1./3., 0.],
                      [ 0.,     0.,     2.,  0.,    0.,    2./3., 6.],
                      [ 0.,     0.,     0.,  1./3., 2./3., 0.,    0.],
                      [ 0.,     0.,    -3.,  0.,    6.,    0.,    0.]], dtype=g_dtype)
        # all omegas == 1
        phase = np.array([1.0, 2.0, 3.0, 2.53, 3.5, 2.0, 4.0], dtype=g_dtype)
    else:
        raise Exception('get_simple_preset: unknown preset name: ' + str(preset))
    return DataPreset(N, k, omega, phase)

def get_preset(config, g_dtype=np.float64):
    if 'N' not in config:
        raise Exception('get_preset: number of oscillators not specified')

    N = int(config['N'])

    k = np.zeros((N,N), dtype=g_dtype)

    connectivity = config.setdefault('connectivity', 'all_to_all')
    if connectivity == 'all_to_all':
        k += 1.0
    elif connectivity == 'linear_unidirectional':
        for i in xrange(0, N-1):
            k[i, i+1] = 1.0
    elif connectivity == 'linear_bidirectional':
        for i in xrange(0, N-1):
            k[i, i+1] = 1.0
            k[i+1, i] = 1.0
    elif connectivity == 'circular_unidirectional':
        for i in xrange(0, N-1):
            k[i, i+1] = 1.0
        k[N-1, 0] = 1.0
    elif connectivity == 'circular_bidirectional':
        for i in xrange(0, N-1):
            k[i, i+1] = 1.0
            k[i+1, i] = 1.0
        k[N-1, 0] = 1.0
        k[0, N-1] = 1.0
    elif connectivity == 'random_remove_unidirectional':
        k += 1.0
        N_conn = N*N - N
        N_remove = int(round(N_conn*float(config.setdefault('conn_remove_coeff', 0.5))))
        print 'Removing %i connections from %i' % (N_remove, N_conn)
        removed_pairs = set()
        for i in xrange(N_remove):
            while True:
                x = np.random.randint(low=0, high=N, size=1)[0]
                y = np.random.randint(low=0, high=N, size=1)[0]
                if x == y:
                    continue
                if (x, y) not in removed_pairs:
                    removed_pairs.add((x,y))
                    break
            k[x, y] = 0.0
    elif connectivity == 'random_remove_bidirectional':
        k += 1.0
        N_conn = (N*N - N) / 2
        N_remove = int(round(N_conn*float(config.setdefault('conn_remove_coeff', 0.5))))
        print 'Removing %i connections from %i' % (N_remove, N_conn)
        removed_pairs = set()
        for i in xrange(N_remove):
            while True:
                x = np.random.randint(low=0, high=N, size=1)[0]
                y = np.random.randint(low=0, high=N, size=1)[0]
                if x == y:
                    continue
                if (x, y) not in removed_pairs and (y, x) not in removed_pairs:
                    removed_pairs.add((x,y))
                    removed_pairs.add((y,x))
                    break
            k[x, y] = 0.0
            k[y, x] = 0.0
    elif connectivity == 'layered_1':
        N_layers = int(config.setdefault('N_layers', 2))
        N_osc_per_layer = int(N / N_layers)
        if N_osc_per_layer * N_layers != N:
            raise Exception('N must be divisable by N_layers')
        # create interconnected clusters
        offset = 0
        for i in xrange(N_layers):
            k[offset:offset + N_osc_per_layer, offset:offset + N_osc_per_layer] = 1.0
            offset += N_osc_per_layer
        # connect clusters
        for layer_index in xrange(N_layers-1):
            for osc_index in xrange(N_osc_per_layer):
                x = layer_index * N_osc_per_layer + osc_index
                y = (layer_index+1) * N_osc_per_layer + osc_index
                k[x, y] = 1.0
                k[y, x] = 1.0
        #print k
    elif connectivity == 'layered_2':
        N_layers = 2
        N_osc_per_layer = int(N / N_layers)
        if N_osc_per_layer * N_layers != N:
            raise Exception('N must be divisable by N_layers')
        for layer_index in xrange(N_layers-1):
            for layer_osc_index in xrange(N_osc_per_layer):
                for osc_index in xrange(N_osc_per_layer):
                    x = layer_index * N_osc_per_layer + layer_osc_index
                    y = (layer_index+1) * N_osc_per_layer + osc_index
                    k[x, y] = 1.0
                    k[y, x] = 1.0
        print k

    for i in xrange(N): # zeros on the diagonal
        k[i, i] = 0.0

    k *= float(config.setdefault('K', 1.0)) # apply final K value
    k /= float(N) # divide by N here, simulation doesn't do this

    freq_dist = config.setdefault('freq_dist', 'lorentz')
    freq_dist_scale = config.setdefault('freq_dist_scale', 1.0)
    freq_dist_location = config.setdefault('freq_dist_location', 0.0)

    if freq_dist == 'lorentz':
        # Cauchy–Lorentz distribution as used by Kuramoto and Daido
        omega = np.random.standard_cauchy(size=N).astype(g_dtype)
        omega *= freq_dist_scale * 0.9
        omega += freq_dist_location
    elif freq_dist == 'normal':
        omega = np.random.normal(loc=freq_dist_location, scale=freq_dist_scale, size=N).astype(g_dtype)
    elif freq_dist == 'uniform':
        omega = np.random.uniform(low=-1.1*freq_dist_scale,
                                  high=1.1*freq_dist_scale,
                                  size=N).astype(g_dtype)
        omega += freq_dist_location

    # NOTE: apply this only to symmetric distributions
    if bool(config.setdefault('freq_dist_no_unary', False)):
        if freq_dist_location == 0.0:
            print 'WARNING: probably applying freq_dist_no_unary to non-symmetric distribution'
        for i in xrange(N):
            if omega[i] < 0.0:
                omega[i] = -omega[i]

    phase = np.random.uniform(low=0.0, high=2.0*math.pi, size=N).astype(g_dtype)

    return DataPreset(N, k, omega, phase)

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "Usage: gen_preset.py preset_name|config_file [K]"
        print "If K is supplied it will overwrite one specified in config_file."
        sys.exit()

    preset_name = sys.argv[1]
    if preset_name.startswith('simple_'):
        preset = get_simple_preset(preset_name)
    else:
        with open(preset_name + '.preset.json', 'r') as f:
            config = json.load(f)
        if len(sys.argv) >= 3:
            config['K'] = float(sys.argv[2])
        preset = get_preset(config)

    preset_out_file_name = preset_name + '.preset'
    try:
        os.remove(preset_out_file_name)
    except Exception:
        pass

    #print preset.k
    #print preset.freq
    #print preset.phase
    write_preset_to_file(preset_name, preset)
    
    if len(sys.argv) < 3:
        import matplotlib.pyplot as plt
        from common import save_plot    
    
        plt.ylabel('Phase histogram')
        range = (0, 2.0 * math.pi)
        plt.hist(preset.phase, bins=50, range=range)
        plt.xlim(range[0], range[1])
        save_plot(preset_name + '.phase')

        plt.ylabel('Frequency histogram')
        range = (np.percentile(preset.freq, 5), np.percentile(preset.freq, 95))
        d = range[1] - range[0]
        range = (range[0]-0.15*d, range[1]+0.15*d)
        plt.hist(preset.freq, bins=50, range=range)
        plt.xlim(range[0], range[1])
        save_plot(preset_name + '.freq')
