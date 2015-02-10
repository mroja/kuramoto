#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import glob
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from common import DataPreset, load_preset_from_file, save_plot


def plot_step(params):
    name = params['name']
    #preset = params['preset']
    step = params['step']
    f_name = params['f_name']
    dir_name = params['dir_name']

    preset = load_preset_from_file(name)

    freq = preset.freq

    with open(f_name, 'r') as f:
        lines = f.readlines()
        
    step_, N, r, mean = (x for x in lines[0].split())
    step_ = int(step_)

    assert(step_ == step)
    
    N = int(N)
    r = float(r)
    mean = float(mean)
    phases = [float(x) for x in lines[1].split()]
    vel = [float(x) for x in lines[2].split()]
    #print len(phases), len(vel)

    print(step)

    #for i in xrange(N):
    #    pos = (phases[i], freq[i])
    #    print pos

    plt.figure()
    plt.suptitle('Step: ' + str(step))
    plt.subplot(2, 1, 1)
    #py.axvline(95)
    #py.axvline(35)
    #plt.xlabel('Phase')
    plt.ylabel('Phase histogram')
    plt.hist(phases, bins=60, range=(0, 2.0 * math.pi))
    plt.xlim(0, 2.0 * math.pi)

    plt.subplot(2, 1, 2)
    #plt.xlabel('Velocity')
    plt.ylabel('Velocity histogram')
    #range = (np.min(vel), np.max(vel))
    range = (-30, 30)
    plt.hist(vel, bins=60, range=range)
    plt.xlim(range[0], range[1])
    save_plot(os.path.join(dir_name, 'hist', str(step)))

    plt.figure()
    plt.title('Step: ' + str(step))
    plt.xlabel('Phase')
    plt.ylabel('Intrinsic frequency')
    plt.xlim(0, 2.0 * math.pi)
    plt.ylim(-3, 3)
    plt.plot(phases, freq, marker='o', ls='')
    save_plot(os.path.join(dir_name, 'phase', str(step)))


def gen_video(dump_dir, subdir_name, framerate):
    pattern = os.path.join(dump_dir, subdir_name, '%d.png')
    out_video = os.path.join(dump_dir, subdir_name + '.avi')
    # TODO: ffmpeg
    cmd = 'avconv -y -start_number 1 -framerate '+str(framerate)+' -i ' + pattern + ' -q:v 1 -vcodec mpeg4 ' + out_video
    #print('Executing: ' + cmd)
    os.system(cmd)


def gen_mean_and_r_plots(dir_name):
    with open(os.path.join(dir_name, 'r.txt')) as f:
        r = [float(x) for x in f.read().split()]
    plt.figure()
    plt.xlabel('Steps')
    plt.ylabel('Order parameter')
    plt.xlim(0, len(r))
    plt.ylim(0, 1)
    plt.plot(range(0, len(r)), r)
    save_plot(os.path.join('dump_' + name, 'r'))

    with open(os.path.join(dir_name, 'mean.txt')) as f:
        mean = [float(x) for x in f.read().split()]
    plt.figure()
    plt.xlabel('Steps')
    plt.ylabel('Mean phase')
    plt.xlim(0, len(mean))
    plt.ylim(0, 2.0 * math.pi)
    plt.plot(range(0, len(mean)), mean)
    save_plot(os.path.join('dump_' + name, 'mean'))

    with open(os.path.join(dir_name, 'mean_vel.txt')) as f:
        mean_vel = [float(x) for x in f.read().split()]
    plt.figure()
    plt.xlabel('Steps')
    plt.ylabel('Mean velocity')
    plt.xlim(0, len(mean_vel))
    plt.plot(range(0, len(mean_vel)), mean_vel)
    save_plot(os.path.join('dump_' + name, 'mean_vel'))


def remove_images(dir_name, remove_dir=True):
    for f in glob.glob(os.path.join(dir_name, '*.png')):
        os.remove(f)
    if remove_dir:
        try:
            os.rmdir(dir_name)
        except OSError as e:
            print('Cannot remove directory: ' + dir_name + ' (' + str(e) + ')')


def remove_step_files(dump_dir):
    for f in glob.glob(os.path.join(dump_dir, '*.txt')):
        os.remove(f)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('Usage: gen_plots.py name')
        sys.exit()

    name = sys.argv[1]
    dir_name = 'dump_' + name
    steps_dir = os.path.join(dir_name, 'steps')

    # read sorted list of states at specific steps
    step_files_all = glob.glob(os.path.join(steps_dir, '*.txt'))
    
    def filter_files(seq):
        for el in seq:
            name = os.path.basename(el).replace('.txt', '')
            if 'r' not in name and 'mean' not in name:
                yield el
    step_files = [f for f in filter_files(step_files_all)]
    input_files = [(int(os.path.basename(f).replace('.txt', '')), f) for f in step_files]
    input_files.sort(key=lambda x: x[0])

    # take every M-th snapshot
    M = 1
    input_files = input_files[::M]

    gen_mean_and_r_plots(steps_dir)

    if 1:
        remove_images(os.path.join(dir_name, 'hist'), remove_dir=False)
        remove_images(os.path.join(dir_name, 'phase'), remove_dir=False)

        ctx = multiprocessing.get_context('spawn')
        pool = ctx.Pool(multiprocessing.cpu_count())

        args = []
        for step, f_name in input_files:
            args.append({
                'name': name,
                'step': step,
                'f_name': f_name,
                'dir_name': dir_name
            })
        #print(args)
        pool.map(plot_step, args)
        pool.close()

        # rename step numbers to consequent integers
        # this is required for video generation step
        plot_num = 1
        for step, f_name in input_files:
            # print plot_num, step
            for x in ['hist', 'phase']:
                os.rename(
                    os.path.join(dir_name, x, str(step) + '.png'),
                    os.path.join(dir_name, x, str(plot_num) + '.png')
                )
            plot_num += 1

        framerate = 8
        gen_video(dir_name, 'hist', framerate)
        gen_video(dir_name, 'phase', framerate)

        remove_images(os.path.join(dir_name, 'hist'), remove_dir=True)
        remove_images(os.path.join(dir_name, 'phase'), remove_dir=True)

    #remove_step_files(dir_name)
