#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
'''
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

def draw_distribution_plot(distribution, image_file):
    ''' Draw distribution plot.
    '''
    x = []
    y = []
    for k, v in distribution.items():
        x.append(k)
        y.append(v) 
    x = np.asarray(x)
    y = np.asarray(y)
    plt.plot(x,y)
    plt.savefig(image_file)
    plt.clf()