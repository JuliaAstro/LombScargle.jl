#!/usr/bin/env python -*- coding: utf-8 -*-

from __future__ import print_function

from timeit import Timer
from numpy import logspace

io = open("python_times.dat", "w")
N = logspace(1, 5, 15).astype(int)
for (i, n) in enumerate(N):
    setup = """
from astropy.timeseries import LombScargle
import numpy as np
t = np.linspace(0.01, 10, {0:d})
s = np.sin(t) + 1.5*np.cos(4 * np.pi * t) + 3
ls = LombScargle(t, s)
""".format(n)
    repeat = 50 if n < 5000 else 20
    print("Iteration {i:d} ({n:d} datapoints)".format(i=i+1, n=n))
    time = min(Timer('ls.autopower()', setup=setup).repeat(repeat, 1))
    print(' Python:', time, 'seconds')
    io.write("{n:d}\t{t:f}\n".format(n=n, t=time))
io.close()
