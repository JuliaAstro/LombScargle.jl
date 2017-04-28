#!/usr/bin/env julia

### images.jl
#
# Copyright (C) 2016 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: periodogram, lomb scargle
#
# This file is a part of LombScargle.jl.
#
# License is BSD 3-clause "New" or "Revised".
#
### Commentary:
#
# This file produces images used in documentation.
#
### Code:

srand(42)

using Measurements, LombScargle, Plots

sz = (1200, 800)

# First periodograms
ntimes = 1001
t = linspace(0.01, 10pi, ntimes)
t += step(t)*rand(ntimes)
s = sinpi(t) + 1.5cospi(2t) + rand(ntimes)
p = lombscargle(t, s)
plot(freqpower(p)..., size = sz, xlim = (0.0, 2.0), xlabel = "Frequency",
     ylabel = "Lomb–Scargle power")
savefig("freq-periodogram.png")
plot(periodpower(p)..., size = sz, xlim = (0.5, 2.5), xlabel = "Period",
     ylabel = "Lomb–Scargle power")
savefig("period-periodogram.png")

# signal with uncertainties
ntimes = 1001
t = linspace(0.01, 10pi, ntimes)
s = sinpi(2t)
errors = rand(0.1:1e-3:4.0, ntimes)
p = lombscargle(t, s, errors, maximum_frequency=1.5)
plot(freqpower(p)..., size = sz, xlim = (0.25, 1.5), xlabel = "Frequency",
     ylabel = "Lomb–Scargle power")
savefig("freq-uncertainties.png")
plot(periodpower(p)..., size = sz, xlim = (0.75, 2.0), xlabel = "Period",
     ylabel = "Lomb–Scargle power")
savefig("period-uncertainties.png")
