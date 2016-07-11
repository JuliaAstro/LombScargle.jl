using LombScargle
using Base.Test

ntimes = 1001
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += (t[2] - t[1])*rand(ntimes)
# The signal
s = sinpi(t) + cospi(2t) + rand(length(t))
# Frequency grid
freqs = linspace(0.01, 3, 1e4)
@test isapprox(power(lombscargle(t, s, freqs, fit_mean=false)),
               power(lombscargle(t, s, freqs, fit_mean=true)),
               atol=5e-3)
