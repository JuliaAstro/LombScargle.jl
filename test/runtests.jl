using LombScargle
using Base.Test

ntimes = 1001
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += (t[2] - t[1])*rand(ntimes)
# The signal
s = sinpi(t) + cospi(2t) + rand(ntimes)
# Frequency grid
freqs = linspace(0.01, 3, 1e4)
@test_approx_eq_eps power(lombscargle(t, s, freqs, fit_mean=false)) power(lombscargle(t, s, freqs, fit_mean=true)) 5e-3

# Test autofrequency function
@test_approx_eq autofrequency(linspace(0.01, 10pi, ntimes))                       0.003184112396292367:0.006368224792584734:79.6824127172165
@test_approx_eq autofrequency(linspace(0.01, 10pi, ntimes), minimum_frequency=0)                   0.0:0.006368224792584734:79.6792286048202
@test_approx_eq autofrequency(linspace(0.01, 10pi, ntimes), maximum_frequency=10) 0.003184112396292367:0.006368224792584734:9.99492881196174
