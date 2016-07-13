using LombScargle
using Base.Test

ntimes = 1001
# Observation times
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += (t[2] - t[1])*rand(ntimes)
# The signal
s = sinpi(t) + cospi(2t) + rand(ntimes)
# Frequency grid
nfreqs = 10000
freqs = linspace(0.01, 3, nfreqs)
# Randomize frequency grid
freqs += (freqs[2] - freqs[1])*rand(nfreqs)
# Use "freqpower" just to call that function and increase code coverage.
# "autofrequency" function is tested below.
@test_approx_eq_eps freqpower(lombscargle(t, s, frequencies=freqs, fit_mean=false))[2] freqpower(lombscargle(t, s, frequencies=freqs, fit_mean=true))[2] 5e-3

# Simple signal, without any randomness
t = linspace(0.01, 10pi, ntimes)
s = sin(t)
pgram1 = lombscargle(t, s, fit_mean=false)
pgram2 = lombscargle(t, s, fit_mean=true)
@test_approx_eq_eps power(pgram1) power(pgram2) 2e-7

# Test findmaxfreq
@test_approx_eq findmaxfreq(pgram1)        [31.997145470342]
@test_approx_eq findmaxfreq(pgram1, 0.965) [0.15602150741832602,31.685102455505348,31.997145470342,63.52622641842902,63.838269433265665]

# Test the values in order to prevent wrong results in both algorithms
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=true))  [0.029886871262324886,0.0005456198989410226,1.912507742056023e-5, 4.54258409531214e-6, 1.0238342782430832e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=false)) [0.02988686776042212, 0.0005456197937194695,1.9125076826683257e-5,4.542583863304549e-6,1.0238340733199874e-5]

# Test autofrequency function
@test_approx_eq LombScargle.autofrequency(t)                       0.003184112396292367:0.006368224792584734:79.6824127172165
@test_approx_eq LombScargle.autofrequency(t, minimum_frequency=0)                   0.0:0.006368224792584734:79.6792286048202
@test_approx_eq LombScargle.autofrequency(t, maximum_frequency=10) 0.003184112396292367:0.006368224792584734:9.99492881196174
