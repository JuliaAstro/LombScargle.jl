using LombScargle
using Measurements
using Base.Test

ntimes = 1001
# Observation times
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += step(t)*rand(ntimes)
# The signal
s = sinpi(t) + cospi(2t) + rand(ntimes)
# Frequency grid
nfreqs = 10000
freqs = linspace(0.01, 3, nfreqs)
# Randomize frequency grid
freqs += step(freqs)*rand(nfreqs)
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
# Test signal with uncertainties
err = linspace(0.5, 1.5, ntimes)
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fit_mean=true))  [0.09230959166317665,0.00156640109976925,  0.0001970465924587832,6.331573873913458e-5,3.794844882537295e-5]
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fit_mean=false)) [0.02988686776042212,0.0005456197937194695,1.9125076826683257e-5,4.542583863304549e-6,1.0238340733199874e-5]
@test power(lombscargle(t, s, err)) ==
    power(lombscargle(t, measurement(s, err)))

# Test autofrequency function
@test_approx_eq LombScargle.autofrequency(t)                       0.003184112396292367:0.006368224792584734:79.6824127172165
@test_approx_eq LombScargle.autofrequency(t, minimum_frequency=0)                   0.0:0.006368224792584734:79.6792286048202
@test_approx_eq LombScargle.autofrequency(t, maximum_frequency=10) 0.003184112396292367:0.006368224792584734:9.99492881196174

### Compare with Astropy
using PyCall
# This test fails on AppVeyor for some strange reason (maybe problems with
# PyCall and Julia 0.5), just skip it on Windows.
if is_unix()
    t = linspace(0.01, 10pi, ntimes)
    t += step(t)*rand(ntimes)
    for f in (x -> sinpi(x), x -> sin(x) + 1.5*cospi(4*x) + 3)
        s = f(t)
        for fitmean in (true, false)
            f_jl, p_jl = freqpower(lombscargle(t, s, fit_mean = fitmean))
            f_py, p_py =
                # Astropy is necessary for this test to be meaningful.  We don't
                # require it, but we have to guard against its absence.
                try
                    (PyCall.@pyimport astropy.stats as ast;
                     ast.LombScargle(t, s, fit_mean = fitmean)[:autopower](method="cython"))
                catch
                    f_jl, p_jl
                end
            # In some cases f_jl and p_jl are one-element longer than f_py and p_py
            @test_approx_eq f_jl[1:length(f_py)] f_py
            @test_approx_eq p_jl[1:length(p_py)] p_py
        end
    end
end
