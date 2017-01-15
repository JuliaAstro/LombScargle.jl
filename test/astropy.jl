### Compare LombScargle.jl with Astropy

using LombScargle, Base.Test, PyCall
PyCall.@pyimport astropy.stats as ast

srand(1)
ntimes = 401

# Un-evenly spaced data.
t = linspace(0.01, 10pi, ntimes)
t += step(t)*rand(ntimes)
for f in (x -> sinpi.(x), x -> sin.(x) + 1.5*cospi.(4*x) + 3)
    s = f(t)
    # "psd" normalization in LombScargle slightly differ from that of
    # Astropy and the test would fail if we includ it.
    for fitmean in (true, false), nrm in ("standard", "model", "log"),
        fast in ((true, "fast"), (false, "cython")), center in (true, false)
        f_jl, p_jl = freqpower(lombscargle(t, s, fit_mean = fitmean,
                                           center_data = center,
                                           normalization=nrm,
                                           maximum_frequency=20, fast = fast[1]))
        f_py, p_py =
            ast.LombScargle(t, s,
                            fit_mean = fitmean,
                            center_data = center)[:autopower](method=fast[2],
                                                              normalization=nrm,
                                                              maximum_frequency=20)
        @test f_jl ≈ f_py
        @test p_jl ≈ p_py
    end
end

# Evenly spaced data.  Use both heteroskedastic and homoskedastic uncertainties.
t = linspace(0.01, 10pi, ntimes)
errors = rand(0.1:1e-3:4.0, ntimes)
for f in (x -> sinpi.(x), x -> sin.(x) + 1.5*cospi.(4*x) + 3), err in (ones(ntimes), errors)
    s = f(t)
    # "psd" normalization in LombScargle slightly differ from that of
    # Astropy and the test would fail if we includ it.
    for fitmean in (true, false), nrm in ("standard", "model", "log"),
        fast in ((true, "fast"), (false, "cython")), center in (true, false)
        f_jl, p_jl = freqpower(lombscargle(t, s, err,
                                           fast = fast[1],
                                           fit_mean = fitmean,
                                           center_data = center,
                                           normalization=nrm,
                                           maximum_frequency=10,
                                           samples_per_peak=10))
        f_py, p_py =
            ast.LombScargle(t, s, dy = err,
                            fit_mean = fitmean,
                            center_data = center)[:autopower](method=fast[2],
                                                              normalization=nrm,
                                                              maximum_frequency=10,
                                                              samples_per_peak=10)
        @test f_jl ≈ f_py
        @test p_jl ≈ p_py
    end
end

# Test heteroskedastic uncertainties with non-fast method.  Test only standard
# normalization, the only one for which LombScargle and Astropy give similar
# results.
for f in (x -> sinpi.(x), x -> sin.(x) + 1.5*cospi.(4*x) + 3)
    s = f(t)
    for fitmean in (true, false), nrm in ("standard", "model"),
        fast in ((true, "fast"), (false, "cython")), center in (true, false)
        f_jl, p_jl = freqpower(lombscargle(t, s, errors,
                                           fit_mean = fitmean,
                                           center_data = center,
                                           fast = fast[1],
                                           normalization = nrm,
                                           maximum_frequency=20))
        f_py, p_py =
            ast.LombScargle(t, s, dy = errors,
                            fit_mean = fitmean,
                            center_data = center)[:autopower](method=fast[2],
                                                              normalization = nrm,
                                                              maximum_frequency=20)
        @test f_jl ≈ f_py
        @test p_jl ≈ p_py
    end
end

# Test the model functions
for f in (x -> sinpi.(x), x -> sin.(x) + 1.5*cospi.(4*x) + 3)
    s = f(t)
    for fm in (true, false), cd in (true, false)
        m_jl = LombScargle.model(t, s, 1/2pi, fit_mean=fm, center_data=cd)
        m_py = ast.LombScargle(t, s, center_data=cd, fit_mean=fm)[:model](t, 1/2pi)
        @test m_jl ≈ m_py
    end
end

info("LombScargle.jl is consistent with Astropy")
