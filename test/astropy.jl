### Compare LombScargle.jl with Astropy

using LombScargle, Base.Test, PyCall
PyCall.@pyimport astropy.stats as ast

ntimes = 201
t = linspace(0.01, 10pi, ntimes)
t += step(t)*rand(ntimes)
for f in (x -> sinpi(x), x -> sin(x) + 1.5*cospi(4*x) + 3)
    s = f(t)
    # "psd" normalization in LombScargle slightly differ from that of
    # Astropy and the test would fail if we includ it.
    for fitmean in (true, false), nrm in ("standard", "model", "log")
        f_jl, p_jl = freqpower(lombscargle(t, s, fit_mean = fitmean,
                                           normalization=nrm,
                                           maximum_frequency=20))
        f_py, p_py =
            ast.LombScargle(t, s,
                            fit_mean = fitmean)[:autopower](method="cython",
                                                            normalization=nrm,
                                                            maximum_frequency=20)
        @test_approx_eq f_jl f_py
        @test_approx_eq p_jl p_py
    end
end

# Test the fast method with evenly spaced data.
t = linspace(0.01, 10pi, ntimes)
for f in (x -> sinpi(x), x -> sin(x) + 1.5*cospi(4*x) + 3)
    s = f(t)
    # "psd" normalization in LombScargle slightly differ from that of
    # Astropy and the test would fail if we includ it.
    for fitmean in (true, ), nrm in ("standard", )
        f_jl, p_jl = freqpower(lombscargle(t, s, fit_mean = fitmean,
                                           normalization=nrm,
                                           maximum_frequency=20))
        f_py, p_py =
            ast.LombScargle(t, s,
                            fit_mean = fitmean)[:autopower](method="fast",
                                                            normalization=nrm,
                                                            maximum_frequency=20)
        @test_approx_eq f_jl f_py
        @test_approx_eq p_jl p_py
    end
end

info("LombScargle.jl is consistent with Astropy")
