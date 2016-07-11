using LombScargle
using Base.Test

x = linspace(0.01, 10pi, 1e3)
y = sinpi(x) + cospi(2x) + rand(length(x))
freqs = linspace(0.01, 3, 1e4)
@test isapprox(power(lombscargle(x, y, freqs, fit_mean=false)),
               power(lombscargle(x, y, freqs, fit_mean=true)),
               atol=5e-3)
