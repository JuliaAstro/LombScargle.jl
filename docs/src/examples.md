Examples
========

```@meta
DocTestSetup = quote
    using LombScargle
end
```

Here is an example of a noisy periodic signal (``\sin(\pi t) + 1.5\cos(2\pi t)``)
sampled at unevenly spaced times.

```@repl example1
using LombScargle
using Random
Random.seed!(42); # make plots reproducible
ntimes = 1001;
t = range(0.01, 10pi, length = ntimes) # Observation times
t += step(t)*rand(ntimes); # Randomize times
s = sinpi.(t) .+ 1.5cospi.(2t) .+ rand(ntimes); # The signal
plan = LombScargle.plan(t, s); # Pre-plan the periodogram
pgram = lombscargle(plan) # Compute the periodogram
```

You can plot the result, for example with
[Plots](https://github.com/tbreloff/Plots.jl) package. Use [`freqpower`](@ref)
function to get the frequency grid and the power of the periodogram as a
2-tuple.

```@example example1
using Plots
gr(size = (1200, 800)) # hide
plot(freqpower(pgram)...)
xlims!(0, 2.0) # hide
xlabel!("Frequency") # hide
ylabel!("Lomb–Scargle power") # hide
```

You can also plot the power vs the period, instead of the frequency, with
[`periodpower`](@ref):

```@example example1
using Plots
gr(size = (1200, 800)) # hide
plot(periodpower(pgram)...)
xlims!(0.5, 2.5) # hide
xlabel!("Period")
ylabel!("Lomb–Scargle power")
```

!!! warning
    If you do not fit for the mean of the signal (`fit_mean=false` keyword to
    [`lombscargle`](@ref) function) without centering the data (`center_data=false`)
    you can get inaccurate results. For example, spurious peaks at low frequencies
    can appear and the real peaks lose power:

    ```@example example1
    plot(freqpower(lombscargle(t, s, fit_mean=false, center_data=false))...)
    xlims!(0, 2.0) # hide
    ```


!!! tip

    You can tune the frequency grid with appropriate keywords to
    [`lombscargle`](@ref) function. For example, in order to increase the sampling
    increase `samples_per_peak`, and set `maximum_frequency` to lower values in
    order to narrow the frequency range:

    ```@example example1
    plot(freqpower(lombscargle(t, s, samples_per_peak=20, maximum_frequency=1.5))...)
    ```


    If you simply want to use your own frequency grid, directly set the
    `frequencies` keyword:

    ```@example example1
    plot(freqpower(lombscargle(t, s, frequencies=0.001:1e-3:1.5))...)
    ```


Signal with Uncertainties
-------------------------

The generalised Lomb--Scargle periodogram is able to handle a signal with
uncertainties, and they will be used as weights in the algorithm.  The
uncertainties can be passed either as the third optional argument `errors` to
[`lombscargle`](@ref) or by providing this function with a `signal` vector of
type `Measurement` (from
[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) package).

```@example uncertainties
using LombScargle # hide
using Measurements, Plots
gr(size = (1200, 800)) # hide
ntimes = 1001
t = range(0.01, stop = 10pi, length = ntimes)
s = sinpi.(2t)
errors = rand(0.1:1e-3:4.0, ntimes)
# Run one of the two following equivalent commands
plot(freqpower(lombscargle(t, s, errors, maximum_frequency=1.5))...)
plot(freqpower(lombscargle(t, measurement.(s, errors), maximum_frequency=1.5))...)
xlims!(0.25, 1.5) # hide
```

This is the plot of the power versus the period:

```julia
# Run one of the two following equivalent commands
plot(periodpower(lombscargle(t, s, errors, maximum_frequency=1.5))...)
plot(periodpower(lombscargle(t, measurement(s, errors), maximum_frequency=1.5))...)
```
```@example uncertainties
plot(periodpower(lombscargle(t, s, errors, maximum_frequency=1.5))...) # hide
xlims!(0.75, 2) # hide
```

We recall that the generalised Lomb--Scargle algorithm is used when the
`fit_mean` optional keyword to [`lombscargle`](@ref) is `true` if no error is
provided, instead it is always used if the signal has uncertainties.

Find Highest Power and Associated Frequencies and Periods
---------------------------------------------------------

[`findmaxfreq`](@ref) function tells you the frequencies with the highest power
in the periodogram (and you can get the period by taking its inverse):

```jldoctest demo1
julia> t = range(0, stop = 10, length = 1001);

julia> s = sinpi.(t);

julia> plan = LombScargle.plan(t, s); # Plan the periodogram

julia> p = lombscargle(plan);

julia> findmaxperiod(p) # Period with highest power
1-element Vector{Float64}:
 0.004987779939149084

julia> findmaxfreq(p) # Frequency with the highest power
1-element Vector{Float64}:
 200.49
```

This peak is at high frequencies, very far from the expected value of the period
of 2. In order to find the real peak, you can either narrow the ranges in order
to exclude higher armonics

```jldoctest demo1
julia> findmaxperiod(p, [1, 10]) # Limit the search to periods in [1, 10]
1-element Vector{Float64}:
 2.0408163265306123

julia> findmaxfreq(p, [0.1, 1]) # Limit the search to frequencies in [0.1, 1]
1-element Vector{Float64}:
 0.49
```

or pass the `threshold` argument to [`findmaxfreq`](@ref) or
[`findmaxperiod`](@ref). You can use [`findmaxpower`](@ref) to discover the
highest power in the periodogram:

```jldoctest demo1
julia> findmaxpower(p)
0.9996235276303144

julia> findmaxperiod(p, 0.95)
10-element Vector{Float64}:
 2.0408163265306123
 1.9607843137254901
 0.010051261433309882
 0.010049241282283187
 0.009951238929246691
 0.009949258780220873
 0.005012782595618828
 0.005012280086211218
 0.004987779939149084
 0.004987282429804

julia> findmaxfreq(p, 0.95)
10-element Vector{Float64}:
   0.49
   0.51
  99.49
  99.51
 100.49
 100.51
 199.49
 199.51
 200.49
 200.51
```

The first peak is the real one, the other double peaks appear at higher
armonics.

!!! tip

    Usually, plotting the periodogram can give you a clearer idea of what's going on.


Significance of the Peaks
-------------------------

The significance of the peaks in the Lomb--Scargle periodogram can be assessed
by measuring the [False-Alarm Probability](#false-alarm-probability). Analytic
expressions of this quantity and its inverse can be obtained with the
[`fap`](@ref) and [`fapinv`](@ref) functions, respectively.

```julia-repl
julia> t = linspace(0.01, 20, samples_per_peak = 10);

julia> s = sinpi.(e.*t).^2 .- cos.(5t).^4;

julia> plan = LombScargle.plan(t, s);

julia> p = lombscargle(plan);

# Find the false-alarm probability for the highest peak.
julia> fap(p, 0.3)
0.028198095962262748
```

Thus, a peak with power ``0.3`` has a probability of ``0.028`` that it is due to
noise only. A quantity that is often used is the inverse of the false-alarm
probability as well: what is the minimum power whose false-alarm probability is
lower than the given probability? For example, if you want to know the minimum
power for which the false-alarm probability is at most ``0.01`` you can use:

```julia-repl
julia> fapinv(p, 0.01)
0.3304696923786712
```

As we already noted, analytic expressions of the false-alarm probability and its
inverse may not be reliable if your data does not satisfy specific
assumptions. A better way to calculate this quantity is to use statistical
methods. One of this is bootstrapping. In `LombScargle.jl`, you can use the
function [`LombScargle.bootstrap`](@ref) to create a bootstrap sample and then
you can calculate the false-alarm probability and its inverse using this sample.

!!! tip

    When applying the bootstrap method you should use the same options you used to
    perform the periodogram on your data. Using the same periodogram plan you used
    to compute the periodogram will ensure that you use the same options. However,
    note that the fast method gives approximate results that for some frequencies
    may not be reliable (they can go outside the range ``[0, 1]`` for the standard
    normalization). More robust results can be obtained with the `fast = false`
    option.

Create a bootstrap sample with 10000 resamplings of the original data,
re-using the same periodogram plan.  The larger the better.
This may take some minutes.
```julia-repl
julia> b = LombScargle.bootstrap(10000, plan)
```

Calculate the false-alarm probability of a peak
with power 0.3 using this bootstrap sample.
```julia-repl
julia> fap(b, 0.3)
0.0209
```

Calculate the lowest power that has probability
less than 0.01 in this bootstrap sample.
```julia-repl
julia> fapinv(b, 0.01)
0.3268290388848437
```

If you query [`fapinv`](@ref) with a too low probability, the corresponding
power cannot be determined and you will get `NaN` as result.

```julia-repl
julia> fapinv(b, 1e-5)
NaN
```

If you want to find the power corresponding to a false-alarm probability of
``\text{prob} = 10^{-5}``, you have to create a new bootstrap sample with ``N``
resamplings so that ``N\cdot\text{prob}`` can be rounded to an integer larger than
or equal to one (for example ``N = 10^{5}``).

Find the Best-Fitting Model
---------------------------

The [`LombScargle.model`](@ref) function can help you to test whether a certain
frequency fits well your data.

```@example
using LombScargle
using Random
using Plots
gr(size = (1200, 800)) # hide
Random.seed!(42)
t = range(0.01, stop = 10pi, length = 1000) # Observation times
s = sinpi.(t) .+ 1.2cospi.(t) .+ 0.3rand(length(t)) # The noisy signal
# Pick-up the best frequency
f = findmaxfreq(lombscargle(t, s, maximum_frequency=10, samples_per_peak=20))[1]
t_fit = range(0, stop = 1, length = 50)
s_fit = LombScargle.model(t, s, f, t_fit/f) # Determine the model
scatter(mod.(t.*f, 1), s, lab="Phased data", title="Best Lomb-Scargle frequency: $f")
plot!(t_fit, s_fit, lab="Best-fitting model", linewidth=4)
```

!!! tip
    If there are more than one dominant frequency you may need to consider more
    models. This task may require some work and patience. Plot the periodogram in
    order to find the best frequencies.

    ```@example
    using LombScargle
    using Random
    using Plots
    gr(size = (1200, 800)) # hide
    Random.seed!(42)
    t = range(0.01, stop = 5, length = 1000) # Observation times
    s = sinpi.(2t) .+ 1.2cospi.(4t) .+ 0.3rand(length(t)) # Noisy signal
    plan = LombScargle.plan(t, s, samples_per_peak=50)
    p = lombscargle(plan)
    # After plotting the periodogram, you discover
    # that it has two prominent peaks around 1 and 2.
    f1 = findmaxfreq(p, [0.8, 1.2])[1] # Get peak frequency around 1
    f2 = findmaxfreq(p, [1.8, 2.2])[1] # Get peak frequency around 2
    fit1 = LombScargle.model(t, s, f1) # Determine the first model
    fit2 = LombScargle.model(t, s, f2) # Determine the second model
    scatter(t, s, lab="Data", title="Best-fitting Lomb-Scargle model")
    plot!(t, fit1 + fit2, lab="Best-fitting model", linewidth=4)
    ```
