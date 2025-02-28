Usage
=====

After installing the package, you can start using it with

```julia
using LombScargle
```

The module defines a new `LombScargle.Periodogram` data type, which,
however, is not exported because you will most probably not need to
directly manipulate such objects. This data type holds both the
frequency and the power vectors of the periodogram.

The main function provided by the package is `lombscargle`:

```@docs
lombscargle(::AbstractVector{<:Real}, rest...)
```

`lombscargle` returns a `LombScargle.Periodogram`. The only two mandatory
arguments are:

- `times`: the vector of observation times
- `signal`: the vector of observations associated with `times`

The optional argument is:

- `errors`: the uncertainties associated to each `signal` point.

All these vectors must have the same length.

!!! tip

    You can pre-plan a periodogram with [`LombScargle.plan`](@ref)
    function, which has the same syntax as [`lombscargle`](@ref)
    described in this section. In this way the actual computation of the
    periodogram is faster and you will save memory. See the [Planning the
    Periodogram](#Planning-the-Periodogram-1) section below.

!!! tip

    `LombScargle.jl` exploits Julia's native
    [multi-threading](http://docs.julialang.org/env1/manual/parallel-computing/#multi-threading-experimental)
    for the non-fast methods (the methods used when you set the keyword
    `fast=false`). Run Julia with ``n`` threads (e.g., `JULIA_NUM_THREADS=4 julia` for
    4 threads, if your machine has 4 physical cores) in order to automatically gain
    an ``n`` -fold scaling.

    Please note that multi-threading is still an experimental feature in Julia, so
    you may encounter issues when running it with more than one thread. For example,
    bug [#17395](https://github.com/JuliaLang/julia/issues/17395) on older versions of Julia
    may prevent the function, on some systems, from effectively scaling.

If the signal has uncertainties, the `signal` vector can also be a vector of
`Measurement` objects (from
[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) package), in
which case you need not to pass a separate `errors` vector for the uncertainties
of the signal. You can create arrays of `Measurement` objects with the
`measurement` function, see `Measurements.jl` manual at
<https://juliaphysics.github.io/Measurements.jl/stable> for more details. The
generalised Lomb--Scargle periodogram by [ZK09](@citet) is always used when the
signal has uncertainties, because the original Lomb--Scargle algorithm cannot
handle them.

!!! tip

    The uncertainties are only used in the generalised Lomb--Scargle algorithm to
    build an
    [inverse-variance](https://en.wikipedia.org/wiki/Inverse-variance_weighting)
    weights vector (see [ZK09](@citet)), that gives more importance to
    datapoints with lower uncertainties. The case where all measurements have the
    same uncertainty (a condition known as
    [homoskedasticity](https://en.wikipedia.org/wiki/Homoscedasticity)) results in a
    constant weights vector, like if there are no uncertainties at all. If you have
    homoskedastic errors, you do not need to provide them to
    [`lombscargle`](@ref).


Planning the Periodogram
------------------------

In a manner similar to planning Fourier transforms with FFTW, it is possible to
speed-up computation of the Lomb--Scargle periodogram by pre-planning it with
[`LombScargle.plan`](@ref) function. It has the same syntax as
[`lombscargle`](@ref), which in the base case is:

```@docs
LombScargle.plan
LombScargle.autofrequency
```

`LombScargle.plan` takes all the same argument as [`lombscargle`](@ref) shown
above and returns a `LombScargle.PeriodogramPlan` object after having
pre-computed certain quantities needed afterwards, and pre-allocated the memory
for the periodogram. It is highly suggested to plan a periodogram before
actually computing it, especially for the fast method. Once you plan a
periodogram, you can pass the `LombScargle.PeriodogramPlan` to
[`lombscargle`](@ref) as the only argument.

```@docs
lombscargle(::LombScargle.PeriodogramPlan)
```

Planning the periodogram has a twofold advantage. First of all, the planning
stage is
[type-unstable](https://docs.julialang.org/en/v1/manual/performance-tips),
because the type of the plan depends on the value of input parameters, and not
on their types. Thus, separating the planning (inherently inefficient) from the
actual computation of the periodogram (completely type-stable) makes overall
computation faster than directly calling [`lombscargle`](@ref). Secondly, the
`LombScargle.PeriodogramPlan` bears the time vector, but the quantities that are
pre-computed in planning stage do not actually depend on it. This is
particularly useful if you want to calculate the [False-Alarm Probability](@ref)
via bootstrapping with [`LombScargle.bootstrap`](@ref): the vector time is
randomly shuffled, but pre-computed quantities will remain the same, saving both
time and memory in each iteration. In addition, you ensure that you will use the
same options you used to compute the periodogram.

Fast Algorithm
--------------

When the frequency grid is evenly spaced, you can compute an approximate
generalised Lomb--Scargle periodogram using a fast algorithm proposed by [PR89](@citet)
that greatly speeds up calculations, as it scales as ``O[N \log(M)]`` for ``N``
data points and ``M`` frequencies. For comparison, the true Lomb--Scargle
periodogram has complexity ``O[NM]``.  The larger the number of datapoints, the
more accurate the approximation.

!!! note

    This method internally performs a [Fast Fourier
    Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) to
    compute some quantities, but it is in no way equivalent to conventional Fourier
    periodogram analysis.

    `LombScargle.jl` uses [FFTW](http://fftw.org/) functions to compute the FFT. You
    can speed-up this task by using multi-threading: call `FFTW.set_num_threads(n)`
    to use ``n`` threads. However, please note that the running time will not scale as
    ``n`` because computation of the FFT is only a part of the algorithm.


The only prerequisite in order to be able to employ this fast method is to
provide a `frequencies` vector as an `AbstractRange` object, which ensures that
the frequency grid is perfectly evenly spaced. This is the default, since
`LombScargle.autofrequency` returns an `AbstractRange` object.

!!! tip


    In Julia, an `AbstractRange` object can be constructed for example with the
    [`range`](https://docs.julialang.org/en/v1/base/math/#Base.range) function
    (you specify the start of the range, and optionally the stop, the length and the
    step of the vector) or with the syntax
    [`start:[step:]stop`](https://docs.julialang.org/en/v1/base/math/#Base.::)
    (you specify the start and the end of the range, and optionally the linear
    step).

Since this fast method is accurate only for large datasets, it is enabled by
default only if the number of output frequencies is larger than 200. You can
override the default choice of using this method by setting the `fast` keyword
to `true` or `false`. We recall that in any case, the `frequencies` vector must
be a `Range` in order to use this method.

To summarize, provided that `frequencies` vector is an `AbstractRange` object,
you can use the fast method:

- by default if the length of the output frequency grid is larger than
  200 points
- in any case with the `fast=true` keyword

Setting `fast=false` always ensures you that this method will not be used,
instead `fast=true` actually enables it only if `frequencies` is an
`AbstractRange`.

Normalization
-------------

By default, the periodogram ``p(f)`` is normalized so that it has values in the
range ``0 \leq p(f) \leq 1``, with ``p = 0`` indicating no improvement of the
fit and ``p = 1`` a "perfect" fit (100% reduction of ``\chi^2`` or
``\chi^2 = 0``). This is the normalization suggested by [LOM76](@citet) and
[ZK09](@citet), and corresponds to the `:standard` normalization in
[`lombscargle`](@ref) function. [ZK09](@citet) wrote the formula for the power of the
periodogram at frequency ``f`` as

```math
p(f) = \frac{1}{YY} \left[
\frac{YC^2_{\tau}}{CC_{\tau}} +
\frac{YS^2_{\tau}}{SS_{\tau}}
\right]
```

See the paper for details. The other normalizations for periodograms
``P(f)`` are calculated from this one. In what follows, ``N`` is the number
of observations.

- `:model`:

  ```math
  P(f) = \frac{p(f)}{1 - p(f)}
  ```

- `:log`:

  ```math
  P(f) = -\log(1 - p(f))
  ```

- `:psd`:

  ```math
  P(f) = \frac{W}{2}\left[\frac{YC^2_{\tau}}{CC_{\tau}} +
  \frac{YS^2_{\tau}}{SS_{\tau}}\right] = p(f) \frac{W*YY}{2}
  ```

  where W is the sum of the inverse of the individual errors,
  ``W = \sum \frac{1}{\sigma_{i}}``, as given in [ZK09](@citet).

- `:Scargle`:

  ```math
  P(f) = \frac{p(f)}{\text{noise level}}
  ```

  This normalization can be used when you know the noise level (expected from
  the a priori known noise variance or population variance), but this isn't
  usually the case. See [SCA82](@citet)

- `:HorneBaliunas`:

  ```math
  P(f) = \frac{N - 1}{2} p(f)
  ```

  This is like the `:Scargle` normalization, where the noise has been estimated
  for Gaussian noise to be ``(N - 1)/2``. See [HB86](@citet)

- If the data contains a signal or if errors are under- or overestimated or if
  intrinsic variability is present, then ``(N-1)/2`` may not be a good
  uncorrelated estimator for the noise level.  [CMB99](@citet) suggested to
  estimate the noise level a posteriori with the residuals of the best fit and
  normalised the periodogram as:

  ```math
  P(f) = \frac{N - 3}{2} \frac{p(f)}{1 - p(f_{\text{best}})}
  ```

  This is the `:Cumming` normalization option

Access Frequency Grid and Power Spectrum of the Periodogram
-----------------------------------------------------------

[`lombscargle`](@ref) returns a `LombScargle.Periodogram` object, but you most
probably want to use the frequency grid and the power spectrum. You can access
these vectors with `freq` and `power` functions, just like in `DSP.jl`
package. If you want to get the 2-tuple `(freq(p), power(p))` use the
`freqpower` function.

```@docs
power
freq
freqpower
```

Access Period Grid
------------------

The following utilities are the analogs of [`freq`](@ref) and
[`freqpower`](@ref), but relative to the periods instead of the
frequencies. Thus `period(p)` returns the vector of periods in the periodogram,
that is `1./freq(p)`, and `periodpower(p)` gives you the 2-tuple `(period(p),
power(p))`.

```@docs
period
periodpower
```

### `findmaxpower`, `findmaxfreq`, and `findmaxperiod` Functions

Once you compute the periodogram, you usually want to know which are the
frequencies or periods with highest power. To do this, you can use the
[`findmaxfreq`](@ref) and [`findmaxperiod`](@ref) functions.

```@docs
findmaxpower
findmaxfreq
findmaxperiod
```

False-Alarm Probability
-----------------------

Noise in the data produce fluctuations in the periodogram that will present
several local peaks, but not all of them related to real periodicities. The
significance of the peaks can be tested by calculating the probability that its
power can arise purely from noise.  The higher the value of the power, the lower
will be this probability.

!!! note

    [CMB99](@citet) showed that the different normalizations result
    in different probability functions. `LombScargle.jl` can calculate the
    probability (and the false-alarm probability) only for the normalizations
    reported by [ZK09](@citet), that are `:standard`, `:Scargle`,
    `:HorneBaliunas`, and `:Cumming`.

The probability ``\Pr(p > p_0)`` that the periodogram power ``p`` can
exceed the value ``p_0`` can be calculated with the [`prob`](@ref) function,
whose first argument is the periodogram and the second one is the ``p_0``
value. The function [`probinv`](@ref) is its inverse: it takes the probability
as second argument and returns the corresponding ``p_0`` value.

```@docs
prob(::LombScargle.Periodogram, ::Real)
probinv(::LombScargle.Periodogram, ::Real)
LombScargle.M
fap(::LombScargle.Periodogram, ::Real)
fapinv(::LombScargle.Periodogram, ::Real)
```

Here are the probability functions for each normalization supported by
`LombScargle.jl`:

- `:standard` (``p \in [0, 1]``):

  ```math
  \Pr(p > p_0) = (1 - p_0)^{(N - 3)/2}
  ```

- `:Scargle` (``p \in [0, \infty)``):

  ```math
  \Pr(p > p_0) = \exp(-p_0)
  ```

- `:HorneBaliunas` (``p \in [0, (N - 1)/2]``):

  ```math
  \Pr(p > p_0) = \left(1 - \frac{2p_0}{N - 1}\right)^{(N - 3)/2}
  ```

- `:Cumming` (``p \in [0, \infty)``):

  ```math
  \Pr(p > p_0) = \left(1 + \frac{2p_0}{N - 3}\right)^{-(N - 3)/2}
  ```

As explained by [SS10](@citet), «the term "false-alarm probability denotes the
probability that at least one out of ``M`` independent power values in a
prescribed search band of a power spectrum computed from a white-noise time
series is expected to be as large as or larger than a given
value». `LombScargle.jl` provides the [`fap`](@ref) function to calculate the
false-alarm probability (FAP) of a given power in a periodogram. Its first
argument is the periodogram, the second one is the value ``p_0`` of the power
of which you want to calculate the FAP. The function [`fap`](@ref) uses the
formula

```math
\mathrm{FAP} = 1 - (1 - \Pr(p > p_0))^M
```

where ``M`` is the number of independent frequencies estimated with ``M = T
\cdot \Delta f``, being ``T`` the duration of the observations and ``\Delta f``
the width of the frequency range in which the periodogram has been calculated
(see [CUM04](@citet)). The function [`fapinv`](@ref) is the inverse of
[`fap`](@ref): it takes as second argument the value of the FAP and returns the
corresponding value ``p_0`` of the power.

The detection threshold ``p_0`` is the periodogram power corresponding to some
(small) value of ``\mathrm{FAP}``, i.e. the value of ``p`` exceeded due to noise
alone in only a small fraction ``\mathrm{FAP}`` of trials. An observed power
larger than ``p_0`` indicates that a signal is likely present (see [CUM04](@citet)).

!!! warning

    Some authors stressed that this method to calculate the false-alarm probability
    is not completely reliable. A different approach to calculate the false-alarm
    probability is to perform Monte Carlo or bootstrap simulations in order to
    determine how often a certain power level ``p_0`` is exceeded just by chance
    (see [CMB99](@citet), [CUM04](@citet), and [ZK09](@citet)). See the
    [Bootstrapping](@ref) section.

### Bootstrapping

One of the possible and simplest statistical methods that you can use to measure
the false-alarm probability and its inverse is
[bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29)
(see section 4.2.2 of [MHC93](@citet)).

!!! note

    We emphasize that you can use this method only if you know your data points are
    [independent and identically
    distributed](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables),
    and they have [white uncorrelated noise](https://en.wikipedia.org/wiki/White_noise).

The recipe of the bootstrap method is very simple to implement:

- repeat the Lomb--Scargle analysis a large number ``N`` of times on the original
  data, but with the signal (and errors, if present) vector randomly
  shuffled. As an alternative, shuffle only the time vector;
- out of all these simulations, store the powers of the highest peaks;
- in order to estimate the false-alarm probability of a given power, count how
  many times the highest peak of the simulations exceeds that power, as a
  fraction of ``N``. If you instead want to find the inverse of the false-alarm
  probability ``\text{prob}``, looks for the ``N\cdot\text{prob}``-th element of the
  highest peaks vector sorted in descending order.

Remember to pass to [`lombscargle`](@ref) function the same options, if any, you
used to compute the Lomb--Scargle periodogram before.

`LombScargle.jl` provides simple methods to perform such analysis. The
[`LombScargle.bootstrap`](@ref) function allows you to create a bootstrap sample
with `N` permutations of the original data.

```@docs
LombScargle.bootstrap
```

The false-alarm probability and its inverse can be calculated with [`fap`](@ref)
and [`fapinv`](@ref) functions respectively.  Their syntax is the same as the
methods introduced above, but with a `LombScargle.Bootstrap` object as first
argument, instead of the `LombScargle.Periodogram` one.

```@docs
fap(::LombScargle.Bootstrap{<:AbstractFloat}, ::Real)
fapinv(::LombScargle.Bootstrap{<:AbstractFloat}, ::Real)
```

`LombScargle.model` Function
----------------------------

For each frequency ``f`` (and hence for the corresponding angular frequency
``\omega = 2\pi f``) the Lomb--Scargle algorithm looks for the
sinusoidal function of the type

```math
a_f\cos(\omega t) + b_f\sin(\omega t) + c_f
```

that best fits the data. In the original Lomb--Scargle algorithm the offset
``c`` is null (see [LOM76](@citet)). In order to find the best-fitting
coefficients ``a_f``, ``b_f``, and ``c_f`` for the given frequency ``f``,
without actually performing the periodogram, you can solve the linear system
``\mathbf{A} \mathbf{x} = \mathbf{y}``, where ``\mathbf{A}`` is the matrix

```math
\begin{bmatrix}
  \cos(\omega t) & \sin(\omega t) & 1
\end{bmatrix} =
\begin{bmatrix}
  \cos(\omega t_1) & \sin(\omega t_1) & 1      \\
  \vdots           & \vdots           & \vdots \\
  \cos(\omega t_n) & \sin(\omega t_n) & 1
\end{bmatrix}
```

``t = [t_1, \dots, t_n]^\mathsf{T}`` is the column vector of observation
times, ``\mathbf{x}`` is the column vector with the unknown coefficients

```math
\begin{bmatrix}
    a_f \\
    b_f \\
    c_f
\end{bmatrix}
```

and ``\mathbf{y}`` is the column vector of the signal. The solution of the matrix
gives the wanted coefficients.

This is what the [`LombScargle.model`](@ref) function does in order to return
the best fitting Lomb--Scargle model for the given signal at the given
frequency.

```@docs
LombScargle.model
```
