# LombScargle.jl

| **Documentation**                       | [**Package Evaluator**][pkgeval-link] | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url]       | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

Introduction
------------

`LombScargle.jl` is a [Julia](http://julialang.org/) package for a fast
multi-threaded estimation of
the [frequency spectrum](https://en.wikipedia.org/wiki/Frequency_spectrum) of a
periodic signal
with
[the Lomb–Scargle periodogram](https://en.wikipedia.org/wiki/The_Lomb–Scargle_periodogram).

Another Julia package that provides tools to perform spectral analysis of
signals is [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl), but its methods
require that the signal has been sampled at equally spaced times.  Instead, the
Lomb–Scargle periodogram enables you to analyze unevenly sampled data as well,
which is a fairly common case in astronomy, a field where this periodogram is
widely used.

The algorithms used in this package are reported in the following papers:

* Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL:
  http://dx.doi.org/10.1086/167197, Bibcode:
  http://adsabs.harvard.edu/abs/1989ApJ...338..277P)
* Townsend, R. H. D. 2010, ApJS, 191, 247 (URL:
  http://dx.doi.org/10.1088/0067-0049/191/2/247, Bibcode:
  http://adsabs.harvard.edu/abs/2010ApJS..191..247T)
* Zechmeister, M., Kürster, M. 2009, A&A, 496, 577 (URL:
  http://dx.doi.org/10.1051/0004-6361:200811296, Bibcode:
  http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)

The package provides facilities to:

* compute the periodogram using different methods (with different speeds) and
  different normalizations.  This is one of the fastest implementations of these
  methods available as free software.  If Julia is run with more than
  one
  [thread](http://docs.julialang.org/en/stable/manual/parallel-computing/#multi-threading-experimental),
  computation is automatically multi-threaded, further speeding up calculations;
* access the frequency and period grid of the resulting periodogram, together
  with the power spectrum;
* find the maximum power in the periodogram and the frequency and period
  corresponding to the peak.  All these queries can be restricted to a specified
  region, in order to search a local maximum, instead of the global one;
* calculate the probability that a peak arises from noise only (false-alarm
  probability) using analytic formulas, in order to assess the significance of
  the peak;
* perform bootstrap resamplings in order to compute the false-alarm probability
  with a statistical method;
* determine the best-fitting Lomb–Scargle model for the given data set at the
  given frequency.

All these features are thoroughly described in the full documentation, see
below.  Here we only give basic information.

### Documentation

The complete manual of `LombScargle.jl` is available at
http://lombscarglejl.readthedocs.io.  It has detailed explanation of all
functions provided by the package and more examples than what you will find
here, also with some plots.  You can also download the PDF version of the manual
from https://media.readthedocs.org/pdf/lombscarglejl/latest/lombscarglejl.pdf.

Installation
------------

`LombScargle.jl` is available for Julia 0.6 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the commands

```julia
julia> Pkg.update()
julia> Pkg.add("LombScargle")
```

Older versions are also available for Julia 0.4 and 0.5.

Usage
-----

After installing the package, you can start using it with

```julia
julia> using LombScargle
```

The module defines a new `LombScargle.Periodogram` data type, which, however, is
not exported because you will most probably not need to directly manipulate such
objects.  This data type holds both the frequency and the power vectors of the
periodogram.

The main function provided by the package is `lombscargle`:

```julia
lombscargle(times, signal[, errors])
```

which returns a `LombScargle.Periodogram`.  The only mandatory arguments are:

* `times`: the vector of observation times
* `signal`: the vector of observations associated with `times`

All these vectors must have the same length.  The only optional argument is:

* `errors`: the uncertainties associated to each `signal` point.  This vector
  must have the same length as `times` and `signal`.

Besides the two arguments introduced above, `lombscargle` has a number of other
optional keywords in order to choose the right algorithm to use and tweak the
periodogram.  For the description of all these arguments see the complete
manual.

If the signal has uncertainties, the `signal` vector can also be a vector of
`Measurement` objects (from
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package), in
which case you need not to pass a separate `errors` vector for the uncertainties
of the signal.  You can create arrays of `Measurement` objects with the
`measurement` function, see `Measurements.jl` manual at
http://measurementsjl.readthedocs.io/ for more details.

With the `LombScargle.plan` function you can pre-plan a periodogram and save
time and memory for the actual computation of the periodogram.  See
the [manual](https://lombscarglejl.readthedocs.io/#planning-the-periodogram) for
details.

Examples
--------

Here is an example of a noisy periodic signal (`sin(π*t) + 1.5*cos(2π*t)`)
sampled at unevenly spaced times.

```julia
julia> using LombScargle

julia> ntimes = 1001
1001

# Observation times
julia> t = linspace(0.01, 10pi, ntimes)
0.01:0.03140592653589793:31.41592653589793

# Randomize times
julia> t += step(t)*rand(ntimes);

# The signal
julia> s = sinpi.(t) .+ 1.5 .* cospi.(2t) .+ rand(ntimes);

# Pre-plan the periodogram (see the documentation)
julia> plan = LombScargle.plan(t, s);

# Compute the periodogram
julia> pgram = lombscargle(plan)
```

You can plot the result, for example with
[`Plots`](https://github.com/tbreloff/Plots.jl) package.  Use `freqpower`
function to get the frequency grid and the power of the periodogram as a
2-tuple.

```julia
using Plots
plot(freqpower(pgram)...)
```

### Signal with Uncertainties

The generalised Lomb–Scargle periodogram (used when the `fit_mean` optional
keyword is `true`) is able to handle a signal with uncertainties, and they will
be used as weights in the algorithm.  The uncertainties can be passed either as
the third optional argument `errors` to `lombscargle` or by providing this
function with a `signal` vector of type `Measurement` (from
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package).

```julia
using Measurements, Plots
ntimes = 1001
t = linspace(0.01, 10pi, ntimes)
s = sinpi.(2t)
errors = rand(0.1:1e-3:4.0, ntimes)
plot(freqpower(lombscargle(t, s, errors, maximum_frequency=1.5))...)
plot(freqpower(lombscargle(t, measurement(s, errors), maximum_frequency=1.5))...)
```

Performance
-----------

A pre-planned periodogram in `LombScargle.jl` computed in single thread mode
with the fast method is more than 2.5 times faster than the implementation of
the same algorithm provided by Astropy, and more than 4 times faster if 4 FFTW
threads are used (on machines with at least 4 physical CPUs).

The following plot shows a comparison between the times needed to compute a
periodogram for a signal with N datapoints using `LombScargle.jl`, with 1 or 4
threads (with `flags = FFTW.MEASURE` for better performance), and the
single-threaded Astropy implementation.  (Julia version: 0.6.0-pre.alpha.242,
commit d694548; `LombScargle.jl` version: 0.3.0; Python version: 3.5.3; Astropy
version: 1.3.  CPU: Intel(R) Core(TM) i7-4700MQ.)

![benchmarks](https://cdn.rawgit.com/giordano/LombScargle.jl/master/perf/benchmarks.svg)

Note that this comparison is unfair, as Astropy doesn’t support pre-planning a
periodogram nor exploiting multi-threading.  A non-planned periodogram in single
thread mode in `LombScargle.jl` is still twice faster than Astropy.

Development
-----------

The package is developed at https://github.com/giordano/LombScargle.jl.  There
you can submit bug reports, make suggestions, and propose pull requests.

### History

The ChangeLog of the package is available in
[NEWS.md](https://github.com/giordano/LombScargle.jl/blob/master/NEWS.md) file
in top directory.

License
-------

The `LombScargle.jl` package is licensed under the BSD 3-clause "New" or
"Revised" License.  The original author is Mosè Giordano.

### Acknowledgemets

This package adapts the implementation in Astropy of the the fast Lomb–Scargle
method by Press & Rybicki (1989).  We claim no endorsement nor promotion by the
Astropy Team.



[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://lombscarglejl.readthedocs.io/en/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://lombscarglejl.readthedocs.io/en/stable/

[pkgeval-link]: http://pkg.julialang.org/?pkg=LombScargle

[pkg-0.5-img]: http://pkg.julialang.org/badges/LombScargle_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/detail/LombScargle.html
[pkg-0.6-img]: http://pkg.julialang.org/badges/LombScargle_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/detail/LombScargle.html

[travis-img]: https://travis-ci.org/giordano/LombScargle.jl.svg?branch=master
[travis-url]: https://travis-ci.org/giordano/LombScargle.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/vv6mho713fuse6qy/branch/master?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/lombscargle-jl

[coveral-img]: https://coveralls.io/repos/github/giordano/LombScargle.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/giordano/LombScargle.jl?branch=master

[codecov-img]: https://codecov.io/gh/giordano/LombScargle.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/giordano/LombScargle.jl
