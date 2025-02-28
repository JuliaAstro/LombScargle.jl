LombScargle.jl
==============

```@meta
DocTestSetup = quote
    using LombScargle
end
```

Introduction
------------

[`LombScargle.jl`](https://github.com/JuliaAstro/LombScargle.jl) is a package for
a fast multi-threaded estimation of the [frequency
spectrum](https://en.wikipedia.org/wiki/Frequency_spectrum) of a periodic signal
with [the Lomb--Scargle
periodogram](https://en.wikipedia.org/wiki/The_Lomb–Scargle_periodogram).  This
is written in [Julia](http://julialang.org/), a modern high-level,
high-performance dynamic programming language designed for technical computing.

Another Julia package that provides tools to perform spectral analysis of
signals is [DSP.jl](https://github.com/JuliaDSP/DSP.jl), but its methods require
that the signal has been sampled at equally spaced times. Instead, the
Lomb--Scargle periodogram enables you to analyze unevenly sampled data as well,
which is a fairly common case in astronomy, a field where this periodogram is
widely used.

The algorithms used in this package are reported in the following papers:

```@bibliography
Pages = []
Canonical = false
PR89
TOW10
ZK09
```

Other relevant papers are:

```@bibliography
Pages = []
Canonical = false
CMB99
CUM04
HB86
LOM76
MHC93
SCA82
SS10
```

The package provides facilities to:

- compute the periodogram using different methods (with different
  speeds) and different normalizations. This is one of the fastest
  implementations of these methods available as free software. If
  Julia is run with more than one
  [thread](http://docs.julialang.org/en/stable/manual/parallel-computing/#multi-threading-experimental),
  computation is automatically multi-threaded, further speeding up
  calculations;
- access the frequency and period grid of the resulting periodogram,
  together with the power spectrum;
- find the maximum power in the periodogram and the frequency and
  period corresponding to the peak. All these queries can be
  restricted to a specified region, in order to search a local
  maximum, instead of the global one;
- calculate the probability that a peak arises from noise only
  (false-alarm probability) using analytic formulas, in order to
  assess the significance of the peak;
- perform bootstrap resamplings in order to compute the false-alarm
  probability with a statistical method;
- determine the best-fitting Lomb--Scargle model for the given data
  set at the given frequency.


Contents
--------
```@contents
Pages = ["index.md"]
Depth = 4
```


Installation
------------

`LombScargle.jl` is available for Julia 1.0 and later versions, and can
be installed with Julia's built-in
[package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the commands

```julia-repl
julia> import Pkg
julia> Pkg.add("LombScargle")
```

Older versions are also available for Julia 0.4-0.6.

For instructions on using the package, please see the [Usage page](@ref Usage).


Performance
-----------

A pre-planned periodogram in `LombScargle.jl` computed in single thread mode
with the fast method is more than 2 times faster than the implementation of the
same algorithm provided by AstroPy, and more than 4 times faster if 4 FFTW
threads are used (on machines with at least 4 physical CPUs).

The following plot shows a comparison between the times needed to compute a
periodogram for a signal with N datapoints using `LombScargle.jl`, with 1 or 4
FFTW threads (with `flags = FFTW.MEASURE` for better performance), and the
single-threaded Astropy implementation.  (Julia version: 1.6.0; `LombScargle.jl`
version: 1.0.0; Python version: 3.8.6; Astropy version: 4.1.  CPU: Intel(R)
Core(TM) i7-4870HQ CPU @ 2.50GHz.)

![image](benchmarks.png)

Note that this comparison is unfair, as Astropy doesn’t support pre-planning a
periodogram nor multi-threading, and it pads vectors for FFT to a length which
is a power of 2, while by default `LombScargle.jl` uses length which are
multiples of 2, 3, 5, 7.  A non-planned periodogram in single thread mode in
`LombScargle.jl` is still twice as fast as Astropy.

Development
-----------

The package is developed at
<https://github.com/JuliaAstro/LombScargle.jl>. There you can submit bug
reports, make suggestions, and propose pull requests.

### History

The ChangeLog of the package is available in
[NEWS.md](https://github.com/JuliaAstro/LombScargle.jl/blob/master/NEWS.md) file
in top directory.

License
-------

The `LombScargle.jl` package is licensed under the BSD 3-clause "New" or
"Revised" License. The original author is Mosè Giordano.

### Acknowledgements

This package adapts the implementation in Astropy of the the fast Lomb--Scargle
method by [PR89](@citet). We claim no endorsement nor promotion by the Astropy Team.

References
----------
```@bibliography
```
