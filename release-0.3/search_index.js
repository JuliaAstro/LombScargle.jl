var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "LombScargle.jl",
    "title": "LombScargle.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#LombScargle.jl-1",
    "page": "LombScargle.jl",
    "title": "LombScargle.jl",
    "category": "section",
    "text": "DocTestSetup = quote\n    using LombScargle\nend"
},

{
    "location": "index.html#Introduction-1",
    "page": "LombScargle.jl",
    "title": "Introduction",
    "category": "section",
    "text": "LombScargle.jl is a package for a fast multi-threaded estimation of the frequency spectrum of a periodic signal with the Lomb–Scargle periodogram.  This is written in Julia, a modern high-level, high-performance dynamic programming language designed for technical computing.Another Julia package that provides tools to perform spectral analysis of signals is DSP.jl, but its methods require that the signal has been sampled at equally spaced times. Instead, the Lomb–Scargle periodogram enables you to analyze unevenly sampled data as well, which is a fairly common case in astronomy, a field where this periodogram is widely used.The algorithms used in this package are reported in the following papers:[PR89] Press, W. H., Rybicki, G. B. 1989, ApJ, 338, 277 (URL: http://dx.doi.org/10.1086/167197, Bibcode: http://adsabs.harvard.edu/abs/1989ApJ...338..277P)\n[TOW10] Townsend, R. H. D. 2010, ApJS, 191, 247 (URL: http://dx.doi.org/10.1088/0067-0049/191/2/247, Bibcode: http://adsabs.harvard.edu/abs/2010ApJS..191..247T)\n[ZK09] Zechmeister, M., Kürster, M. 2009, A&A, 496, 577 (URL: http://dx.doi.org/10.1051/0004-6361:200811296, Bibcode: http://adsabs.harvard.edu/abs/2009A%26A...496..577Z)Other relevant papers are:[CMB99] Cumming, A., Marcy, G. W., & Butler, R. P. 1999, ApJ, 526, 890 (URL: http://dx.doi.org/10.1086/308020, Bibcode: http://adsabs.harvard.edu/abs/1999ApJ...526..890C)\n[CUM04] Cumming, A. 2004, MNRAS, 354, 1165 (URL: http://dx.doi.org/10.1111/j.1365-2966.2004.08275.x, Bibcode: http://adsabs.harvard.edu/abs/2004MNRAS.354.1165C)\n[HB86] Horne, J. H., & Baliunas, S. L. 1986, ApJ, 302, 757 (URL: http://dx.doi.org/10.1086/164037, Bibcode: http://adsabs.harvard.edu/abs/1986ApJ...302..757H)\n[LOM76] Lomb, N. R. 1976, Ap&SS, 39, 447 (URL: http://dx.doi.org/10.1007/BF00648343, Bibcode: http://adsabs.harvard.edu/abs/1976Ap%26SS..39..447L)\n[MHC93] Murdoch, K. A., Hearnshaw, J. B., & Clark, M. 1993, ApJ, 413, 349 (URL: http://dx.doi.org/10.1086/173003, Bibcode: http://adsabs.harvard.edu/abs/1993ApJ...413..349M)\n[SCA82] Scargle, J. D. 1982, ApJ, 263, 835 (URL: http://dx.doi.org/10.1086/160554, Bibcode: http://adsabs.harvard.edu/abs/1982ApJ...263..835S)\n[SS10] Sturrock, P. A., & Scargle, J. D. 2010, ApJ, 718, 527 (URL: http://dx.doi.org/10.1088/0004-637X/718/1/527, Bibcode: http://adsabs.harvard.edu/abs/2010ApJ...718..527S)The package provides facilities to:compute the periodogram using different methods (with different speeds) and different normalizations. This is one of the fastest implementations of these methods available as free software. If Julia is run with more than one thread, computation is automatically multi-threaded, further speeding up calculations;\naccess the frequency and period grid of the resulting periodogram, together with the power spectrum;\nfind the maximum power in the periodogram and the frequency and period corresponding to the peak. All these queries can be restricted to a specified region, in order to search a local maximum, instead of the global one;\ncalculate the probability that a peak arises from noise only (false-alarm probability) using analytic formulas, in order to assess the significance of the peak;\nperform bootstrap resamplings in order to compute the false-alarm probability with a statistical method;\ndetermine the best-fitting Lomb–Scargle model for the given data set at the given frequency."
},

{
    "location": "index.html#Installation-1",
    "page": "LombScargle.jl",
    "title": "Installation",
    "category": "section",
    "text": "LombScargle.jl is available for Julia 0.6 and later versions, and can be installed with Julia built-in package manager. In a Julia session run the commandsjulia> Pkg.update()\njulia> Pkg.add(\"LombScargle\")Older versions are also available for Julia 0.4 and 0.5."
},

{
    "location": "index.html#LombScargle.lombscargle-Tuple{AbstractArray{#s6,1} where #s6<:Real,Vararg{Any,N} where N}",
    "page": "LombScargle.jl",
    "title": "LombScargle.lombscargle",
    "category": "method",
    "text": "lombscargle(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},\n            errors::AbstractVector{<:Real}=ones(signal); keywords...)\n\nCompute the Lomb–Scargle periodogram of the signal vector, observed at times.  You can also specify the uncertainties for each signal point with errors argument.  All these vectors must have the same length.\n\nAll optional keywords are described in the docstring of LombScargle.plan.\n\nIf the signal has uncertainties, the signal vector can also be a vector of Measurement objects (from Measurements.jl package), in which case you don’t need to pass a separate errors vector for the uncertainties of the signal.\n\n\n\n"
},

{
    "location": "index.html#Usage-1",
    "page": "LombScargle.jl",
    "title": "Usage",
    "category": "section",
    "text": "After installing the package, you can start using it withusing LombScargleThe module defines a new LombScargle.Periodogram data type, which, however, is not exported because you will most probably not need to directly manipulate such objects. This data type holds both the frequency and the power vectors of the periodogram.The main function provided by the package is lombscargle:lombscargle(::AbstractVector{<:Real}, rest...)lombscargle returns a LombScargle.Periodogram. The only two mandatory arguments are:times: the vector of observation times\nsignal: the vector of observations associated with timesThe optional argument is:errors: the uncertainties associated to each signal point.All these vectors must have the same length.tip: Tip\nYou can pre-plan a periodogram with LombScargle.plan function, which has the same syntax as lombscargle described in this section. In this way the actual computation of the periodogram is faster and you will save memory. See the Planning the Periodogram section below.tip: Tip\nLombScargle.jl exploits Julia\'s native multi-threading for the non-fast methods (the methods used when you set the keyword fast=false). Run Julia with n threads (e.g., JULIA_NUM_THREADS=4 julia for 4 threads, if your machine has 4 physical cores) in order to automatically gain an n -fold scaling.Please note that multi-threading is still an experimental feature in Julia, so you may encounter issues when running it with more than one thread. For example, bug #17395 (if still open) may prevent the function, on some systems, from effectively scaling.If the signal has uncertainties, the signal vector can also be a vector of Measurement objects (from Measurements.jl package), in which case you need not to pass a separate errors vector for the uncertainties of the signal. You can create arrays of Measurement objects with the measurement function, see Measurements.jl manual at https://juliaphysics.github.io/Measurements.jl/stable for more details. The generalised Lomb–Scargle periodogram by [ZK09] is always used when the signal has uncertainties, because the original Lomb–Scargle algorithm cannot handle them.tip: Tip\nThe uncertainties are only used in the generalised Lomb–Scargle algorithm to build an inverse-variance weights vector (see [ZK09]), that gives more importance to datapoints with lower uncertainties. The case where all measurements have the same uncertainty (a condition known as homoskedasticity) results in a costant weights vector, like if there are no uncertainties at all. If you have homoskedastic errors, you do not need to provide them to lombscargle."
},

{
    "location": "index.html#LombScargle.plan",
    "page": "LombScargle.jl",
    "title": "LombScargle.plan",
    "category": "function",
    "text": "LombScargle.plan(times::AbstractVector{<:Real}, signal::AbstractVector{<:Real},\n                 errors::AbstractVector{<:Real}=ones(signal);\n                 normalization::Symbol=:standard,\n                 noise_level::Real=1,\n                 center_data::Bool=true,\n                 fit_mean::Bool=true,\n                 fast::Bool=true,\n                 flags::Integer=FFTW.ESTIMATE,\n                 timelimit::Real=Inf,\n                 oversampling::Integer=5,\n                 Mfft::Integer=4,\n                 samples_per_peak::Integer=5,\n                 nyquist_factor::Integer=5,\n                 minimum_frequency::Real=NaN,\n                 maximum_frequency::Real=NaN,\n                 frequencies::AbstractVector{Real}=\n                 autofrequency(times,\n                               samples_per_peak=samples_per_peak,\n                               nyquist_factor=nyquist_factor,\n                               minimum_frequency=minimum_frequency,\n                               maximum_frequency=maximum_frequency))\n\nPre-plan the Lomb–Scargle periodogram of the signal vector, observed at times.  The periodogram can then be computed by passing the result of this function to lombscargle.\n\nYou can also specify the uncertainties for each signal point with errors argument.  All these vectors must have the same length.\n\nOptional keywords arguments are:\n\nnormalization: how to normalize the periodogram.  Valid choices are: :standard, :model, :log, :psd, :Scargle, :HorneBaliunas, :Cumming\nnoise_level: the noise level used to normalize the periodogram when normalization is set to :Scargle\nfit_mean: if true, fit for the mean of the signal using the Generalised Lomb–Scargle algorithm (see Zechmeister & Kürster paper below).  If this is false and no uncertainty on the signal is provided, the original algorithm by Lomb and Scargle will be employed (see Townsend paper below)\ncenter_data: if true, subtract the weighted mean of signal from signal itself before performing the periodogram.  This is especially important if fit_mean is false\nfrequencies: the frequecy grid (not angular frequencies) at which the periodogram will be computed, as a vector.  If not provided, it is an evenly spaced grid of type Range, automatically determined with LombScargle.autofrequency function, which see.  See below for other available keywords that can be used to affect the frequency grid without directly setting frequencies\n\nYou can explicitely require to use or not the fast method by Press & Rybicki, overriding the default choice, by setting the fast keyword.  In any case, frequencies must be a Range object (this is the default) in order to actually use this method.  A few other keywords are available to adjust the settings of the periodogram when the fast method is used (otherwise they are ignored):\n\nfast: whether to use the fast method.\nflags: this integer keyword is a bitwise-or of FFTW planner flags, defaulting to FFTW.ESTIMATE.  Passing FFTW.MEASURE or FFTW.PATIENT will instead spend several seconds (or more) benchmarking different possible FFT algorithms and picking the fastest one; see the FFTW manual for more information on planner flags.\ntimelimit: specifies a rough upper bound on the allowed planning time, in seconds.\noversampling: oversampling the frequency factor for the approximation; roughly the number of time samples across the highest-frequency sinusoid. This parameter contains the tradeoff between accuracy and speed.\nMfft: the number of adjacent points to use in the FFT approximation.\n\nIn addition, you can use all optional keyword arguments of LombScargle.autofrequency function in order to tune the frequencies.\n\nIf the signal has uncertainties, the signal vector can also be a vector of Measurement objects (from Measurements.jl package), in which case you don’t need to pass a separate errors vector for the uncertainties of the signal.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.autofrequency",
    "page": "LombScargle.jl",
    "title": "LombScargle.autofrequency",
    "category": "function",
    "text": "autofrequency(times::AbstractVector{Real};\n              samples_per_peak::Integer=5,\n              nyquist_factor::Integer=5,\n              minimum_frequency::Real=NaN,\n              maximum_frequency::Real=NaN)\n\nDetermine a suitable frequency grid for the given vector of times.\n\nOptional keyword arguments are:\n\nsamples_per_peak: the approximate number of desired samples across the typical peak\nnyquist_factor: the multiple of the average Nyquist frequency used to choose the maximum frequency if maximum_frequency is not provided\nminimum_frequency: if specified, then use this minimum frequency rather than one chosen based on the size of the baseline\nmaximum_frequency: if specified, then use this maximum frequency rather than one chosen based on the average Nyquist frequency\n\nThis is based on prescription given at https://jakevdp.github.io/blog/2015/06/13/lomb-scargle-in-python/ and uses the same keywords names adopted in Astropy.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.lombscargle-Tuple{LombScargle.PeriodogramPlan}",
    "page": "LombScargle.jl",
    "title": "LombScargle.lombscargle",
    "category": "method",
    "text": "lombscargle(plan::PeriodogramPlan)\n\nCompute the Lomb–Scargle periodogram for the given plan.  This method has no other arguments.  See documentation of LombScargle.plan for how to plan a Lomb–Scargle periodogram.\n\n\n\n"
},

{
    "location": "index.html#Planning-the-Periodogram-1",
    "page": "LombScargle.jl",
    "title": "Planning the Periodogram",
    "category": "section",
    "text": "In a manner similar to planning Fourier transforms with FFTW, it is possible to speed-up computation of the Lomb–Scargle periodogram by pre-planning it with LombScargle.plan function. It has the same syntax as lombscargle, which in the base case is:LombScargle.plan\nLombScargle.autofrequencyLombScargle.plan takes all the same argument as lombscargle shown above and returns a LombScargle.PeriodogramPlan object after having pre-computed certain quantities needed afterwards, and pre-allocated the memory for the periodogram. It is highly suggested to plan a periodogram before actually computing it, especially for the fast method. Once you plan a periodogram, you can pass the LombScargle.PeriodogramPlan to lombscargle as the only argument.lombscargle(::LombScargle.PeriodogramPlan)Planning the periodogram has a twofold advantage. First of all, the planning stage is type-unstable, because the type of the plan depends on the value of input parameters, and not on their types. Thus, separating the planning (inherently inefficient) from the actual computation of the periodogram (completely type-stable) makes overall computation faster than directly calling lombscargle. Secondly, the LombScargle.PeriodogramPlan bears the time vector, but the quantities that are pre-computed in planning stage do not actually depend on it. This is particularly useful if you want to calculate the false-alarm probability via bootstrapping with LombScargle.bootstrap function: the vector time is randomly shuffled, but pre-computed quantities will remain the same, saving both time and memory in each iteration. In addition, you ensure that you will use the same options you used to compute the periodogram."
},

{
    "location": "index.html#Fast-Algorithm-1",
    "page": "LombScargle.jl",
    "title": "Fast Algorithm",
    "category": "section",
    "text": "When the frequency grid is evenly spaced, you can compute an approximate generalised Lomb–Scargle periodogram using a fast algorithm proposed by [PR89] that greatly speeds up calculations, as it scales as ON log(M) for N data points and M frequencies. For comparison, the true Lomb–Scargle periodogram has complexity ONM.  The larger the number of datapoints, the more accurate the approximation.note: Note\nThis method internally performs a Fast Fourier Transform (FFT) to compute some quantities, but it is in no way equivalento to conventional Fourier periodogram analysis.LombScargle.jl uses FFTW functions to compute the FFT. You can speed-up this task by using multi-threading: call FFTW.set_num_threads(n) to use n threads. However, please note that the running time will not scale as n because computation of the FFT is only a part of the algorithm.The only prerequisite in order to be able to employ this fast method is to provide a frequencies vector as a Range object, which ensures that the frequency grid is perfectly evenly spaced. This is the default, since LombScargle.autofrequency returns a Range object.tip: Tip\nIn Julia, a Range object can be constructed for example with the linspace function (you specify the start and the end of the range, and optionally the length of the vector) or with the syntax start:stop (you specify the start and the end of the range, and optionally the linear step; a related function is colon). Somewhere in the middle is the range function: you specify the start of the range and the length of the vector, and optionally the linear step.Since this fast method is accurate only for large datasets, it is enabled by default only if the number of output frequencies is larger than 200. You can override the default choice of using this method by setting the fast keyword to true or false. We recall that in any case, the frequencies vector must be a Range in order to use this method.To summarize, provided that frequencies vector is a Range object, you can use the fast method:by default if the length of the output frequency grid is larger than 200 points\nin any case with the fast=true keywordSetting fast=false always ensures you that this method will not be used, instead fast=true actually enables it only if frequencies is a Range."
},

{
    "location": "index.html#Normalization-1",
    "page": "LombScargle.jl",
    "title": "Normalization",
    "category": "section",
    "text": "By default, the periodogram p(f) is normalized so that it has values in the range 0 leq p(f) leq 1, with p = 0 indicating no improvement of the fit and p = 1 a \"perfect\" fit (100% reduction of chi^2 or chi^2 = 0). This is the normalization suggested by [LOM76] and [ZK09], and corresponds to the :standard normalization in lombscargle function. [ZK09] wrote the formula for the power of the periodogram at frequency f asp(f) = frac1YYleftfracYC^2_tauCC_tau +\nfracYS^2_tauSS_taurightSee the paper for details. The other normalizations for periodograms P(f) are calculated from this one. In what follows, N is the number of observations.:model:\nP(f) = fracp(f)1 - p(f)\n:log:\nP(f) = -log(1 - p(f))\n:psd:\nP(f) = frac12leftfracYC^2_tauCC_tau +\nfracYS^2_tauSS_tauright = p(f) fracYY2\n:Scargle:\nP(f) = fracp(f)textnoise level\nThis normalization can be used when you know the noise level (expected from the a priori known noise variance or population variance), but this isn\'t usually the case. See [SCA82]\n:HorneBaliunas:\nP(f) = fracN - 12 p(f)\nThis is like the :Scargle normalization, where the noise has been estimated for Gaussian noise to be (N - 1)2. See [HB86]\nIf the data contains a signal or if errors are under- or overestimated or if intrinsic variability is present, then (N-1)2 may not be a good uncorrelated estimator for the noise level.  [CMB99] suggested to estimate the noise level a posteriori with the residuals of the best fit and normalised the periodogram as:\nP(f) = fracN - 32 fracp(f)1 - p(f_textbest)\nThis is the :Cumming normalization option"
},

{
    "location": "index.html#LombScargle.power",
    "page": "LombScargle.jl",
    "title": "LombScargle.power",
    "category": "function",
    "text": "power(p::Periodogram)\n\nReturn the power vector of Lomb–Scargle periodogram p.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.freq",
    "page": "LombScargle.jl",
    "title": "LombScargle.freq",
    "category": "function",
    "text": "freq(p::Periodogram)\n\nReturn the frequency vector of Lomb–Scargle periodogram p.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.freqpower",
    "page": "LombScargle.jl",
    "title": "LombScargle.freqpower",
    "category": "function",
    "text": "freqpower(p::Periodogram)\n\nReturn the 2-tuple (freq(p), power(p)), where freq(p) and power(p) are the frequency vector and the power vector of Lomb–Scargle periodogram p respectively.\n\n\n\n"
},

{
    "location": "index.html#Access-Frequency-Grid-and-Power-Spectrum-of-the-Periodogram-1",
    "page": "LombScargle.jl",
    "title": "Access Frequency Grid and Power Spectrum of the Periodogram",
    "category": "section",
    "text": "power\nfreq\nfreqpowerlombscargle function returns a LombScargle.Periodogram object, but you most probably want to use the frequency grid and the power spectrum. You can access these vectors with freq and power functions, just like in DSP.jl package. If you want to get the 2-tuple (freq(p), power(p)) use the freqpower function."
},

{
    "location": "index.html#LombScargle.period",
    "page": "LombScargle.jl",
    "title": "LombScargle.period",
    "category": "function",
    "text": "power(p::Periodogram)\n\nReturn the period vector of Lomb–Scargle periodogram p.  It is equal to 1 ./ freq(p).\n\n\n\n"
},

{
    "location": "index.html#LombScargle.periodpower",
    "page": "LombScargle.jl",
    "title": "LombScargle.periodpower",
    "category": "function",
    "text": "periodpower(p::Periodogram)\n\nReturn the 2-tuple (period(p), power(p)), where period(p) and power(p) are the period vector and the power vector of Lomb–Scargle periodogram p respectively.\n\n\n\n"
},

{
    "location": "index.html#Access-Periods-and-their-and-Power-in-the-Periodogram-1",
    "page": "LombScargle.jl",
    "title": "Access Periods and their and Power in the Periodogram",
    "category": "section",
    "text": "period\nperiodpowerThese utilities are the analogs of freq and freqpower, but relative to the periods instead of the frequencies. Thus period(p) returns the vector of periods in the periodogram, that is 1./freq(p), and periodpower(p) gives you the 2-tuple (period(p), power(p))."
},

{
    "location": "index.html#LombScargle.findmaxpower",
    "page": "LombScargle.jl",
    "title": "LombScargle.findmaxpower",
    "category": "function",
    "text": "findmaxpower(p::Periodogram)\n\nReturn the highest power of the periodogram p.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.findmaxfreq",
    "page": "LombScargle.jl",
    "title": "LombScargle.findmaxfreq",
    "category": "function",
    "text": "findmaxfreq(p::Periodogram, [interval::AbstractVector{Real}], threshold::Real=findmaxpower(p))\n\nReturn the array of frequencies with the highest power in the periodogram p. If a scalar real argument threshold is provided, return the frequencies with power larger than or equal to threshold.  If you want to limit the search to a narrower frequency range, pass as second argument a vector with the extrema of the interval.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.findmaxperiod",
    "page": "LombScargle.jl",
    "title": "LombScargle.findmaxperiod",
    "category": "function",
    "text": "findmaxperiod(p::Periodogram, [interval::AbstractVector{Real}], threshold::Real=findmaxpower(p))\n\nReturn the array of period with the highest power in the periodogram p.  If a scalar real argument threshold is provided, return the period with power larger than or equal to threshold.  If you want to limit the search to a narrower period range, pass as second argument a vector with the extrema of the interval.\n\n\n\n"
},

{
    "location": "index.html#findmaxpower,-findmaxfreq,-and-findmaxperiod-Functions-1",
    "page": "LombScargle.jl",
    "title": "findmaxpower, findmaxfreq, and findmaxperiod Functions",
    "category": "section",
    "text": "findmaxpower\nfindmaxfreq\nfindmaxperiodOnce you compute the periodogram, you usually want to know which are the frequencies or periods with highest power. To do this, you can use the findmaxfreq and findmaxperiod functions.  They return the vector of frequencies and periods, respectively, with the highest power in the periodogram p. If a scalar real argument threshold is provided, return the frequencies with power larger than or equal to threshold. If you want to limit the search to a narrower frequency or period range, pass as second argument a vector with the extrema of the interval.The value of the highest power of a periodogram can be calculated with the findmaxpower function."
},

{
    "location": "index.html#LombScargle.prob-Tuple{LombScargle.Periodogram,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.prob",
    "category": "method",
    "text": "prob(P::Periodogram, pow::Real)\n\nReturn the probability that the periodogram power can exceed the value pow.\n\nIts inverse is the probinv function.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.probinv-Tuple{LombScargle.Periodogram,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.probinv",
    "category": "method",
    "text": "probinv(P::Periodogram, prob::Real)\n\nReturn the power value of the periodogram power whose probability is prob.\n\nThis is the inverse of prob function.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.M",
    "page": "LombScargle.jl",
    "title": "LombScargle.M",
    "category": "function",
    "text": "LombScargle.M(P::Periodogram)\n\nEstimates the number of independent frequencies in the periodogram P.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.fap-Tuple{LombScargle.Periodogram,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.fap",
    "category": "method",
    "text": "fap(P::Periodogram, pow::Real)\n\nReturn the false-alarm probability for periodogram P and power value pow.\n\nIts inverse is the fapinv function.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.fapinv-Tuple{LombScargle.Periodogram,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.fapinv",
    "category": "method",
    "text": "fapinv(P::Periodogram, prob::Real)\n\nReturn the power value of the periodogram whose false-alarm probability is prob.\n\nThis is the inverse of fap function.\n\n\n\n"
},

{
    "location": "index.html#False-Alarm-Probability-1",
    "page": "LombScargle.jl",
    "title": "False-Alarm Probability",
    "category": "section",
    "text": "prob(::LombScargle.Periodogram, ::Real)\nprobinv(::LombScargle.Periodogram, ::Real)\nLombScargle.M\nfap(::LombScargle.Periodogram, ::Real)\nfapinv(::LombScargle.Periodogram, ::Real)Noise in the data produce fluctuations in the periodogram that will present several local peaks, but not all of them related to real periodicities. The significance of the peaks can be tested by calculating the probability that its power can arise purely from noise.  The higher the value of the power, the lower will be this probability.note: Note\n[CMB99] showed that the different normalizations result in different probability functions. LombScargle.jl can calculate the probability (and the false-alarm probability) only for the normalizations reported by [ZK09], that are :standard, :Scargle, :HorneBaliunas, and :Cumming.The probability textProb(p  p_0) that the periodogram power p can exceed the value p_0 can be calculated with the prob function, whose first argument is the periodogram and the second one is the p_0 value. The function probinv is its inverse: it takes the probability as second argument and returns the corresponding p_0 value.Here are the probability functions for each normalization supported by LombScargle.jl::standard (p in 0 1):\ntextProb(p  p_0) = (1 - p_0)^(N - 3)2\n:Scargle (p in 0 infty)):\ntextProb(p  p_0) = exp(-p_0)\n:HorneBaliunas (p in 0 (N - 1)2):\ntextProb(p  p_0) = left(1 - frac2p_0N - 1right)^(N - 3)2\n:Cumming (p in 0 infty)):\ntextProb(p  p_0) = left(1 + frac2p_0N - 3right)^-(N - 3)2As explained by [SS10], «the term \"false-alarm probability denotes the probability that at least one out of M independent power values in a prescribed search band of a power spectrum computed from a white-noise time series is expected to be as large as or larger than a given value». LombScargle.jl provides the fap function to calculate the false-alarm probability (FAP) of a given power in a periodogram. Its first argument is the periodogram, the second one is the value p_0 of the power of which you want to calculate the FAP. The function fap uses the formulatextFAP = 1 - (1 - textProb(p  p_0))^Mwhere M is the number of independent frequencies estimated with M = T cdot Delta f, being T the duration of the observations and Delta f the width of the frequency range in which the periodogram has been calculated (see [CUM04]). The function fapinv is the inverse of fap: it takes as second argument the value of the FAP and returns the corresponding value p_0 of the power.The detection threshold p_0 is the periodogram power corresponding to some (small) value of textFAP, i.e. the value of p exceeded due to noise alone in only a small fraction textFAP of trials. An observed power larger than p_0 indicates that a signal is likely present (see [CUM04]).warning: Warning\nSome authors stressed that this method to calculate the false-alarm probability is not completely reliable. A different approach to calculate the false-alarm probability is to perform Monte Carlo or bootstrap simulations in order to determine how often a certain power level p_0 is exceeded just by chance (see [CMB99], [CUM04], and [ZK09]). See next section."
},

{
    "location": "index.html#LombScargle.bootstrap",
    "page": "LombScargle.jl",
    "title": "LombScargle.bootstrap",
    "category": "function",
    "text": "LombScargle.bootstrap(N::Integer,\n                      times::AbstractVector{Real},\n                      signal::AbstractVector{Real},\n                      errors::AbstractVector{Real}=ones(signal); ...)\n\nCreate N bootstrap samples, perform the Lomb–Scargle analysis on them, and store all the highest peaks for each one in a LombScargle.Bootstrap object. All the arguments after N are passed around to lombscargle, which see.\n\n\n\nLombScargle.bootstrap(N::Integer, plan::PeriodogramPlan)\n\nCreate N bootstrap samples, perform the Lomb–Scargle analysis on them for the given plan, and store all the highest peaks for each one in a LombScargle.Bootstrap object.\n\nSee documentation of LombScargle.plan for how to plan a Lomb–Scargle periodogram.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.fap-Tuple{LombScargle.Bootstrap{#s6} where #s6<:AbstractFloat,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.fap",
    "category": "method",
    "text": "fap(b::Bootstrap, power::Real)\n\nReturn the false-alarm probability for power in the bootstrap sample b.\n\nIts inverse is the fapinv function.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.fapinv-Tuple{LombScargle.Bootstrap{#s6} where #s6<:AbstractFloat,Real}",
    "page": "LombScargle.jl",
    "title": "LombScargle.fapinv",
    "category": "method",
    "text": "fapinv(b::Bootstrap, prob::Real)\n\nReturn the power value whose false-alarm probability is prob in the bootstrap sample b.\n\nIt returns NaN if the requested probability is too low and the power cannot be determined with the bootstrap sample b.  In this case, you should enlarge your bootstrap sample so that N*fap can be rounded to an integer larger than or equal to 1.\n\nThis is the inverse of fap function.\n\n\n\n"
},

{
    "location": "index.html#Bootstrapping-1",
    "page": "LombScargle.jl",
    "title": "Bootstrapping",
    "category": "section",
    "text": "LombScargle.bootstrap\nfap(::LombScargle.Bootstrap{<:AbstractFloat}, ::Real)\nfapinv(::LombScargle.Bootstrap{<:AbstractFloat}, ::Real)One of the possible and most simple statistical methods that you can use to measure the false-alarm probability and its inverse is bootstrapping (see section 4.2.2 of [MHC93]).note: Note\nWe emphasize that you can use this method only if you know your data points are independent and identically distributed, and they have white uncorrelated noise.The recipe of the bootstrap method is very simple to implement:repeat the Lomb–Scargle analysis a large number N of times on the original data, but with the signal (and errors, if present) vector randomly shuffled. As an alternative, shuffle only the time vector;\nout of all these simulations, store the powers of the highest peaks;\nin order to estimate the false-alarm probability of a given power, count how many times the highest peak of the simulations exceeds that power, as a fraction of N. If you instead want to find the inverse of the false-alarm probability textprob, looks for the Ncdottextprob-th element of the highest peaks vector sorted in descending order.Remember to pass to lombscargle function the same options, if any, you used to compute the Lomb–Scargle periodogram before.LombScargle.jl provides simple methods to perform such analysis. The LombScargle.bootstrap function allows you to create a bootstrap sample with N permutations of the original data. All the arguments after the first one are passed around to lombscargle. The output is a LombScargle.Bootstrap object.You can also pass to LombScargle.bootstrap a pre-computed LombScargle.PeriodogramPlan as second argument (this method takes no other argument nor keyword). In this way you will be sure to use exactly the same options you used before for computing the periodogram with the same periodogram plan.The false-alarm probability and its inverse can be calculated with fap and fapinv functions respectively.  Their syntax is the same as the methods introduced above, but with a LombScargle.Bootstrap object as first argument, instead of the LombScargle.Periodogram one."
},

{
    "location": "index.html#LombScargle.model",
    "page": "LombScargle.jl",
    "title": "LombScargle.model",
    "category": "function",
    "text": "LombScargle.model(times::AbstractVector{Real},\n                  signal::AbstractVector{R2},\n                  [errors::AbstractVector{R3},]\n                  frequency::Real,\n                  [times_fit::AbstractVector{R4}];\n                  center_data::Bool=true,\n                  fit_mean::Bool=true)\n\nReturn the best fitting Lomb–Scargle model for the given signal at the given frequency.\n\nMandatory arguments are:\n\ntimes: the observation times\nsignal: the signal, sampled at times (must have the same length as times)\nfrequency: the frequency at which to calculate the model\n\nOptional arguments are:\n\nerrors: the vector of uncertainties of the signal.  If provided, it must have the same length as signal and times, and be the third argument\ntimes_fit: the vector of times at which the model will be calculated.  It defaults to times.  If provided, it must come after frequency\n\nOptional keyword arguments center_data and fit_mean have the same meaning as in lombscargle function, which see.\n\n\n\n"
},

{
    "location": "index.html#LombScargle.model-Function-1",
    "page": "LombScargle.jl",
    "title": "LombScargle.model Function",
    "category": "section",
    "text": "LombScargle.modelFor each frequency f (and hence for the corresponding angular frequency omega = 2pi f) the Lomb–Scargle algorithm looks for the sinusoidal function of the typea_fcos(omega t) + b_fsin(omega t) + c_fthat best fits the data. In the original Lomb–Scargle algorithm the offset c is null (see [LOM76]). In order to find the best-fitting coefficients a_f, b_f, and c_f for the given frequency f, without actually performing the periodogram, you can solve the linear system mathbfAx = mathbfy, where mathbfA is the matrixbeginaligned\nbeginbmatrix\n  cos(omega t)  sin(omega t)  1\nendbmatrix =\nbeginbmatrix\n  cos(omega t_1)  sin(omega t_1)  1      \n  vdots              vdots              vdots \n  cos(omega t_n)  sin(omega t_n)  1\nendbmatrix\nendalignedt = t_1 dots t_n^textT is the column vector of observation times, x is the column vector with the unknown coefficientsbeginaligned\nbeginbmatrix\n  a_f \n  b_f \n  c_f\nendbmatrix\nendalignedand textbfy is the column vector of the signal. The solution of the matrix gives the wanted coefficients.This is what the LombScargle.model function does in order to return the best fitting Lomb–Scargle model for the given signal at the given frequency.Mandatory arguments are:times: the observation times\nsignal: the signal, sampled at times (must have the same length as times)\nfrequency: the frequency at which to calculate the modelThe optional arguments are:errors: the vector of uncertainties of the signal. If provided, it must have the same length as signal and times, and be the third argument. Like for lombscargle, if the signal has uncertainties, the signal vector can also be a vector of Measurement objects, and this argument should be omitted\ntimes_fit: the vector of times at which the model will be calculated. It defaults to times. If provided, it must come after frequencyOptional boolean keywords center_data and fit_mean have the same meaning as in lombscargle function:fit_mean: whether to fit for the mean. If this is false, like in the original Lomb–Scargle periodogram, mathbfA does not have the third column of ones, c_f is set to 0 and the unknown vector to be determined becomes x = a_f b_f^textT\ncenter_data: whether the data should be pre-centered before solving the linear system. This is particularly important if fit_mean=false"
},

{
    "location": "index.html#Examples-1",
    "page": "LombScargle.jl",
    "title": "Examples",
    "category": "section",
    "text": "Here is an example of a noisy periodic signal (sin(pi t) + 15cos(2pi t)) sampled at unevenly spaced times.julia> using LombScargle\n\njulia> ntimes = 1001\n1001\n\njulia> t = linspace(0.01, 10pi, ntimes) # Observation times\n0.01:0.03140592653589793:31.41592653589793\n\njulia> t += step(t)*rand(ntimes) # Randomize times\n\njulia> s = sinpi.(t) .+ 1.5cospi.(2t) .+ rand(ntimes) # The signal\n\njulia> plan = LombScargle.plan(t, s); # Pre-plan the periodogram\n\njulia> pgram = lombscargle(plan) # Compute the periodogram\nLombScargle.Periodogram{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Array{Float64,1}}([0.000472346, 0.000461633, 0.000440906, 0.000412717, 0.000383552, 0.000355828, 0.000289723, 0.000154585, 3.44734e-5, 5.94437e-7  …  3.15125e-5, 0.000487391, 0.0018939, 0.00367003, 0.00484181, 0.00495189, 0.00453233, 0.00480968, 0.00619657, 0.0074052], 0.003185690706734265:0.00637138141346853:79.72190993602499, [0.0295785, 0.0540516, 0.0780093, 0.122759, 0.15685, 0.192366, 0.206601, 0.252829, 0.265771, 0.315443  …  31.1512, 31.1758, 31.2195, 31.2342, 31.2752, 31.293, 31.3517, 31.3761, 31.4148, 31.4199], :standard)You can plot the result, for example with Plots package. Use freqpower function to get the frequency grid and the power of the periodogram as a 2-tuple.using Plots\nplot(freqpower(pgram)...)(Image: image)You can also plot the power vs the period, instead of the frequency, with periodpower:using Plots\nplot(periodpower(pgram)...)(Image: image)warning: Warning\nIf you do not fit for the mean of the signal (fit_mean=false keyword to lombscargle function) without centering the data (center_data=false) you can get inaccurate results. For example, spurious peaks at low frequencies can appear and the real peaks lose power:plot(freqpower(lombscargle(t, s, fit_mean=false, center_data=false))...)(Image: image)tip: Tip\nYou can tune the frequency grid with appropriate keywords to lombscargle function. For example, in order to increase the sampling increase samples_per_peak, and set maximum_frequency to lower values in order to narrow the frequency range:plot(freqpower(lombscargle(t, s, samples_per_peak=20, maximum_frequency=1.5))...)(Image: image)If you simply want to use your own frequency grid, directly set the frequencies keyword:plot(freqpower(lombscargle(t, s, frequencies=0.001:1e-3:1.5))...)(Image: image)"
},

{
    "location": "index.html#Signal-with-Uncertainties-1",
    "page": "LombScargle.jl",
    "title": "Signal with Uncertainties",
    "category": "section",
    "text": "The generalised Lomb–Scargle periodogram is able to handle a signal with uncertainties, and they will be used as weights in the algorithm.  The uncertainties can be passed either as the third optional argument errors to lombscargle or by providing this function with a signal vector of type Measurement (from Measurements.jl package).using Measurements, Plots\nntimes = 1001\nt = linspace(0.01, 10pi, ntimes)\ns = sinpi.(2t)\nerrors = rand(0.1:1e-3:4.0, ntimes)\n# Run one of the two following equivalent commands\nplot(freqpower(lombscargle(t, s, errors, maximum_frequency=1.5))...)\nplot(freqpower(lombscargle(t, measurement(s, errors), maximum_frequency=1.5))...)(Image: image)This is the plot of the power versus the period:# Run one of the two following equivalent commands\nplot(periodpower(lombscargle(t, s, errors, maximum_frequency=1.5))...)\nplot(periodpower(lombscargle(t, measurement(s, errors), maximum_frequency=1.5))...)(Image: image)We recall that the generalised Lomb–Scargle algorithm is used when the fit_mean optional keyword to lombscargle is true if no error is provided, instead it is always used if the signal has uncertainties."
},

{
    "location": "index.html#Find-Highest-Power-and-Associated-Frequencies-and-Periods-1",
    "page": "LombScargle.jl",
    "title": "Find Highest Power and Associated Frequencies and Periods",
    "category": "section",
    "text": "findmaxfreq function tells you the frequencies with the highest power in the periodogram (and you can get the period by taking its inverse):julia> t = linspace(0, 10, 1001);\n\njulia> s = sinpi.(t);\n\njulia> plan = LombScargle.plan(t, s); # Plan the periodogram\n\njulia> p = lombscargle(plan);\n\njulia> findmaxperiod(p) # Period with highest power\n1-element Array{Float64,1}:\n 0.00498778\n\njulia> findmaxfreq(p) # Frequency with the highest power\n1-element Array{Float64,1}:\n 200.49This peak is at high frequencies, very far from the expected value of the period of 2. In order to find the real peak, you can either narrow the ranges in order to exclude higher armonicsjulia> findmaxperiod(p, [1, 10]) # Limit the search to periods in [1, 10]\n1-element Array{Float64,1}:\n 2.04082\n\njulia> findmaxfreq(p, [0.1, 1]) # Limit the search to frequencies in [0.1, 1]\n1-element Array{Float64,1}:\n 0.49or pass the threshold argument to findmaxfreq or findmaxperiod. You can use findmaxpower to discover the highest power in the periodogram:julia> findmaxpower(p)\n0.9958310178312316\n\njulia> findmaxperiod(p, 0.95)\n10-element Array{Float64,1}:\n 2.04082\n 1.96078\n 0.0100513\n 0.0100492\n 0.00995124\n 0.00994926\n 0.00501278\n 0.00501228\n 0.00498778\n 0.00498728\n\njulia> findmaxfreq(p, 0.95)\n10-element Array{Float64,1}:\n   0.49\n   0.51\n  99.49\n  99.51\n 100.49\n 100.51\n 199.49\n 199.51\n 200.49\n 200.51The first peak is the real one, the other double peaks appear at higher armonics.tip: Tip\nUsually, plotting the periodogram can give you a clue of what\'s going on."
},

{
    "location": "index.html#Significance-of-the-Peaks-1",
    "page": "LombScargle.jl",
    "title": "Significance of the Peaks",
    "category": "section",
    "text": "The significance of the peaks in the Lomb–Scargle periodogram can be assessed by measuring the False-Alarm Probability. Analytic expressions of this quantity and its inverse can be obtained with the fap and fapinv functions, respectively.julia> t = linspace(0.01, 20, samples_per_peak = 10)\n\njulia> s = sinpi.(e.*t).^2 .- cos.(5t).^4\n\njulia> plan = LombScargle.plan(t, s);\n\njulia> p = lombscargle(plan)\n\n# Find the false-alarm probability for the highest peak.\njulia> fap(p, 0.3)\n0.028198095962262748Thus, a peak with power 03 has a probability of 0028 that it is due to noise only. A quantity that is often used is the inverse of the false-alarm probability as well: what is the minimum power whose false-alarm probability is lower than the given probability? For example, if you want to know the minimum power for which the false-alarm probability is at most 001 you can use:julia> fapinv(p, 0.01)\n0.3304696923786712As we already noted, analytic expressions of the false-alarm probability and its inverse may not be reliable if your data does not satisfy specific assumptions. A better way to calculate this quantity is to use statistical methods. One of this is bootstrapping. In LombScargle.jl, you can use the function LombScargle.bootstrap to create a bootstrap sample and then you can calculate the false-alarm probability and its inverse using this sample.tip: Tip\nWhen applying the bootstrap method you should use the same options you used to perform the periodogram on your data. Using the same periodogram plan you used to compute the periodogram will ensure that you use the same options. However, note that the fast method gives approximate results that for some frequencies may not be reliable (they can go outside the range 0 1 for the standard normalization). More robust results can be obtained with the fast = false option.# Create a bootstrap sample with 10000\n# resamplings of the original data, re-using the\n# same periodogram plan.  The larger the better.\n# This may take some minutes.\njulia> b = LombScargle.bootstrap(10000, plan)\n\n# Calculate the false-alarm probability of a peak\n# with power 0.3 using this bootstrap sample.\njulia> fap(b, 0.3)\n0.0209\n\n# Calculate the lowest power that has probability\n# less than 0.01 in this bootstrap sample.\njulia> fapinv(b, 0.01)\n0.3268290388848437If you query fapinv with a too low probability, the corresponding power cannot be determined and you will get NaN as result.julia> fapinv(b, 1e-5)\nNaNIf you want to find the power corresponding to a false-alarm probability of textprob = 10^-5, you have to create a new bootstrap sample with N resamplings so that Ncdottextprob can be rounded to an integer larger than or equal to one (for example N = 10^5)."
},

{
    "location": "index.html#Find-the-Best-Fitting-Model-1",
    "page": "LombScargle.jl",
    "title": "Find the Best-Fitting Model",
    "category": "section",
    "text": "The LombScargle.model function can help you to test whether a certain frequency fits well your data.using Plots\nt = linspace(0.01, 10pi, 1000) # Observation times\ns = sinpi.(t) .+ 1.2cospi.(t) .+ 0.3rand(length(t)) # The noisy signal\n# Pick-up the best frequency\nf = findmaxfreq(lombscargle(t, s, maximum_frequency=10, samples_per_peak=20))[1]\nt_fit = linspace(0, 1)\ns_fit = LombScargle.model(t, s, f, t_fit/f) # Determine the model\nscatter(mod.(t.*f, 1), s, lab=\"Phased data\", title=\"Best Lomb-Scargle frequency: $f\")\nplot!(t_fit, s_fit, lab=\"Best-fitting model\", linewidth=4)(Image: image)tip: Tip\nIf there are more than one dominant frequency you may need to consider more models. This task may require some work and patience. Plot the periodogram in order to find the best frequencies.using Plots\nt = linspace(0.01, 5, 1000) # Observation times\ns = sinpi.(2t) .+ 1.2cospi.(4t) .+ 0.3rand(length(t)) # Noisy signal\nplan = LombScargle.plan(t, s, samples_per_peak=50)\np = lombscargle(plan)\n# After plotting the periodogram, you discover\n# that it has two prominent peaks around 1 and 2.\nf1 = findmaxfreq(p, [0.8, 1.2])[1] # Get peak frequency around 1\nf2 = findmaxfreq(p, [1.8, 2.2])[1] # Get peak frequency around 2\nfit1 = LombScargle.model(t, s, f1) # Determine the first model\nfit2 = LombScargle.model(t, s, f2) # Determine the second model\nscatter(t, s, lab=\"Data\", title=\"Best-fitting Lomb-Scargle model\")\nplot!(t, fit1 + fit2, lab=\"Best-fitting model\", linewidth=4)(Image: image)"
},

{
    "location": "index.html#Performance-1",
    "page": "LombScargle.jl",
    "title": "Performance",
    "category": "section",
    "text": "A pre-planned periodogram in LombScargle.jl computed in single thread mode with the fast method is more than 2.9 times faster than the implementation of the same algorithm provided by AstroPy, and more than 4.5 times faster if 4 FFTW threads are used (on machines with at least 4 physical CPUs).The following plot shows a comparison between the times needed to compute a periodogram for a signal with N datapoints using LombScargle.jl, with 1 or 4 threads (with flags = FFTW.MEASURE for better performance), and the single-threaded AstroPy implementation.  (Julia version: 0.7.0-DEV.2309, commit 7ae9955c93; LombScargle.jl version: 0.3.1; Python version: 3.5.4; Astropy version: 2.0.2. CPU: Intel(R) Core(TM) i7-4700MQ.)(Image: image)Note that this comparison is unfair, as AstroPy doesn\'t support pre-planning a periodogram nor exploiting multi-threading. A non-planned periodogram in single thread mode in LombScargle.jl is still twice faster than AstroPy."
},

{
    "location": "index.html#Development-1",
    "page": "LombScargle.jl",
    "title": "Development",
    "category": "section",
    "text": "The package is developed at https://github.com/giordano/LombScargle.jl. There you can submit bug reports, make suggestions, and propose pull requests."
},

{
    "location": "index.html#History-1",
    "page": "LombScargle.jl",
    "title": "History",
    "category": "section",
    "text": "The ChangeLog of the package is available in NEWS.md file in top directory."
},

{
    "location": "index.html#License-1",
    "page": "LombScargle.jl",
    "title": "License",
    "category": "section",
    "text": "The LombScargle.jl package is licensed under the BSD 3-clause \"New\" or \"Revised\" License. The original author is Mosè Giordano."
},

{
    "location": "index.html#Acknowledgements-1",
    "page": "LombScargle.jl",
    "title": "Acknowledgements",
    "category": "section",
    "text": "This package adapts the implementation in Astropy of the the fast Lomb–Scargle method by [PR89]. We claim no endorsement nor promotion by the Astropy Team."
},

]}
