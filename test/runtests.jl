using LombScargle
using Measurements
using Base.Test

Threads.nthreads() > 1 && info("Julia running in multi-threading mode")

# Fix random seed to make test reproducible.
srand(1)

ntimes = 1001
# Observation times
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += step(t)*rand(ntimes)
# The signal
s = sinpi.(t) + cospi.(2t) + rand(ntimes)
# Frequency grid
nfreqs = 10000
freqs = linspace(0.01, 3, nfreqs)
# Randomize frequency grid
freqs += step(freqs)*rand(nfreqs)
# Use "freqpower" and "periodpower" just to call that function and increase code
# coverage.  "autofrequency" function is tested below.
pgram1 = lombscargle(t, s, frequencies=freqs, fit_mean=false)
pgram2 = lombscargle(t, s, frequencies=freqs, fit_mean=true)
@test freqpower(pgram1)[2] ≈ periodpower(pgram2)[2] atol = 5e-3
# Make sure there are no infinities in `pgram1'.  It seems to work only on
# 64-bit systems.
Sys.WORD_SIZE == 64 && @test(!(Inf in power(pgram1)) && !(Inf in power(pgram2)))

# Simple signal, without any randomness
t = collect(linspace(0.01, 10pi, ntimes))
s = sin.(t)
pgram1 = lombscargle(t, s, fast = false, fit_mean=false)
pgram2 = lombscargle(t, s, fast = false, fit_mean=true)
@test power(pgram1) ≈ power(pgram2) atol = 4e-7
pgram3 = lombscargle(t, s, fast = false, center_data=false, fit_mean=false)
pgram4 = lombscargle(t, s, fast = false, center_data=false, fit_mean=true)
@test power(pgram3) ≈ power(pgram4) atol = 4e-7

# Test the values in order to prevent wrong results in both algorithms
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=true)) ≈ [0.029886871262324886,0.0005456198989410226,1.912507742056023e-5, 4.54258409531214e-6, 1.0238342782430832e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=false)) ≈ [0.02988686776042212, 0.0005456197937194695,1.9125076826683257e-5,4.542583863304549e-6,1.0238340733199874e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, center_data=false, fit_mean=true)) ≈ [0.029886871262325004,0.0005456198989536703,1.9125077421448458e-5,4.5425840956285145e-6,1.023834278337881e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, center_data=false, fit_mean=false)) ≈ [0.029886868328967767,0.0005456198924872134,1.9125084251687147e-5,4.542588504467314e-6,1.0238354525870936e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:model)) ≈ [0.030807614469885718,0.0005459177625354441,1.9125443196143085e-5,4.54260473047638e-6,1.0238447607164715e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:log)) ≈ [0.030342586720560734,0.0005457688036440774,1.9125260307148152e-5,4.542594412890309e-6,1.0238395194654036e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:psd)) ≈ [7.474096700871138,0.1364484040771917,0.004782791641128195,0.0011360075968541799,0.002560400630125523]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:Scargle)) ≈ [0.029886871262324904,0.0005456198989410194,1.912507742056126e-5,4.54258409531238e-6,1.0238342782428552e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:HorneBaliunas)) ≈ [14.943435631162451,0.2728099494705097,0.009562538710280628,0.00227129204765619,0.005119171391214276]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:Cumming)) ≈ [15.372999620472974,0.2806521440709115,0.009837423440787873,0.0023365826071340815,0.005266327088140394]
@test_throws ErrorException lombscargle(t, s, frequencies=0.2:0.2:1, normalization=:foo)

# Test signal with uncertainties
err = collect(linspace(0.5, 1.5, ntimes))
@test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=true)) ≈ [0.06659683848818691,0.09230959166317589,0.006625919314669043,0.0015664010997692612,0.0005085442118408477,0.00019704659245878378,9.658452525613897e-5,6.331573873913433e-5,4.903871967643573e-5,3.7948448825374025e-5]
@test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=false)) ≈ [0.0664002483305464,0.09219168665786254,0.006625915010614472,0.0015663421089042564,0.0005085109569237008,0.00019703233981948823,9.6577091433651e-5,6.33101344670203e-5,4.9033581990442793e-5,3.7944076990210425e-5]
@test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=false, center_data=false)) ≈ [0.06920814049261209,0.09360344864985352,0.006634919960009565,0.0015362072871144769,0.0004858250831632676,0.00018179850370583626,8.543727416919218e-5,5.379994730581837e-5,4.0107232867763e-5,2.9784059487535237e-5]
@test power(lombscargle(t, s, err)) ==
    power(lombscargle(t, measurement.(s, err)))

# Test fast method
t = collect(linspace(0.01, 10pi, ntimes))
s = sin.(t)
pgram5 = lombscargle(t, s, maximum_frequency=30, fast=true)
pgram6 = lombscargle(t, s, maximum_frequency=30, fast=false)
@test power(pgram5) ≈ power(pgram6) atol = 3e-6
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, fit_mean=true)) ≈ [0.029886871262324963,0.0005453105325981758,1.9499330722168046e-5,2.0859593514888897e-6,1.0129019249708592e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, fit_mean=false)) ≈ [0.029886867760422008,0.0005453104197620392,1.9499329579010375e-5,2.085948496002562e-6,1.0128073536975395e-5]
@test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, center_data=false, fit_mean=false)) ≈ [0.029886868328967718,0.0005453105220405559,1.949931928224576e-5,2.0859802347505357e-6,1.0127777365273726e-5]
@test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, fit_mean=true)) ≈ [0.09230959166317655,0.0015654929813132702,0.00019405185108843607,6.0898671943944786e-5,6.0604505038256276e-5]
@test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, fit_mean=false)) ≈ [0.09219168665786258,0.0015654342453078724,0.00019403694017215876,6.088944186950046e-5,6.05771360378885e-5]
@test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, center_data=false, fit_mean=false)) ≈ [0.09360344864985332,0.0015354489715019735,0.0001784388515190763,4.744247354697125e-5,3.240223498703448e-5]
@test power(lombscargle(t, s, err)) ==
    power(lombscargle(t, measurement.(s, err)))

# Compare result of uncertainties with both methods (fast and non-fast).
srand(1)
errors = rand(0.1:1e-3:4.0, ntimes)
@test power(lombscargle(t, s, errors)) ≈ power(lombscargle(t, s, errors, fast = false)) atol = 0.3
@test power(lombscargle(t, s, errors, fit_mean = false)) ≈ power(lombscargle(t, s, errors, fit_mean = false, fast = false)) atol = 0.2

##################################################
### Testing utilities

# Test findmaxfreq and findmaxpower
@test findmaxfreq(pgram1) ≈ [31.997145470342]
@test findmaxfreq(pgram1, 0.965) ≈ [0.15602150741832602,31.685102455505348,31.997145470342,63.52622641842902,63.838269433265665]
@test findmaxperiod(pgram1) ≈ 1./findmaxfreq(pgram1)
@test findmaxperiod(pgram1, 0.965) ≈ 1./findmaxfreq(pgram1, 0.965)
t = linspace(0.01, 10pi, 1001)
s = sinpi.(2t) + cospi.(4t)
p = lombscargle(t, s, maximum_frequency=4)
@test findmaxfreq(p, [0.9, 1.1]) ≈ [1.0029954048320957]
@test findmaxfreq(p, [1.9, 2.1]) ≈ [2.002806697267899]
@test findmaxperiod(p, [1/0.9, 1/1.1]) ≈ 1./findmaxfreq(p, [0.9, 1.1])
@test findmaxperiod(p, [1/1.9, 1/2.1]) ≈ 1./findmaxfreq(p, [1.9, 2.1])
@test findmaxpower(pgram1) ≈ 0.9695017551608017

# Test autofrequency function
@test LombScargle.autofrequency(t) ≈ 0.003184112396292367:0.006368224792584734:79.6824127172165
@test LombScargle.autofrequency(t, minimum_frequency=0) ≈ 0.0:0.006368224792584734:79.6792286048202
@test LombScargle.autofrequency(t, maximum_frequency=10) ≈ 0.003184112396292367:0.006368224792584734:9.99492881196174
# This last test also makes sure that `freq` and `power` fields of a Periodogram
# object can have different type.
@test freq(lombscargle(1:11, big.(sin.(1:11)))) ≈ 0.01:0.02:2.75

# Test probabilities and FAP
t = collect(linspace(0.01, 10pi, 101))
s = sin.(t)
for norm in (:standard, :Scargle, :HorneBaliunas, :Cumming)
    P = lombscargle(t, s, normalization = norm)
    for z_0 in (0.1, 0.5, 0.9)
        @test prob(P, probinv(P, z_0)) ≈ z_0
        @test fap(P,  fapinv(P, z_0)) ≈ z_0
    end
end
P = lombscargle(t, s, normalization = :log)
@test_throws ErrorException prob(P, 0.5)
@test_throws ErrorException probinv(P, 0.5)

# Test model function
@test s ≈ LombScargle.model(t, s, 1/2pi, center_data=false, fit_mean=false)
s = sinpi.(t) + pi*cospi.(t) + e
@test s ≈ LombScargle.model(t, measurement.(s, ones(s)), 0.5)

# Test add_at!
a = ones(Int, 3)
LombScargle.add_at!(a, [3, 1, 3, 1, 2], 1:5)
@test a ≈ [7, 6, 5]

# Test extirpolation function
x = collect(linspace(0, 10))
y = sin.(x)
vec13 = Vector{Complex{Float64}}(13)
LombScargle.extirpolate!(vec13, x, y, 13)
@test vec13 ≈ [0.39537718210649553,3.979484140636793,4.833090108345013,0.506805556164743,-3.828112427525919,-4.748341359084166,-1.3022050566901917,3.3367666084342256,5.070478111668922,1.291245296032218,-0.8264466821981216,0.0,0.0]
x = collect(linspace(0, 10))
y = sin.(x)
@test LombScargle.extirpolate!(Vector{Complex{Float64}}(11), x, y, 11) ≈ vec13[1:11]

# Test trig_sum!
FFTW.set_num_threads(2)
info("Multi-threading in FFTW enabled")
x = collect(linspace(0, 10))
y = sin.(x)
N = 10
Nfft = nextpow2(5N)
bfft_vec = Vector{Complex{Float64}}(Nfft)
p = plan_bfft(bfft_vec)
grid = similar(bfft_vec)
C, S = LombScargle.trig_sum!(grid, bfft_vec, p, x, y, 1, N, Nfft)
@test S ≈ [0.0,0.3753570125888358,0.08163980192703546,-0.10139634351774979,-0.4334223744905633,-2.7843373311769777,0.32405810159838055,0.05729663600471602,-0.13191736591325876,-0.5295781583202946]
@test C ≈ [8.708141477890015,-0.5402668064176129,-0.37460815054027985,-0.3793457539084364,-0.5972351546196192,14.612204307982497,-0.5020253140297526,-0.37724493022381034,-0.394096831764578,-0.6828241623474718]

# Test bootstrap
srand(1)
@test LombScargle.bootstrap(5, x, y).p ≈
    [0.25956949678034225, 0.2360115683328911, 0.22016267001066891, 0.1665406952388801,
     0.12516095308735742]
b = LombScargle.bootstrap(50, x, y, fast = true)
@test fap(b, fapinv(b, 0.02)) ≈ 0.02
err = collect(linspace(0.5, 1.5))
b = LombScargle.bootstrap(50, x, measurement.(y, err), fast = true)
@test fap(b, fapinv(b, 0.02)) ≈ 0.02
