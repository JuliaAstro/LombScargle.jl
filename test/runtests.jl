using LombScargle
using Measurements
using Base.Test

ntimes = 1001
# Observation times
t = linspace(0.01, 10pi, ntimes)
# Randomize times
t += step(t)*rand(ntimes)
t = collect(t)
# The signal
s = sinpi(t) + cospi(2t) + rand(ntimes)
# Frequency grid
nfreqs = 10000
freqs = linspace(0.01, 3, nfreqs)
# Randomize frequency grid
freqs += step(freqs)*rand(nfreqs)
freqs = collect(freqs)
# Use "freqpower" just to call that function and increase code coverage.
# "autofrequency" function is tested below.
@test_approx_eq_eps freqpower(lombscargle(t, s, frequencies=freqs, fit_mean=false))[2] freqpower(lombscargle(t, s, frequencies=freqs, fit_mean=true))[2] 5e-3

# Simple signal, without any randomness
t = collect(linspace(0.01, 10pi, ntimes))
s = sin(t)
pgram1 = lombscargle(t, s, fit_mean=false)
pgram2 = lombscargle(t, s, fit_mean=true)
@test_approx_eq_eps power(pgram1) power(pgram2) 2e-7
pgram3 = lombscargle(t, s, center_data=false, fit_mean=false)
pgram4 = lombscargle(t, s, center_data=false, fit_mean=true)
@test_approx_eq_eps power(pgram3) power(pgram4) 3e-7

# Test the values in order to prevent wrong results in both algorithms
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=true))  [0.029886871262324886,0.0005456198989410226,1.912507742056023e-5, 4.54258409531214e-6, 1.0238342782430832e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, fit_mean=false)) [0.02988686776042212, 0.0005456197937194695,1.9125076826683257e-5,4.542583863304549e-6,1.0238340733199874e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, center_data=false, fit_mean=true))  [0.029886871262325004,0.0005456198989536703,1.9125077421448458e-5,4.5425840956285145e-6,1.023834278337881e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, center_data=false, fit_mean=false)) [0.029886868328967767,0.0005456198924872134,1.9125084251687147e-5,4.542588504467314e-6,1.0238354525870936e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="model")) [0.030807614469885718,0.0005459177625354441,1.9125443196143085e-5,4.54260473047638e-6,1.0238447607164715e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="log")) [0.030342586720560734,0.0005457688036440774,1.9125260307148152e-5,4.542594412890309e-6,1.0238395194654036e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="psd")) [7.474096700871138,0.1364484040771917,0.004782791641128195,0.0011360075968541799,0.002560400630125523]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="Scargle")) [0.029886871262324904,0.0005456198989410194,1.912507742056126e-5,4.54258409531238e-6,1.0238342782428552e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="HorneBaliunas")) [14.943435631162451,0.2728099494705097,0.009562538710280628,0.00227129204765619,0.005119171391214276]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, normalization="Cumming")) [15.372999620472974,0.2806521440709115,0.009837423440787873,0.0023365826071340815,0.005266327088140394]
@test_throws ErrorException lombscargle(t, s, frequencies=0.2:0.2:1, normalization="foo")

# Test signal with uncertainties
err = collect(linspace(0.5, 1.5, ntimes))
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=true))  [0.06659683848818687,0.09361438921056377,0.006815919926284516,0.0016347568319229223,0.0005385706045724484,0.00021180745624003642,0.00010539881897690343,7.01610752020905e-5,5.519295593372065e-5,4.339157565349008e-5]
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=false)) [0.0692080444168825,0.09360343748833044,0.006634919855866448,0.0015362074096358814,0.000485825178683968,0.000181798596773626,8.543735242380012e-5,5.380000205539795e-5,4.010727072660524e-5,2.97840883747593e-5]
@test power(lombscargle(t, s, err)) ==
    power(lombscargle(t, measurement(s, err)))

# Test fast method
t = linspace(0.01, 10pi, ntimes)
s = sin(t)
pgram5 = lombscargle(t, s, maximum_frequency=30, use_fft=true)
pgram6 = lombscargle(t, s, maximum_frequency=30, use_fft=false)
@test_approx_eq_eps power(pgram5) power(pgram6) 2e-6
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, use_fft=true, fit_mean=true))  [0.029886871262324963,0.0005453105325981758,1.9499330722168046e-5,2.0859593514888897e-6,1.0129019249708592e-5]
@test_approx_eq power(lombscargle(t, s, frequencies=0.2:0.2:1, use_fft=true, fit_mean=false)) [0.02988686776042208,0.000545619793719469,1.912507682668341e-5,4.542583863304592e-6,1.0238340733199892e-5]
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.2:0.2:1, use_fft=true, fit_mean=true))  [0.09230959166317655,0.0015654929813132702,0.00019405185108843607,6.0898671943944786e-5,6.0604505038256276e-5]
@test_approx_eq power(lombscargle(t, s, err, frequencies=0.2:0.2:1, use_fft=true, fit_mean=false)) [0.09219168665786258,0.0015654342453078724,0.00019403694017215876,6.088944186950046e-5,6.05771360378885e-5]

##################################################
### Testing utilities

# Test findmaxfreq and findmaxpower
@test_approx_eq findmaxfreq(pgram1)        [31.997145470342]
@test_approx_eq findmaxfreq(pgram1, 0.965) [0.15602150741832602,31.685102455505348,31.997145470342,63.52622641842902,63.838269433265665]
@test_approx_eq findmaxpower(pgram1) 0.9695017551608017

# Test autofrequency function
@test_approx_eq LombScargle.autofrequency(t)                       0.003184112396292367:0.006368224792584734:79.6824127172165
@test_approx_eq LombScargle.autofrequency(t, minimum_frequency=0)                   0.0:0.006368224792584734:79.6792286048202
@test_approx_eq LombScargle.autofrequency(t, maximum_frequency=10) 0.003184112396292367:0.006368224792584734:9.99492881196174

# Test probabilities and FAP
t = collect(linspace(0.01, 10pi, 101))
s = sin(t)
for norm in ("standard", "Scargle", "HorneBaliunas", "Cumming")
    P = lombscargle(t, s, normalization = norm)
    for z_0 in (0.1, 0.5, 0.9)
        @test_approx_eq prob(P, probinv(P, z_0)) z_0
        @test_approx_eq fap(P,  fapinv(P, z_0))  z_0
    end
end
P = lombscargle(t, s, normalization = "log")
@test_throws ErrorException prob(P, 0.5)
@test_throws ErrorException probinv(P, 0.5)

# Test extirpolation function
x = linspace(0, 10)
y = sin(x)
@test_approx_eq LombScargle.extirpolate(x, y) [0.39537718210649553,3.979484140636793,4.833090108345013,0.506805556164743,-3.828112427525919,-4.748341359084166,-1.3022050566901917,3.3367666084342256,5.070478111668922,1.291245296032218,-0.8264466821981216,0.0,0.0]
@test_approx_eq LombScargle.extirpolate(x, y, 11) LombScargle.extirpolate(x, y)[1:11]

# Test trig_sum
C, S = LombScargle.trig_sum(x, y, 1, 10)
@test_approx_eq S [0.0,0.3753570125888358,0.08163980192703546,-0.10139634351774979,-0.4334223744905633,-2.7843373311769777,0.32405810159838055,0.05729663600471602,-0.13191736591325876,-0.5295781583202946]
@test_approx_eq C [8.708141477890015,-0.5402668064176129,-0.37460815054027985,-0.3793457539084364,-0.5972351546196192,14.612204307982497,-0.5020253140297526,-0.37724493022381034,-0.394096831764578,-0.6828241623474718]

# Test choose_use_fft
@test LombScargle.choose_use_fft(300, true,  :maybe) == true
@test LombScargle.choose_use_fft(300, true,  true)   == true
@test LombScargle.choose_use_fft(300, true,  false)  == false
@test LombScargle.choose_use_fft(300, false, :maybe) == false
@test LombScargle.choose_use_fft(300, false, true)   == false
@test LombScargle.choose_use_fft(300, false, false)  == false

@test LombScargle.choose_use_fft(100, true,  :maybe) == false
@test LombScargle.choose_use_fft(100, true,  true)   == true
@test LombScargle.choose_use_fft(100, true,  false)  == false
@test LombScargle.choose_use_fft(100, false, :maybe) == false
@test LombScargle.choose_use_fft(100, false, true)   == false
@test LombScargle.choose_use_fft(100, false, false)  == false
