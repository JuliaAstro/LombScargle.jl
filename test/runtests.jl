using LombScargle
using Measurements, FFTW, StableRNGs
using Test

Threads.nthreads() > 1 && (FFTW.set_num_threads(2); @info("Multi-threading enabled"))

ntimes = 1001
# Observation times
t = range(0.01, stop=10pi, length=ntimes)
# Signal
s = sin.(t)

pgram1 = @inferred(lombscargle(LombScargle.plan(t, s, fast = false, fit_mean=false)))
pgram2 = @inferred(lombscargle(LombScargle.plan(t, s, fast = false, fit_mean=true)))
pgram3 = @inferred(lombscargle(LombScargle.plan(t, s, fast = false, center_data=false, fit_mean=false)))
pgram4 = @inferred(lombscargle(LombScargle.plan(t, s, fast = false, center_data=false, fit_mean=true)))

@testset "Periodograms" begin
    @testset "Random stuff" begin
        rng = StableRNG(1)
        # Randomized times
        trandom = t .+ step(t) .* rand(rng, ntimes)
        # Randomized signal
        srandom = sinpi.(trandom) .+ cospi.(2trandom) .+ rand(rng, ntimes)
        # Frequency grid
        nfreqs = 10000
        freqs = range(0.01, stop=3, length=nfreqs)
        # Randomize frequency grid
        freqs += step(freqs) * rand(rng, nfreqs)
        # Use "freqpower" and "periodpower" just to call that function and increase code
        # coverage.  "autofrequency" function is tested below.
        prandom1 = lombscargle(trandom, srandom, frequencies=freqs, fit_mean=false)
        prandom2 = lombscargle(trandom, srandom, frequencies=freqs, fit_mean=true)
        @test freqpower(pgram1)[2] ≈ periodpower(pgram2)[2] atol = 6e-7
        @testset "Infinities" begin
            # Make sure there are no infinities in `pgram1'.  It seems to work only on
            # 64-bit systems.
            Sys.WORD_SIZE == 64 && @test(!(Inf in power(prandom1)) && !(Inf in power(prandom2)))
        end
    end

    @testset "fit_mean and center_data" begin
        @test power(pgram1) ≈ power(pgram2) atol = 5e-7
        @test power(pgram3) ≈ power(pgram4) atol = 5e-7
    end

    # Test the values in order to prevent wrong results in both algorithms
    @testset "Test values" begin
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
    end

    err = collect(range(0.5, stop=1.5, length=ntimes))

    @testset "Signal with uncertainties" begin
        @test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=true)) ≈ [0.06659683848818691,0.09230959166317589,0.006625919314669043,0.0015664010997692612,0.0005085442118408477,0.00019704659245878378,9.658452525613897e-5,6.331573873913433e-5,4.903871967643573e-5,3.7948448825374025e-5]
        @test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=false)) ≈ [0.0664002483305464,0.09219168665786254,0.006625915010614472,0.0015663421089042564,0.0005085109569237008,0.00019703233981948823,9.6577091433651e-5,6.33101344670203e-5,4.9033581990442793e-5,3.7944076990210425e-5]
        @test power(lombscargle(t, s, err, frequencies=0.1:0.1:1, fit_mean=false, center_data=false)) ≈ [0.06920814049261209,0.09360344864985352,0.006634919960009565,0.0015362072871144769,0.0004858250831632676,0.00018179850370583626,8.543727416919218e-5,5.379994730581837e-5,4.0107232867763e-5,2.9784059487535237e-5]
        @test power(lombscargle(t, s, err)) ==
            power(lombscargle(t, measurement.(s, err)))
        @test power(lombscargle(t, s, err, frequencies = 0.1:0.1:1, fast = false, normalization = :psd)) ≈ [21.851224476672318,30.287888352835566,2.1740438975167593,0.5139550589572747,0.16685947834022155,0.06465335925734642,0.031690545531213095,0.020774656147387098,0.01609019430987704,0.012451342926314715]
    end

    pgram5 = lombscargle(t, s, maximum_frequency=30, fast=true)
    pgram6 = lombscargle(t, s, maximum_frequency=30, fast=false)

    @testset "Fast method" begin
        @test power(pgram5) ≈ power(pgram6) atol = 3e-4
        @test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, fit_mean=true, flags = FFTW.MEASURE, timelimit = 5)) ≈
            [0.029886871262325053, 0.0005447325913220627, 2.1246300058375996e-5, 4.1259517049745417e-7, 5.04610747916143e-5]
        @test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, fit_mean=false, flags = FFTW.MEASURE, timelimit = 5)) ≈
            [0.0298868677604221, 0.0005447324642349108, 2.1246261449681898e-5, 4.125675324400738e-7, 5.042481545188997e-5]
        @test power(lombscargle(t, s, frequencies=0.2:0.2:1, fast=true, center_data=false, fit_mean=false, flags = FFTW.MEASURE, timelimit = 5)) ≈
            [0.029886868328967798, 0.0005447325727792022, 2.1246201600761156e-5, 4.1251689656274853e-7, 5.0422982073668846e-5]
        @test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, fit_mean=true)) ≈ [0.09230959166317658, 0.0015636510410697796, 0.00019210902295493893, 8.221668511639927e-5, 0.00021747947231895047]
        @test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, fit_mean=false)) ≈ [0.09219168665786256, 0.0015635926899499144, 0.00019209076832480172, 8.2184576575774e-5, 0.00021723818031507325]
        @test power(lombscargle(t, s, err, frequencies=0.2:0.2:1, fast=true, center_data=false, fit_mean=false)) ≈ [0.09360344864985316, 0.001533862558849818, 0.00017415674032714024, 5.1945371781883936e-5, 0.00011977935090727627]
        @test power(lombscargle(t, s, err)) ==
            power(lombscargle(t, measurement.(s, err)))
    end

    # Compare result of uncertainties with both methods (fast and non-fast).
    @testset "Fast and non-fast methods" begin
        rng = StableRNG(1)
        errors = rand(rng, 0.1:1e-3:4.0, ntimes)
        @test power(lombscargle(t, s, errors)) ≈ power(lombscargle(t, s, errors, fast = false)) atol = 0.3
        @test power(lombscargle(t, s, errors, fit_mean = false)) ≈ power(lombscargle(t, s, errors, fit_mean = false, fast = false)) atol = 0.2
    end

   @testset "Float32" begin
       @test typeof(period(lombscargle(Vector{Float32}(t), Vector{Float32}(s), fast = false))) == Vector{Float32}
   end
end

@testset "Utilities" begin
    @testset "findmaxpower, findmaxfreq, findmaxperiod" begin
        @test findmaxfreq(pgram1) ≈ [31.997145470342]
        @test findmaxfreq(pgram1, 0.965) ≈ [0.15602150741832602,31.685102455505348,31.997145470342,63.52622641842902,63.838269433265665]
        @test findmaxperiod(pgram1) ≈ 1 ./ findmaxfreq(pgram1)
        @test findmaxperiod(pgram1, 0.965) ≈ 1 ./ findmaxfreq(pgram1, 0.965)
        global s = sinpi.(2t) .+ cospi.(4t)
        p = lombscargle(t, s, maximum_frequency=4)
        @test findmaxfreq(p, [0.9, 1.1]) ≈ [1.0029954048320957]
        @test findmaxfreq(p, [1.9, 2.1]) ≈ [2.002806697267899]
        @test findmaxperiod(p, [1/0.9, 1/1.1]) ≈ 1 ./ findmaxfreq(p, [0.9, 1.1])
        @test findmaxperiod(p, [1/1.9, 1/2.1]) ≈ 1 ./ findmaxfreq(p, [1.9, 2.1])
        @test findmaxpower(pgram1) ≈ 0.9695017551608017
    end

    @testset "LombScargle.autofrequency" begin
        @test LombScargle.autofrequency(t) ≈ 0.003184112396292367:0.006368224792584734:79.6824127172165
        @test LombScargle.autofrequency(t, minimum_frequency=0) ≈ 0.0:0.006368224792584734:79.6792286048202
        @test LombScargle.autofrequency(t, maximum_frequency=10) ≈ 0.003184112396292367:0.006368224792584734:9.99492881196174
        # This last test also makes sure that `freq` and `power` fields of a Periodogram
        # object can have different type.
        if Threads.nthreads() > 1 && Sys.isapple()
            # This would make Julia crash, skip it.  See
            # https://github.com/JuliaLang/julia/issues/35702
            @test_skip freq(lombscargle(1:11, big.(sin.(1:11)))) ≈ 0.01:0.02:2.75
        else
            @test freq(lombscargle(1:11, big.(sin.(1:11)))) ≈ 0.01:0.02:2.75
        end
    end

    @testset "Probabilities and FAP" begin
        global t = collect(range(0.01, stop = 10pi, length = 101))
        global s = sin.(t)
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
    end

    @testset "LombScargle.model" begin
        @test s ≈ LombScargle.model(t, s, 1/2pi, center_data=false, fit_mean=false)
        global s = sinpi.(t) .+ pi .* cospi.(t) .+ ℯ
        @test s ≈ LombScargle.model(t, measurement.(s, fill(1, size(s))), 0.5)
    end

    @testset "LombScargle.add_at!" begin
        a = ones(Int, 3)
        LombScargle.add_at!(a, [3, 1, 3, 1, 2], 1:5)
        @test a ≈ [7, 6, 5]
    end

    @testset "LombScargle.extirpolate!" begin
        x = collect(range(0, stop = 10, length = 50))
        y = sin.(x)
        vec13 = Vector{Complex{Float64}}(undef, 13)
        LombScargle.extirpolate!(vec13, x, y, 13)
        @test vec13 ≈ [0.39537718210649553,3.979484140636793,4.833090108345013,0.506805556164743,-3.828112427525919,-4.748341359084166,-1.3022050566901917,3.3367666084342256,5.070478111668922,1.291245296032218,-0.8264466821981216,0.0,0.0]
        x = collect(range(0, stop = 10, length = 50))
        y = sin.(x)
        @test LombScargle.extirpolate!(Vector{Complex{Float64}}(undef, 11), x, y, 11) ≈ vec13[1:11]
    end

    x = collect(range(0, stop = 10, length = 50))
    y = sin.(x)

    @testset "LombScargle.trig_sum!" begin
        N = 10
        Nfft = nextpow(2, 5N)
        fftgrid = Vector{Complex{Float64}}(undef, N)
        bfft_vec = Vector{Complex{Float64}}(undef, Nfft)
        p = plan_bfft(bfft_vec)
        grid = similar(bfft_vec)
        C, S = LombScargle.trig_sum!(grid, fftgrid, bfft_vec, p, x, y, 1, N, Nfft, minimum(x))
        @test S ≈ [0.0,0.3753570125888358,0.08163980192703546,-0.10139634351774979,-0.4334223744905633,-2.7843373311769777,0.32405810159838055,0.05729663600471602,-0.13191736591325876,-0.5295781583202946]
        @test C ≈ [8.708141477890015,-0.5402668064176129,-0.37460815054027985,-0.3793457539084364,-0.5972351546196192,14.612204307982497,-0.5020253140297526,-0.37724493022381034,-0.394096831764578,-0.6828241623474718]
    end

    @testset "Bootstrap" begin
        rng = StableRNG(1)
        plan = LombScargle.plan(x, y)
        # Fill the periodogram in the plan with random numbers, to remove
        # possible NaNs, which would make the check below harder.  Zeroing the
        # vector would be uninteresting.
        plan.P .= rand.()
        P = copy(plan.P)
        @test LombScargle.bootstrap(rng, 5, plan).p ≈
            [0.2583163570869385, 0.24972129609003385, 0.23092196417031927, 0.18502993723773883, 0.1789661587851332]
        # Make sure the periodogram in the plan didn't change.
        @test plan.P == P
        rng = StableRNG(1)
        plan = LombScargle.plan(x, y; fit_mean = false)
        plan.P .= rand.()
        P = copy(plan.P)
        @test LombScargle.bootstrap(rng, 5, plan).p ≈
            [0.2580212594813987, 0.24897444742007394, 0.23090538831280463, 0.18492853793841382, 0.17895144706266247]
        @test plan.P == P
        err = collect(range(0.5, stop = 1.5, length = 50))
        rng = StableRNG(1)
        plan = LombScargle.plan(x, measurement.(y, err); fast = true)
        plan.P .= rand.()
        P = copy(plan.P)
        b = LombScargle.bootstrap(rng, 50, plan)
        @test fap(b, fapinv(b, 0.02)) ≈ 0.02
        @test plan.P == P
        rng = StableRNG(1)
        plan = LombScargle.plan(x, measurement.(y, err); fast = false)
        plan.P .= rand.()
        P = copy(plan.P)
        b = LombScargle.bootstrap(rng, 50, plan)
        @test fap(b, fapinv(b, 0.02)) ≈ 0.02
        @test plan.P == P
        rng = StableRNG(1)
        plan = LombScargle.plan(x, measurement.(y, err), fast = false, fit_mean = false)
        plan.P .= rand.()
        P = copy(plan.P)
        @test fapinv(LombScargle.bootstrap(rng, 50, plan),
                     0.2) ≈ 0.23917691901908134
        @test plan.P == P
        rng = StableRNG(1)
        plan = LombScargle.plan(x, y; fast=false, fit_mean=false)
        plan.P .= rand.()
        P = copy(plan.P)
        @test fapinv(LombScargle.bootstrap(rng, 1000, plan), 0.2) ≈
            0.2157617143004672
        @test plan.P == P
    end
end
