# test/types.jl --- type definitions, constructors

struct AbsorbanceSpectrum{xT<:Number,yT<:Number} <: EvenSpacedSpectrum{xT,yT}
    epoints::Tuple{xT,xT}
    absorbances::Vector{yT}
end
SpectrumBase.endpoints(s::AbsorbanceSpectrum) = s.epoints
SpectrumBase.intensities(s::AbsorbanceSpectrum) = s.absorbances

struct FastPulseSpectrum{xT<:Number,yT<:Number} <: AbstractSpectrum{xT,yT}
    xvalues::Vector{xT}
    intensities::Vector{yT}
end
SpectrumBase.range(s::FastPulseSpectrum) = s.xvalues
SpectrumBase.intensities(s::FastPulseSpectrum) = s.intensities
SpectrumBase.endpoints(s::FastPulseSpectrum) = begin
    v = range(s)
    (first(v), last(v))
end

@testset "Spectrum types" begin
    # Write test code here
    
    xv = (700.0, 250.0); rg = 700.0:-0.5:250.0
    yvals = randn(length(rg))
    testspec1d = AbsorbanceSpectrum(xv, yvals)

    @test AbsorbanceSpectrum <: EvenSpacedSpectrum
    @test AbsorbanceSpectrum <: AbstractSpectrum
    @test range(testspec1d) == rg
    @test step(testspec1d) == -0.5
    @test length(testspec1d) == length(rg) == length(yvals)
    @test SpectrumBase.isevenspaced(testspec1d) == true
    @test SpectrumBase.isevenspaced(AbsorbanceSpectrum) == true

    xv2 = randn(100); yvals2 = randn(100)
    testspec1d2 = FastPulseSpectrum(xv2, yvals2)

    @test FastPulseSpectrum <: AbstractSpectrum
    @test !SpectrumBase.isevenspaced(FastPulseSpectrum)
    @test !SpectrumBase.isevenspaced(testspec1d2)
    @test length(testspec1d2) == 100
    @test_throws MethodError step(testspec1d2)
end

@testset "Lineshapes" begin
    # create the lineshapes and define spectra based on them
    ## Gaussian
    @test_throws ArgumentError Gaussian(0.0, 0.0, 1.0)
    @test_throws ArgumentError Gaussian(0.0, -1.0, 1.0)
    gs = Gaussian(0.0, 2.5, 3.0)
    @test lscenter(gs) == 0.0
    @test scale(gs) == 2.5

    evalpoints = -5.0:0.01:5.0
    fp = (first(evalpoints), last(evalpoints))
    gaussian_y = map(evalpoints) do x
        evaluate(gs, x)
    end
    gs_absspec = AbsorbanceSpectrum((-5.0, 5.0), gaussian_y)
    @test length(gs_absspec) == length(evalpoints)
    @test endpoints(gs_absspec) == fp

    ## Cauchy
    @test_throws ArgumentError Cauchy(0.0, 0.0, 1.0)
    @test_throws ArgumentError Cauchy(0.0, -1.0, 1.0)
    lt = Cauchy(0.0, 2.5, 3.0)
    @test lscenter(lt) == 0.0
    @test scale(lt) == 2.5

    cauchy_y = map(evalpoints) do x
        evaluate(lt, x)
    end
    lt_absspec = AbsorbanceSpectrum((-5.0, 5.0), cauchy_y)
    @test length(lt_absspec) == length(evalpoints)
    @test endpoints(lt_absspec) == fp

    ## RaisedCosine
    @test_throws ArgumentError RaisedCosine(0.0, 0.0, 1.0)
    @test_throws ArgumentError RaisedCosine(0.0, -1.0, 1.0)
    rc = RaisedCosine(0.0, 2.5, 3.0)
    @test lscenter(rc) == 0
    @test scale(rc) == 2.5

    rcosine_y = map(evalpoints) do x
        evaluate(rc, x)
    end
    rc_absspec = AbsorbanceSpectrum((-5.0, 5.0), rcosine_y)
    @test length(rc_absspec) == length(evalpoints)
    @test endpoints(rc_absspec) == fp
end
