using FE_Julia
using Test

include("validation_tests_parameters.jl")
include("validation_tests.jl")
include("validation_tests_aliasing.jl")

println("Starting Order of Accuracy Tests")

@testset "Order of Accuracy Tests" begin

    @testset "DGStd Tests" begin

        println("DGStd LinAdv (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_LinAdv_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.01

        println("DGStd LinAdv 1D (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_LinAdv_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

        println("DGStd Burgers 1D (Gassner manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_Burgers_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.01

        println("DGStd Burgers 1D (Gassner manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_Burgers_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

        println("DGStd Euler 1D (Euler manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_Chand_Euler_1D_b(4)-GLL_q(4)-GL"), 8)[end,end] ≈ -4.0 atol=0.02

        println("DGStd Euler 1D (Euler manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGStd_Chand_Euler_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

    end

    @testset "DGFluxDiff Tests" begin
        
        println("DGFluxDiff Burgers 1D (Gassner manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGFluxDiff_Burgers_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.01

        println("DGFluxDiff Burgers 1D (Gassner manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGFluxDiff_Burgers_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

        println("DGFluxDiff Euler 1D (Euler manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGFluxDiff_Chand_Euler_1D_b(4)-GLL_q(4)-GL"), 8)[end,end] ≈ -4.0 atol=0.02

        println("DGFluxDiff Euler 1D (Euler manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGFluxDiff_Chand_Euler_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.04    
    end

    @testset "DGAddRes Tests" begin

        println("DGAddRes Burgers 1D (Gassner manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGAddRes_Burgers_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.01

        println("DGAddRes Burgers 1D (Gassner manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGAddRes_Burgers_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

        println("DGAddRes Euler 1D (Euler manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGAddRes_Chand_Euler_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.06

        println("DGAddRes Euler 1D (Euler manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGAddRes_Chand_Euler_1D_b(5)-GLL_q(5)-GL"), 6)[end,end] ≈ -5.0 atol=0.06   
    end

    @testset "DGArtVisc Tests" begin

        println("DGArtVisc Burgers 1D (Gassner manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGArtVisc_Burgers_1D_b(4)-GLL_q(4)-GL"), 7)[end,end] ≈ -4.0 atol=0.01

        println("DGArtVisc Burgers 1D (Gassner manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGArtVisc_Burgers_1D_b(5)-GLL_q(5)-GL"), 7)[end,end] ≈ -5.0 atol=0.01

        println("DGArtVisc Euler 1D (Euler manufactured) (p3):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGArtVisc_Chand_Euler_1D_b(4)-GLL_q(4)-GL"), 8)[end,end] ≈ -4.0 atol=0.02

        println("DGArtVisc Euler 1D (Euler manufactured) (p4):")
        @test test_OOA(make_validation_tests_parameters("test_OOA_DGArtVisc_Chand_Euler_1D_b(5)-GLL_q(5)-GL"), 6)[end,end] ≈ -5.0 atol=0.06   
    end
end


println("Starting Entropy Stability Tests")

@testset "Entropy Stability Tests" begin

    @testset "DGFluxDiff Tests" begin

        println("DGFluxDiff Burgers 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGFluxDiff_Burgers_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-15

        println("DGFluxDiff Euler 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGFluxDiff_Chand_Euler_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-14
    end

    @testset "DGAddRes Tests" begin

        println("DGAddRes Burgers 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGAddRes_Burgers_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-14

        println("DGAddRes Euler 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGAddRes_Chand_Euler_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-14
    end

    @testset "DGArtVisc Tests" begin

        println("DGArtVisc Burgers 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGArtVisc_Burgers_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-14

        println("DGArtVisc Euler 1D (p4)")
        @test test_entropy(make_validation_tests_parameters("test_ent_DGArtVisc_Chand_Euler_1D_b(5)-GLL_q(5)-GL")) ≈ 0.0 atol=1e-14
    end
end


println("Starting Aliasing Tensor Tests")

@testset "Aliasing Tensor Tests" begin
    nodes = make_nodes("(9)-GLL"); tpflux = FE_Julia.SplitTPFlux(pi/4, 1-pi/4);
    @test test_Delta(nodes, tpflux, 6, 5) ≈ 0 atol=1e-13

    nodes = make_nodes("(12)-GLL"); tpflux = FE_Julia.SplitTPFlux(pi/4, 1-pi/4);
    @test test_Delta(nodes, tpflux, 11, 5) ≈ 0 atol=1e-13

    nodes = make_nodes("(9)-GL"); tpflux = FE_Julia.SplitTPFlux(pi/4, 1-pi/4);
    @test test_Delta(nodes, tpflux, 6, 5) ≈ 0 atol=1e-13

    nodes = make_nodes("(12)-GL"); tpflux = FE_Julia.SplitTPFlux(pi/4, 1-pi/4);
    @test test_Delta(nodes, tpflux, 11, 5) ≈ 0 atol=1e-13

    nodes = make_nodes("(9)-GL"); tpflux = FE_Julia.SplitTPFlux(11/19, 8/19); # Dealiasing flux for GL
    @test Delta_empirical(nodes, tpflux, 7,4,8) ≈ 0 atol=1e-13

    nodes = make_nodes("(12)-GL"); tpflux = FE_Julia.SplitTPFlux(14/25, 11/25); # Dealiasing flux for GL
    @test Delta_empirical(nodes, tpflux, 11,3,11) ≈ 0 atol=1e-13

    nodes = make_nodes("(9)-GLL"); tpflux = FE_Julia.SplitTPFlux(9/17, 8/17); # Dealiasing flux for GLL
    @test Delta_empirical(nodes, tpflux, 7,2,8) ≈ 0 atol=1e-13

    nodes = make_nodes("(12)-GLL"); tpflux = FE_Julia.SplitTPFlux(12/23, 11/23); # Dealiasing flux for GLL
    @test Delta_empirical(nodes, tpflux, 11,1,11) ≈ 0 atol=1e-13
end