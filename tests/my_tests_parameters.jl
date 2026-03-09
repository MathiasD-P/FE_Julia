#####################################################################
# This is just a file to store the parameters for each testcase
# We merely list everything that's needed in each case.
# The parameters can be instantiated via make_my_tests_parameters.
# MAKE SURE TO CHOOSE A "GOOD" NAME WHEN ADDING A TESTCASE.
#####################################################################

using FE_Julia

function make_my_tests_parameters(testname::String)
    if occursin(r"test_OOA_DGStd_LinAdv_1D_p.*$", testname)
        p = parse(Int, match(r"test_OOA_DGStd_LinAdv_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="LinAdv",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        enodes="(10)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="upwind",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        a=1.2)

    elseif occursin(r"test_OOA_DGStd_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_OOA_DGStd_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        enodes="(10)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="ES_Chandrashekar_dissip",
                        ICname="GassnerEuler",
                        BCname="periodic",
                        sourcename = "GassnerEuler",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)

    elseif occursin(r"test_OOA_DGFluxDiff_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_OOA_DGFluxDiff_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="ES_Chandrashekar_dissip",
                        twoptfluxtype="EC_Chandrashekar",
                        ICname="GassnerEuler",
                        BCname="periodic",
                        sourcename = "GassnerEuler",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)

    elseif occursin(r"test_OOA_DGArtVisc_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_OOA_DGArtVisc_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+1)*")-GL",
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="ES_Chandrashekar_dissip",
                        AVcoeff="AVEC",
                        ICname="GassnerEuler",
                        BCname="periodic",
                        sourcename = "GassnerEuler",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)
    
    elseif occursin(r"test_ent_DGStd_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_ent_DGStd_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_Chandrashekar",
                        ICname="GaussianVelocity",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=17000,
                        dt=0.00001,
                        calc_entropy=true,
                        gamma=1.4)
                    
    elseif occursin(r"test_ent_DGFluxDiff_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_ent_DGFluxDiff_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_Chandrashekar",
                        twoptfluxtype="EC_Chandrashekar",
                        ICname="GaussianVelocity",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=17000,
                        dt=0.00001,
                        calc_entropy=true,
                        gamma=1.4)

    elseif occursin(r"test_ent_DGFluxDiff_Burgers_1D_p.*$", testname)
        p = parse(Int, match(r"test_ent_DGFluxDiff_Burgers_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_split",
                        twoptfluxtype="EC_split",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.00002,
                        calc_entropy=true)
    
    elseif occursin(r"test_ent_DGArtVisc_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_ent_DGArtVisc_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+2)*")-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=10,
                        numfluxtype="ES_Chandrashekar_dissip",
                        AVcoeff="AVdissip",
                        ICname="GaussianVelocity",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=17000,
                        dt=0.00001,
                        calc_entropy=true,
                        gamma=1.4)
    
    elseif occursin(r"test_visc_DGArtVisc_Chand_Euler_1D_p.*$", testname)
        p = parse(Int, match(r"test_visc_DGArtVisc_Chand_Euler_1D_p(.*)$", testname)[1])
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes="("*string(p+1)*")-GL",
                        qnodes="("*string(p+1)*")-GL",
                        enodes="(15)-GL", # WE USE THESE NODES FOR L2 PROJ OF ICS
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=3,
                        numfluxtype="ES_Chandrashekar_dissip",
                        AVcoeff="AVdissip",
                        ICname="ChanWave",
                        BCname="periodic",
                        gamma=1.4)

    else
        error("Undefined test case name!")
    end

    return param
end