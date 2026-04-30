#####################################################################
# This is just a file to store the parameters for each testcase
# We merely list everything that's needed in each case.
# The parameters can be instantiated via make_my_tests_parameters.
# MAKE SURE TO CHOOSE A "GOOD" NAME WHEN ADDING A TESTCASE.
#####################################################################

using FE_Julia

function make_my_tests_parameters(testname::String)
    if occursin(r"test_OOA_DGStd_LinAdv_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGStd_LinAdv_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="LinAdv",
                        dim=1,
                        dgtype="DGStd",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
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

    elseif occursin(r"test_OOA_DGStd_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGStd_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGStd",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="LF",
                        ICname="GassnerBurgers",
                        BCname="periodic",
                        sourcename = "GassnerBurgers",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true)

    elseif occursin(r"test_OOA_DGStd_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGStd_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
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
    
    elseif occursin(r"test_OOA_DGFluxDiff_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGFluxDiff_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="LF",
                        twoptfluxtype="EC_split",
                        ICname="GassnerBurgers",
                        BCname="periodic",
                        sourcename = "GassnerBurgers",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true)

    elseif occursin(r"test_OOA_DGFluxDiff_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGFluxDiff_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
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
    
    elseif occursin(r"test_OOA_DGArtVisc_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGArtVisc_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="LF",
                        AVcoeff="AVEC",
                        ICname="GassnerBurgers",
                        BCname="periodic",
                        sourcename = "GassnerBurgers",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true)

    elseif occursin(r"test_OOA_DGArtVisc_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGArtVisc_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
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
    
    elseif occursin(r"test_OOA_DGAddRes_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGAddRes_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGAddRes",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="LF",
                        Rescorr="RescorrEC",
                        ICname="GassnerBurgers",
                        BCname="periodic",
                        sourcename = "GassnerBurgers",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)
    
    elseif occursin(r"test_OOA_DGAddRes_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_OOA_DGAddRes_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGAddRes",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="ES_Chandrashekar_dissip",
                        Rescorr="RescorrEC",
                        ICname="GassnerEuler",
                        BCname="periodic",
                        sourcename = "GassnerEuler",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        gamma=1.4)
    
    elseif occursin(r"test_ent_DGStd_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGStd_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
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
                    
    elseif occursin(r"test_ent_DGFluxDiff_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGFluxDiff_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
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

    elseif occursin(r"test_ent_DGFluxDiff_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGFluxDiff_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
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
    
    elseif occursin(r"test_ent_DGArtVisc_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGArtVisc_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_Chandrashekar",
                        AVcoeff="AVEC",
                        ICname="GaussianVelocity",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=17000,
                        dt=0.00001,
                        calc_entropy=true,
                        gamma=1.4)
    
    elseif occursin(r"test_ent_DGArtVisc_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGArtVisc_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_split",
                        AVcoeff="AVEC",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.00002,
                        calc_entropy=true)
    
    elseif occursin(r"test_ent_DGAddRes_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGAddRes_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGAddRes",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_Chandrashekar",
                        Rescorr="RescorrEC",
                        ICname="GaussianVelocity",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=17000,
                        dt=0.00001,
                        calc_entropy=true,
                        gamma=1.4)
    
    elseif occursin(r"test_ent_DGAddRes_Burgers_1D_.*$", testname)
        nodeparse = match(r"test_ent_DGAddRes_Burgers_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGAddRes",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=20,
                        numfluxtype="EC_split",
                        Rescorr="RescorrEC",
                        ICname="exp_1state",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.00002,
                        calc_entropy=true)
    
    elseif occursin(r"test_visc_DGArtVisc_Chand_Euler_1D_.*$", testname)
        nodeparse = match(r"test_visc_DGArtVisc_Chand_Euler_1D_b(.*)_q(.*)$", testname)
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGArtVisc",
                        bnodes=nodeparse[1],
                        qnodes=nodeparse[2],
                        enodes="(12)-GL", # WE USE THESE NODES FOR L2 PROJ OF ICS
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