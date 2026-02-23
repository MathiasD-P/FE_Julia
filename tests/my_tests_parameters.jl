#####################################################################
# This is just a file to store the parameters for each testcase
# We merely list everything that's needed in each case.
# The parameters can be instantiated via make_my_tests_parameters.
# MAKE SURE TO CHOOSE A "GOOD" NAME WHEN ADDING A TESTCASE.
#####################################################################

using FE_Julia

function make_my_tests_parameters(testname::String)
    if testname == "test_OOA_DGStd_LinAdv_1D"
        param = parameters(
                        pdetype="LinAdv",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="(6)-GL",
                        qnodes="(7)-GL",
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

    elseif testname == "test_OOA_DGStd_Chand_Euler_1D"
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGStd",
                        bnodes="(5)-GL",
                        qnodes="(6)-GL",
                        enodes="(10)-GL",
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=2,
                        numfluxtype="EC_Chandrashekar",
                        ICname="IsentropicDensityWave",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=10000,
                        dt=0.0001,
                        OOAtest=true,
                        calc_entropy = true,
                        gamma=1.4)

    elseif testname == "test_OOA_DGFluxDiff_Chand_Euler_1D"
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="(5)-GL",
                        qnodes="(6)-GL",
                        enodes="(10)-GL",
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
                        Nsteps=4000,
                        dt=0.00025,
                        OOAtest=true,
                        calc_entropy = true,
                        gamma=1.4)
                    
    elseif testname == "test_ent_DGFluxDiff_Chand_Euler_1D"
        param = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="(4)-GL",
                        qnodes="(5)-GL",
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

    elseif testname == "test_ent_DGFluxDiff_Burgers_1D"
        param = parameters(
                        pdetype="Burgers",
                        dim=1,
                        dgtype="DGFluxDiff",
                        bnodes="(7)-GL",
                        qnodes="(8)-GL",
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
    end

    return param
end