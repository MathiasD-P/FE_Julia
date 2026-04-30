using FE_Julia
using LinearAlgebra
using Plots

function make_parameters(bnodes, qnodes)
    params = parameters(
                        pdetype="EulerPerfGas",
                        dim=1,
                        bnodes=bnodes,
                        qnodes=qnodes,
                        enodes=bnodes[1:3] * "-GL", # just enough for exact integration
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=1,
                        numfluxtype="ES_Chandrashekar_dissip",
                        twoptfluxtype="EC_Chandrashekar",
                        Rescorr="RescorrEC",
                        AVcoeff="AVEC",
                        ICname="GassnerEuler",
                        sourcename="GassnerEuler",
                        BCname="periodic",
                        gamma=1.4)

    return params
end

function L2norm(u, chie, we, dg)
    norm = (FE_Julia.block_matmul(chie, u, dg.mesh.Nel)).^2

    for ielem = 1:dg.mesh.Nel
        index = 1+length(we)*(ielem-1):length(we)*ielem
        @views norm[index,:] .= norm[index,:] .* dg.mesh.J[ielem]
    end

    return sum(sqrt.(sum(FE_Julia.block_matmul(Diagonal(we), norm, dg.mesh.Nel), dims=1)))
end

function research_test(bnodes, qnodes, DGnames, colors, Nrefinements)

    params = make_parameters(bnodes, qnodes)
    chie, we, _ = FE_Julia.extract_volume_quadrature(make_nodes(params.bnodes), make_nodes(params.bnodes))

    reserrors = Dict(mydg => zeros(Nrefinements, 3) for mydg in DGnames)

    for i in 1:Nrefinements

        # Standard DG residual
        params.dgtype = "DGArtVisc"
        params.AVcoeff = "NoAV"
        u0, BChandler, dg = set_up_problem(params)
        stdresidual = similar(u0)
        build_residual!(stdresidual, u0, 0.0, BChandler, dg, params)
        params.AVcoeff = "AVEC"

        # iterate through all dg schemes
        for mydg in DGnames
            # Standard DG residual
            params.dgtype = mydg
            u0, BChandler, dg = set_up_problem(params)
            residual = similar(u0)
            build_residual!(residual, u0, 0.0, BChandler, dg, params)

            reserrors[mydg][i,2] = L2norm(residual .- stdresidual, chie, we, dg)
            reserrors[mydg][i,1] = dg.DOF
        end

        params.Neldim = round(params.Neldim * 1.5)
    end

    # Let's compute the OOA
    for mydg in DGnames
        reserrors[mydg][2:end,3] = (log10.(reserrors[mydg][2:end,2]) .- log10.(reserrors[mydg][1:end-1,2])) ./ (log10.(reserrors[mydg][2:end,1]) .- log10.(reserrors[mydg][1:end-1,1]))
    end

    # Now, we can plot the convergence of the residuals
    plt = plot(title= "PDE: " * params.pdetype * "; Basis: " * bnodes * "; Quad: " * qnodes, xlabel="log10(DOF)", ylabel="log10(L2)")
    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = myDG * " (" * params.AVcoeff * ")"
        elseif myDG == "DGAddRes"
            mylabel = myDG * " (" * params.Rescorr * ")"
        else
            mylabel = myDG
        end
        scatter!(plt, log10.(reserrors[myDG][:,1]), log10.(reserrors[myDG][:,2]), color=colors[colori], label=mylabel)

        println(myDG)
        display(reserrors[myDG])
    end
    display(plt)

    return reserrors
end

#%% Collocated GLL
bnodes = "(4)-GLL"
qnodes = "(4)-GLL"
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"] # colors used for plotting
Nrefinements = 15

reserrors = research_test(bnodes, qnodes, DGnames, colors, Nrefinements);

#%% Collocated GLL (DOESN'T WORK BECAUSE OF NUMFLUX)
bnodes = "(4)-GL"
qnodes = "(5)-GL"
DGnames = ["DGFluxDiff", "DGArtVisc", "DGAddRes"]
colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"] # colors used for plotting
Nrefinements = 10

reserrors = research_test(bnodes, qnodes, DGnames, colors, Nrefinements);