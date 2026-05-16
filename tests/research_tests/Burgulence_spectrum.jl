using FE_Julia
using FFTW
using FINUFFT
using LaTeXStrings
using Plots

# Default plotting parameters
default(
    fontfamily="Computer Modern",
    titlefont = font(12, "Computer Modern"),
    guidefont = font(12, "Computer Modern"),
    tickfont = font(10, "Computer Modern"),
    legendfont = font(10, "Computer Modern")
)

function make_parameters(bnodes, qnodes, numflux)
    params = parameters(
                        pdetype="Burgers",
                        dim=1,
                        bnodes=bnodes,
                        qnodes=qnodes,
                        fnodes="(1)-GL",
                        refelem="interval",
                        domain="unit_interval_linear",
                        Neldim=50,
                        numfluxtype=numflux,
                        twoptfluxtype="EC_split",
                        Rescorr="RescorrEC",
                        AVcoeff="AVdissip",
                        ICname="Burgulence",
                        BCname="periodic",
                        ODE_solver="LSERK45",
                        Nsteps=3000, #2950
                        dt=0.0001,
                        kmax=50)

    return params
end

function research_test(bnodes, qnodes, DGnames, numflux, colors, filenames=nothing)

    # Initialize parameters
    params = make_parameters(bnodes, qnodes, numflux)

    # Make our Burgulence object
    myBurgulence = FE_Julia.Burgulence(make_nodes(params.bnodes), params.Neldim, params.kmax)

    # First we compute our results
    results = Dict()
    for myDG in DGnames
        params.dgtype = myDG
        results[myDG] = set_up_and_solve(params)
    end

    # For completeness, we show initial conditions in the spatial and the frequency domain
    k = fftfreq(length(myBurgulence.FFTpts), 1/(myBurgulence.FFTpts[1] - myBurgulence.FFTpts[2]))
    uhatIC = fft(FE_Julia.block_matmul(myBurgulence.chimap, myBurgulence.IC, params.Neldim))

    plt_IC_spat = plot(results[DGnames[1]]["dg"].bpts, myBurgulence.IC, xlabel=L"x", ylabel=L"u_0", legend=false, color=:red)

    plt_IC_freq = plot(fftshift(abs.(k)[k .!= 0]), fftshift(abs.(uhatIC[k .!= 0]).^2), xlabel=L"k", ylabel=L"E(k)", xscale = :log10, yscale = :log10, label=L"\hat{u}_0", color=:red)
    plot!(plt_IC_freq, fftshift(abs.(k[k .!= 0])), fftshift(1e4 .* abs.(k[k .!= 0]).^-2), color=:gray, linestyle=:dash, label=L"1/k^2", ylims=(1e-5,1e5))

    # Post-processing and plotting
    
    plt_all = plot(xscale = :log10,
               yscale = :log10,
               title= "Numflux: " * params.numfluxtype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               xlabel=L"k",
               ylabel=L"E(k)",
               legend=:bottomleft)

    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = "Artificial Viscosity"
        elseif myDG == "DGAddRes"
            mylabel = "Residual Correction"
        elseif myDG == "DGStd"
            mylabel = "DG"
        elseif myDG == "DGFluxDiff"
            mylabel = "Flux Differencing"
        end

        u = results[myDG]["solution"]
        uhat = fft(FE_Julia.block_matmul(myBurgulence.chimap, u, params.Neldim))
        plot!(plt_all, fftshift(abs.(k)[k .!= 0]), fftshift(abs.(uhat[k .!= 0]).^2), label=mylabel, color=colors[colori])
    end
    plot!(plt_all, fftshift(abs.(k[k .!= 0])), fftshift(1e4 .* abs.(k[k .!= 0]).^-2), color=:gray, linestyle=:dash, label=L"1/k^2", ylims=(1e-5,1e5))

    plt_zoom = plot(xscale = :log10,
               yscale = :log10,
               title= "Numflux: " * params.numfluxtype * "; Basis: " * bnodes * "; Quad: " * qnodes,
               xlabel=L"k",
               ylabel=L"k^2E(k)",
               legend=:bottomleft)

    for (colori, myDG) in enumerate(DGnames)
        if myDG == "DGArtVisc"
            mylabel = "Artificial Viscosity"
        elseif myDG == "DGAddRes"
            mylabel = "Residual Correction"
        elseif myDG == "DGStd"
            mylabel = "DG"
        elseif myDG == "DGFluxDiff"
            mylabel = "Flux Differencing"
        end

        u = results[myDG]["solution"]
        uhat = fft(FE_Julia.block_matmul(myBurgulence.chimap, u, params.Neldim))
        plot!(plt_zoom, fftshift(abs.(k)[k .!= 0]), fftshift(abs.(k)[k .!= 0].^2 .* abs.(uhat[k .!= 0]).^2), label=mylabel, color=colors[colori], xlims=(50,130))
    end

    if isnothing(filenames)
        display(plt_IC_spat)
        display(plt_IC_freq)
        display(plt_all)
        display(plt_zoom)

    else
        savefig(plt_IC_spat, filenames[1] * ".png")
        savefig(plt_IC_freq, filenames[2] * ".png")
        savefig(plt_all, filenames[3] * ".png")
        savefig(plt_zoom, filenames[4] * ".png")
    end

    return results
end

# # function for plotting only
# function plot_research_test(results, bnodes, qnodes, DGnames, numflux, colors)

#     # Initialize parameters
#     params = make_parameters(bnodes, qnodes, numflux)

#     # Make our Burgulence object
#     myBurgulence = FE_Julia.Burgulence(make_nodes(params.bnodes), params.Neldim, params.kmax)



#     # For completeness, we show initial conditions in the spatial and the frequency domain
#     k = fftfreq(length(myBurgulence.FFTpts), 1/(myBurgulence.FFTpts[1] - myBurgulence.FFTpts[2]))
#     uhatIC = fft(FE_Julia.block_matmul(myBurgulence.chimap, myBurgulence.IC, params.Neldim))

#     plt_IC_spat = plot(results[DGnames[1]]["dg"].bpts, myBurgulence.IC, xlabel=L"x", ylabel=L"u_0", legend=false, color=:black)
#     display(plt_IC_spat)
#     plt_IC_freq = plot(abs.(k)[k .!= 0], abs.(uhatIC[k .!= 0]).^2, xlabel=L"k", ylabel=L"E(k)", xscale = :log10, yscale = :log10, label=L"\hat{u}_0", color=:black)
#     plot!(plt_IC_freq, fftshift(abs.(k[k .!= 0])), fftshift(1e6 .* abs.(k[k .!= 0]).^-2), color=:gray, linestyle=:dash, label=L"1/k^2", xlims=(1,1000))
#     display(plt_IC_freq)


#     # Post-processing and plotting
#     plt_all = plot(xscale = :log10,
#                yscale = :log10,
#                title= "Numflux: " * params.numfluxtype * "; Basis: " * bnodes * "; Quad: " * qnodes,
#                xlabel=L"k",
#                ylabel=L"k^2E(k)",
#                legend=:bottomleft)

#     for (colori, myDG) in enumerate(DGnames)
#         u = results[myDG]["solution"]
#         uhat = fft(FE_Julia.block_matmul(myBurgulence.chimap, u, params.Neldim))

#         # compare with initial conditions
#         plt_IC_spat = plot(results[DGnames[1]]["dg"].bpts, myBurgulence.IC, xlabel=L"x", ylabel=L"u", label=L"u_0", color=:black, title=myDG)
#         plot!(plt_IC_spat, results[DGnames[1]]["dg"].bpts, u, label=L"u(t)", color=colors[colori])
#         display(plt_IC_spat)

#         plt_IC_freq = plot(abs.(k)[k .!= 0], abs.(uhatIC[k .!= 0]).^2, xlabel=L"k", ylabel=L"E(k)", xscale = :log10, yscale = :log10, label=L"\hat{u}_0", color=:black, title=myDG)
#         plot!(plt_IC_freq, abs.(k)[k .!= 0], abs.(uhat[k .!= 0]).^2, label=L"u(t)", color=colors[colori])
#         plot!(plt_IC_freq, fftshift(abs.(k[k .!= 0])), fftshift(1e6 .* abs.(k[k .!= 0]).^-2), color=:gray, linestyle=:dash, label=L"1/k^2", xlims=(1,1000), ylims=(1e-4,1e6))
#         display(plt_IC_freq)

#         # Now plot all DG results
#         plot!(plt_all, abs.(k)[k .!= 0], abs.(uhat[k .!= 0]).^2 ./ abs.(uhatIC)[k .!= 0].^2, label=myDG, color=colors[colori], xlims=(10,150), ylims=(100, 10000))
#     end

#     display(plt_all)

#     return
# end

# Collocated GL
bnodes = "(5)-GL"
qnodes = "(5)-GL"
numflux = "LF"
DGnames = ["DGAddRes", "DGArtVisc", "DGFluxDiff", "DGStd"]
colors = ["#009E73", "#56B4E9", "#E69F00", "#F0E442"] # colors used for plotting

results = research_test(bnodes, qnodes, DGnames, numflux, colors,
                        ("outputs/fig_BurgulenceIC", "outputs/fig_BurgulenceICfreq", "outputs/fig_Burgulencespectrum", "outputs/fig_Burgulencespectrumzoom"))

