# THIS ASSUMES CONSTANT LENGTH ELEMENTS

mutable struct Burgulence
    FFTpts::Matrix{Float64} # points for FFT evaluation
    FFTrefpts::Matrix{Float64} # ref points for FFT evaluation
    IC::Matrix{Float64} # Burgulence at basis nodes
    chimap::Matrix{Float64} # mapping from basis nodes to FFPts

    function Burgulence(bnodes::Tensorprod_nodes, Nel::Integer, kmax::Real, seednum=1234)
        if bnodes.dim != 1
            error("Basis nodes must be of dimension 1 to initialize Burgulence!")
        end

        inodes = deepcopy(bnodes)
        inodes.name = "GL" # We always use GL points for interpolation nodes
        irefpts = evaluate(inodes)
        brefpts = evaluate(bnodes)
        chi_itob = Lagrange_Vandermonde1D(brefpts, irefpts)

        # We assemble equidistant points for FFT evaluation
        FFTptsinit = range(-0.5, 0.5, numnodes(bnodes) * Nel + 1)
        FFTptsinit = FFTptsinit[1:end-1]
        
        # We make brown noise from white noise
        Random.seed!(seednum)
        uwhite = rand(Float64, size(FFTptsinit))
        k = fftfreq(length(FFTptsinit), 1/(FFTptsinit[1] - FFTptsinit[2]))
        uhat = fft(uwhite)
        uhat .= uhat ./ abs.(k)
        uhat .= uhat ./ sqrt(length(uhat))
        uhat[abs.(k) .> kmax] .= 0.0
        uhat[1] = 1.0

        # Evaluate bpts and ipts
        dx = 1 / Nel
        endpts = range(-0.5, 0.5, Nel+1)
        endpts = endpts[1:end-1]
        ipts = reduce(vcat, dx .* irefpts .+ endpts[ielem] for ielem in 1:Nel)

        # We sample at the basis points and normalize
        IC = nufft1d2(2*pi .* ipts, 1, 1e-9, uhat, modeord=1)
        IC = real(IC)
        IC = reshape(IC, (length(IC),1))
        IC = block_matmul(chi_itob, IC, Nel)

        # We assemble final FFTpts and mapping
        FFTrefpts = range(0.0, 1.0, numnodes(bnodes) + 1)
        FFTrefpts = FFTrefpts[1:end-1]
        FFTrefpts = reshape(FFTrefpts .+ 0.5 .* (FFTrefpts[2] - FFTrefpts[1]), (length(FFTrefpts),1))

        chimap = Lagrange_Vandermonde1D(FFTrefpts, brefpts)
        FFTpts = reduce(vcat, dx .* FFTrefpts .+ endpts[ielem] for ielem in 1:Nel)

        new(
            FFTpts,
            FFTrefpts,
            IC,
            chimap
        )
    end
end