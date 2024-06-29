#----------------------------------------------------------------------------------------------#
#                                 11-ORIG-julia-BaseVector.jl                                  #
#----------------------------------------------------------------------------------------------#
# A "from the ground up", pure, and plain julia  implementation  of  the  Taylor-Green  vortex #
# decay problem solver in LBM.                                                                 #
#----------------------------------------------------------------------------------------------#

"""
`init(ℙ::Type{<:AbstractFloat}, l::Int)::Dict{Symbol, Dict}`\n
Computes (i) types, (ii) case, (iii) lattice, and (iv) properties simulation parameters, and
returns as a `Dict{Symbol, Dict}`.
"""
function init(ℙ::Type{𝕋} where 𝕋<:AbstractFloat,    # The floating point precision
              l::Int)                               # log2(scale)
    𝕀, 𝕌    = ℙ == Float64 ? (Int64, UInt64) : (Int32, UInt32)
    scale   = 𝕌(1) << l
    chunk   = 𝕌(32)
    maxIt   = 𝕌(204800)
    NY = NX = scale * chunk
    nu      = ℙ(1.0/6.0)
    w0, w1, w2 = ℙ(4.0/9.0), ℙ(1.0/9.0), ℙ(1.0/36.0)
    return Dict{Symbol, Dict}(
        # Types
        :typ => Dict{Symbol, DataType}(
            :i   => 𝕀,
            :u   => 𝕌,
            :p   => ℙ,
        ),
        # Case parameters
        :cas => Dict{Symbol, 𝕌}(
            :sca => scale,
            :NX  => NX,
            :NY  => NY,
            :IT  => 𝕌(round(maxIt / scale / scale)),
        ),
        # Lattice Stencil
        :lat => Dict{Symbol, Dict}(
            :int => Dict{Symbol, 𝕌}(:dim => 𝕌(2), :vel => 𝕌(9)),
            :flo => Dict{Symbol, ℙ}(:a => ℙ(√3/2), ),
            :vec => Dict{Symbol, Vector{ℙ}}(
                :w   => ℙ[w0, w1, w1, w1, w1, w2, w2, w2, w2],
                :ξx  => ℙ[+0, +1, +0, -1, +0, +1, -1, -1, +1],  # ℙ*ℙ ⋗ ℙ*𝕀 1.43×
                :ξy  => ℙ[+0, +0, +1, +0, -1, +1, +1, -1, -1],
            ),
        ),
        # Properties
        :pro => Dict{Symbol, ℙ}(
            :ν      => nu,
            :τ      => ℙ(3nu + 0.5),
            :u_max  => ℙ(0.04 / scale),
            :ρ₀     => one(ℙ),
        ),
    )
end


