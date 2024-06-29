#----------------------------------------------------------------------------------------------#
#                                 11-ORIG-julia-BaseVector.jl                                  #
#----------------------------------------------------------------------------------------------#
# A "from the ground up", pure, and plain julia  implementation  of  the  Taylor-Green  vortex #
# decay problem solver in LBM.                                                                 #
#----------------------------------------------------------------------------------------------#

"""
"""
function init(ℙ::Type{𝕋} where 𝕋<:AbstractFloat,    # The floating point precision
              l::Int)                               # log2(scale)
    𝕀, 𝕌 = ℙ == Float64 ? Int64, UInt64 : Int32, UInt32
    scale   = 𝕌(1) << l
    chunk   = 𝕌(32)
    NY = NX = scale * chunk
    w0, w1, w2 = ℙ(4.0/9.0), ℙ(1.0/9.0), ℙ(1.0/36.0)
    PAR = Dict{Symbol, Dict}(
        # Domain parameters
        :dom => Dict{Symbol, 𝕌}(
            :sca => scale,
            :NX  => NX,
            :NY  => NY,
        ),
        # Lattice stencil
        :lat => Dict{Symbol, Union{𝕌, ℙ, Vector{Union{𝕀, ℙ}}}}(
            :dim => 𝕌(2),   # D2...
            :vel => 𝕌(9),   # ...Q9
            :w   => ℙ[w0, w1, w1, w1, w1, w2, w2, w2, w2]
            :ξx  => 𝕀[+0, +1, +0, -1, +0, +1, -1, -1, +1]
            :ξy  => 𝕀[+0, +0, +1, +0, -1, +1, +1, -1, -1]
        ),
    )
