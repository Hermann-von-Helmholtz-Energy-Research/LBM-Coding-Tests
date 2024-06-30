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

```julia-REPL
julia> @benchmark init(Float64, 0)
BenchmarkTools.Trial: 10000 samples with 10 evaluations.
 Range (min … max):  1.478 μs … 323.095 μs  ┊ GC (min … max):  0.00% … 97.87%
 Time  (median):     1.558 μs (↓)           ┊ GC (median):     0.00%
 Time  (mean ± σ):   1.798 μs ±   6.372 μs  ┊ GC (mean ± σ):  10.97% ±  3.24%
              ↓
        ▄▆▇██▇┊▅▄▄▂▁▁
  ▂▃▃▅▆███████┊██████▇▇▅▄▄▃▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  1.48 μs         Histogram: frequency by time        1.87 μs <

 Memory estimate: 4.91 KiB, allocs estimate: 50.
```
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
            :flo => Dict{Symbol, ℙ}(:a => √ℙ(3.0/2.0), :cs => inv(√ℙ(3.0))),
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

"""
```
taylor_green(t::𝕋, x::𝕌, y::𝕌;
             cas::Dict{Symbol, 𝕌},
             pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}```\n
Function to compute the exact solution for Taylor-Green vortex decay.

```julia-REPL
julia> par = init(Float64, 0)
[...]
julia> @benchmark taylor_green(par[:typ][:p](0.0), UInt64(17), UInt64(17), cas=par[:cas], pro=par[:pro])
BenchmarkTools.Trial: 10000 samples with 197 evaluations.
 Range (min … max):  452.980 ns …  13.235 μs  ┊ GC (min … max): 0.00% … 95.39%
 Time  (median):     455.132 ns (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   460.005 ns ± 188.156 ns  ┊ GC (mean ± σ):  0.74% ±  1.85%
         ↓
      ▁▅██▄                
  ▂▂▃▅██████▆▅▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▂▁▁▁▁▂▁▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▃
  453 ns           Histogram: frequency by time          471 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.

```
"""
function taylor_green(t::𝕋, x::𝕌, y::𝕌;
                      cas::Dict{Symbol, 𝕌},
                      pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}
    𝐍𝐱  = cas[:NX]
    𝐍𝐲  = cas[:NY]
    𝐔   = pro[:u_max]
    ϱ   = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐𝛑  = 𝟐 * π
    kx  = 𝟐𝛑 / 𝐍𝐱     # promote_type(UInt32, Float##) -> Float##
    ky  = 𝟐𝛑 / 𝐍𝐲
    td  = pro[:ν] * (kx*kx + ky*ky)
    X   = 𝕋(x) - 𝐍𝐱 / 𝟐     # Centered vortex
    Y   = 𝕋(y) - 𝐍𝐲 / 𝟐     # Centered vortex
    𝚞   = - 𝐔 * √(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-t * td)
    𝚟   = + 𝐔 * √(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-t * td)
    P   = - 𝕋(0.25) * ϱ * 𝐔 * 𝐔 * ((ky / kx) * cos(𝟐 * kx * X) + (kx / ky) * sin(𝟐 * ky * Y))
    ρ   = ϱ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

"""
```
taylor_green_sq(t::𝕋, x::𝕌, y::𝕌;
                cas::Dict{Symbol, 𝕌},
                pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}```\n
Function to compute the exact solution for Taylor-Green vortex decay in a square domain.

```julia-REPL
julia> par = init(Float64, 0)
[...]
julia> @benchmark taylor_green_sq(par[:typ][:p](0.0), UInt64(17), UInt64(17), cas=par[:cas], pro=par[:pro])
BenchmarkTools.Trial: 10000 samples with 200 evaluations.
 Range (min … max):  404.690 ns …  11.110 μs  ┊ GC (min … max): 0.00% … 94.99%
 Time  (median):     406.475 ns (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   410.443 ns ± 154.835 ns  ┊ GC (mean ± σ):  0.66% ±  1.81%
         ↓
    ▃▅▇██▇▆▅▄▂▂▂▂▂▂▂▁▁                                          ▂
  ▆████████████████████▇▆▆▅▆▃▄▃▃▃▅▁▁▁▁▄▁▁▁▃▁▁▁▁▃▄▇█▇███▇▇▆▇▄▅▅▆ █
  405 ns        Histogram: log(frequency) by time        421 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.

```
"""
function taylor_green_sq(t::𝕋, x::𝕌, y::𝕌;
                         cas::Dict{Symbol, 𝕌},
                         pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}
    𝐍   = cas[:NX]
    𝐔   = pro[:u_max]
    ϱ   = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐𝛑  = 𝟐 * π
    k   = 𝟐𝛑 / 𝐍                # promote_type(UInt32, Float##) -> Float##
    td  = pro[:ν] * k * k * 𝟐
    X   = 𝕋(x) - 𝐍 / 𝟐          # Centered vortex
    Y   = 𝕋(y) - 𝐍 / 𝟐          # Centered vortex
    𝚞   = - 𝐔 * cos(k * X) * sin(k * Y) * exp(-t * td)
    𝚟   = + 𝐔 * sin(k * X) * cos(k * Y) * exp(-t * td)
    P   = - 𝕋(0.25) * ϱ * 𝐔 * 𝐔 * (cos(𝟐 * k * X) + sin(𝟐 * k * Y))
    ρ   = ϱ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end


