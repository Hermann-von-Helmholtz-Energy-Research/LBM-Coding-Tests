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
        ▄▆▇██▇▆▅▄▄▂▁▁
  ▂▃▃▅▆██████████████▇▇▅▄▄▃▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
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
BenchmarkTools.Trial: 10000 samples with 200 evaluations.
 Range (min … max):  403.895 ns …  10.905 μs  ┊ GC (min … max): 0.00% … 95.08%
 Time  (median):     405.905 ns (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   409.369 ns ± 151.088 ns  ┊ GC (mean ± σ):  0.64% ±  1.81%
          ↓
    ▁▃▅▇███▇▆▅▄▂▁▂▂▂▁▂▂▁▁                               ▁       ▃
  ▅▇███████████████████████▇▆▆▅▁▄▄▄▁▃▁▁▁▁▁▁▁▁▃▃▁▁▁▁▄▃▅▇▇██▆██▇▇ █
  404 ns        Histogram: log(frequency) by time        419 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.
```
"""
function taylor_green(t::𝕋, x::𝕌, y::𝕌;
                      cas::Dict{Symbol, 𝕌},
                      pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}
    𝐍𝐱  = cas[:NX]
    𝐍𝐲  = cas[:NY]
    ϱ   = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐𝛑  = 𝟐 * π
    kx  = 𝟐𝛑 / 𝐍𝐱       # promote_type(UInt32, Float##) -> Float##
    ky  = 𝟐𝛑 / 𝐍𝐲
    td  = pro[:ν] * (kx*kx + ky*ky)
    𝐔𝐞  = pro[:u_max] * exp(-t * td)
    X   = x - 𝐍𝐱 / 𝟐    # Centered vortex
    Y   = y - 𝐍𝐲 / 𝟐    # Centered vortex
    sx, cx  = sincos(kx * X)
    sy, cy  = sincos(ky * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * √(ky / kx) * cx * sy
    𝚟   = + 𝐔𝐞 * √(kx / ky) * sx * cy
    P   = - 𝕋(0.25) * ϱ * 𝐔𝐞 * 𝐔𝐞 * ((ky / kx) * c2x + (kx / ky) * c2y)
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
 Range (min … max):  402.195 ns …  11.641 μs  ┊ GC (min … max): 0.00% … 95.25%
 Time  (median):     404.585 ns (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   408.274 ns ± 162.089 ns  ┊ GC (mean ± σ):  0.69% ±  1.83%
           ↓
      ▃▅▇███▇▆▆▅▄▃▃▂▂▁ ▁▁▁▁                              ▁      ▃
  ▄▆▆█████████████████████████▇▆▇▅▆▅▃▄▁▁▅▁▄▁▁▃▃▁▃▃▅▅▆▇▇█▇█▇▇▇▇▇ █
  402 ns        Histogram: log(frequency) by time        418 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.

```
"""
function taylor_green_sq(t::𝕋, x::𝕌, y::𝕌;
                         cas::Dict{Symbol, 𝕌},
                         pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕌}
    𝐍   = cas[:NX]
    ϱ   = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐𝛑  = 𝟐 * π
    k   = 𝟐𝛑 / 𝐍        # promote_type(UInt32, Float##) -> Float##
    td  = pro[:ν] * 𝟐 * k * k   # 𝟐*k*k ⋗ (k*k+k*k) 2.81x
    𝐔𝐞  = pro[:u_max] * exp(-t * td)
    X   = x - 𝐍 / 𝟐     # Centered vortex
    Y   = y - 𝐍 / 𝟐     # Centered vortex
    sx, cx  = sincos(k * X)
    sy, cy  = sincos(k * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * cx * sy
    𝚟   = + 𝐔𝐞 * sx * cy
    P   = - 𝕋(0.25) * ϱ * 𝐔𝐞 * 𝐔𝐞 * (c2x + c2y)
    ρ   = ϱ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

"""
```
init_equilibrium()
```\n
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.
"""
function init_equilibrium(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋};
                          cas::Dict{Symbol, 𝕌},
                          lat::Dict{Symbol, Dict})::Nothing where {𝕋, 𝕌}
    vec = lat[:vec]
    ξx  = vec[:ξx]
    ξy  = vec[:ξy]
    for 𝑦 in Base.OneTo(cas[:NY])
        for 𝑥 in Base.OneTo(cas[:NX])
            ϱ, 𝚞, 𝚟 = ρ[𝑥, 𝑦], 𝑢[𝑥, 𝑦], 𝑣[𝑥, 𝑦]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟                      # OP1
            for 𝑖 in Base.OneTo(lat[:int][:vel])
                ξ𝘂 = ξx[𝑖] * 𝚞 + ξy[𝑖] * 𝚟
                𝑓[𝑥, 𝑦, 𝑖] = lat[:vec][:wi][𝑖] * ϱ * (
                    + 𝕋(1.0)
                    + 𝕋(3.0) * ξ𝘂
                    + 𝕋(4.5) * ξ𝘂 * ξ𝘂
                    - 𝕋(1.5) * 𝘂𝘂
                )
            end
        end
    end
end


