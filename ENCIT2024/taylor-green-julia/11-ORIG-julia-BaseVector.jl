#==============================================================================================#
#                                 11=ORIG=julia=BaseVector.jl                                  #
#==============================================================================================#
# A "from the ground up", pure, and plain julia  implementation  of  the  Taylor-Green  vortex #
# decay problem solver in LBM.                                                                 #
#==============================================================================================#

#----------------------------------------------------------------------------------------------#
#                                  Simulation Initialization                                   #
#----------------------------------------------------------------------------------------------#

"""
```
function init(_::Type{𝕋}, l::Int)::NamedTuple where 𝕋<:AbstractFloat
```
Computes (i) types, (ii) case, (iii) lattice, and (iv) properties simulation parameters, and
returns as a `Dict{Symbol, Dict}`.

```
> using BenchmarkTools
> include("./11-ORIG-julia-BaseVector.jl");
> @benchmark init(Float64, 0)
BenchmarkTools.Trial: 10000 samples with 949 evaluations.
 Range (min … max):   97.034 ns …  1.807 μs  ┊ GC (min … max): 0.00% … 92.40%
 Time  (median):     106.483 ns (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   108.382 ns ± 26.571 ns  ┊ GC (mean ± σ):  1.32% ±  5.03%
              ↓
          ▁▄▇██▇▄▂▁▁                                           ▂
  ▆▅▁▄▇▆▇▇████████████▇▆▄▃▁▄▄▅▅▆▅▄▁▁▄▄▁▅▄▄▄▁▄▄▁▁▃▅▄▁▄▃▄▃▁▄▁▆▃▄ █
  97 ns         Histogram: log(frequency) by time       145 ns <

 Memory estimate: 384 bytes, allocs estimate: 3.
```
"""
function init(_::Type{𝕋}, l::Int)::NamedTuple where 𝕋<:AbstractFloat
    𝕀           = 𝕋 == Float64 ? Int64 : Int32
    scale       = 𝕀(1) << l
    chunk       = 𝕀(32)
    maxIt       = 𝕀(204800)
    NY = NX     = scale * chunk
    nu          = 𝕋(1.0/6.0)
    w0, w1, w2  = 𝕋(4.0/9.0), 𝕋(1.0/9.0), 𝕋(1.0/36.0)
    return (
        typ=(i=𝕀, f=𝕋),
        sup=(sca=scale, IT=𝕀(round(maxIt / scale / scale))),
        cas=(NX=NX, NY=NY),
        lat=(int=(dim=𝕀(2), vel=𝕀(9)), flo=(a=√𝕋(3.0/2.0), cs=inv(√𝕋(3.0))),
            vec=(w=𝕋[w0, w1, w1, w1, w1, w2, w2, w2, w2],
                 ξx=𝕋[+0, +1, +0, -1, +0, +1, -1, -1, +1],  # 𝕋*𝕋 ⋗ 𝕋*𝕀 1.43×
                 ξy=𝕋[+0, +0, +1, +0, -1, +1, +1, -1, -1])),
        pro=(ν=nu, τ=𝕋(3.0nu + 0.5), u_max=𝕋(0.04 / scale), ρ₀=one(𝕋))
    )
end


#----------------------------------------------------------------------------------------------#
#                                         Taylor-Green                                         #
#----------------------------------------------------------------------------------------------#

"""
```
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀, NY::𝕀,
                      ν::𝕋, τ::𝕋, u_max::𝕋, ρ₀::𝕋)::NTuple{3, 𝕋} where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex decay in a rectangular domain.
"""
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀, NY::𝕀,
                      ν::𝕋, τ::𝕋, u_max::𝕋, ρ₀::𝕋)::NTuple{3, 𝕋} where {𝕋, 𝕀}
    𝟐   = 𝕋(2.0)
    𝟐⁻¹ = inv(𝟐)
    𝟐𝛑  = 𝟐 * π
    kx  = 𝟐𝛑 / NX       # promote_type(UInt32, Float##) -> Float##
    ky  = 𝟐𝛑 / NY
    td  = ν * (kx*kx + ky*ky)
    𝐔𝐞  = u_max * exp(-t * td)
    X   = x + 𝟐⁻¹
    Y   = y + 𝟐⁻¹
    sx, cx  = sincos(kx * X)
    sy, cy  = sincos(ky * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * √(ky / kx) * cx * sy
    𝚟   = + 𝐔𝐞 * √(kx / ky) * sx * cy
    P   = - 𝕋(0.25) * ρ₀ * 𝐔𝐞 * 𝐔𝐞 * ((ky / kx) * c2x + (kx / ky) * c2y)
    ρ   = ρ₀ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

# Dictionary argument version
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀, NY::𝕀,
                      pro::Dict{Symbol,𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕀}
    ν       = pro[:ν]
    τ       = pro[:τ]
    u_max   = pro[:u_max]
    ρ₀      = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐⁻¹ = inv(𝟐)
    𝟐𝛑  = 𝟐 * π
    kx  = 𝟐𝛑 / NX       # promote_type(UInt32, Float##) -> Float##
    ky  = 𝟐𝛑 / NY
    td  = ν * (kx*kx + ky*ky)
    𝐔𝐞  = u_max * exp(-t * td)
    X   = x + 𝟐⁻¹
    Y   = y + 𝟐⁻¹
    sx, cx  = sincos(kx * X)
    sy, cy  = sincos(ky * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * √(ky / kx) * cx * sy
    𝚟   = + 𝐔𝐞 * √(kx / ky) * sx * cy
    P   = - 𝕋(0.25) * ρ₀ * 𝐔𝐞 * 𝐔𝐞 * ((ky / kx) * c2x + (kx / ky) * c2y)
    ρ   = ρ₀ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

"""
```
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀,
                      ν::𝕋, τ::𝕋, u_max::𝕋, ρ₀::𝕋)::NTuple{3, 𝕋} where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex decay in a square domain.
"""
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀,
                      ν::𝕋, τ::𝕋, u_max::𝕋, ρ₀::𝕋)::NTuple{3, 𝕋} where {𝕋, 𝕀}
    𝟐   = 𝕋(2.0)
    𝟐⁻¹ = inv(𝟐)
    𝟐𝛑  = 𝟐 * π
    k   = 𝟐𝛑 / NX       # promote_type(UInt32, Float##) -> Float##
    td  = ν * (k*k + k*k)
    𝐔𝐞  = u_max * exp(-t * td)
    X   = x + 𝟐⁻¹
    Y   = y + 𝟐⁻¹
    sx, cx  = sincos(k * X)
    sy, cy  = sincos(k * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * cx * sy
    𝚟   = + 𝐔𝐞 * sx * cy
    P   = - 𝕋(0.25) * ρ₀ * 𝐔𝐞 * 𝐔𝐞 * (c2x + c2y)
    ρ   = ρ₀ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

# Dictionary argument version
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀,
                      pro::Dict{Symbol,𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕀}
    ν       = pro[:ν]
    τ       = pro[:τ]
    u_max   = pro[:u_max]
    ρ₀      = pro[:ρ₀]
    𝟐   = 𝕋(2.0)
    𝟐⁻¹ = inv(𝟐)
    𝟐𝛑  = 𝟐 * π
    k   = 𝟐𝛑 / NX       # promote_type(UInt32, Float##) -> Float##
    td  = ν * (k*k + k*k)
    𝐔𝐞  = u_max * exp(-t * td)
    X   = x + 𝟐⁻¹
    Y   = y + 𝟐⁻¹
    sx, cx  = sincos(k * X)
    sy, cy  = sincos(k * Y)
    c2x = cx * cx - sx * sx
    c2y = cy * cy - sy * sy
    𝚞   = - 𝐔𝐞 * cx * sy
    𝚟   = + 𝐔𝐞 * sx * cy
    P   = - 𝕋(0.25) * ρ₀ * 𝐔𝐞 * 𝐔𝐞 * (c2x + c2y)
    ρ   = ρ₀ + 𝕋(3.0) * P
    return ρ, 𝚞, 𝚟
end

"""
```
function taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                      pro::NamedTuple)::Nothing where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex in FIELDS (ρ, 𝑢, 𝑣).

```
> include("./11-ORIG-julia-BaseVector.jl");
> using BenchmarkTools
> par = init(Float64, 0);
> 𝕀, 𝕋 = par.typ;
> NX, NY = par.cas;
> f = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> g = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> ρ = Array{𝕋, 2}(undef, NX, NY);
> 𝑢 = Array{𝕋, 2}(undef, NX, NY);
> 𝑣 = Array{𝕋, 2}(undef, NX, NY);
> @benchmark taylor_green(zero(𝕋), ρ, 𝑢, 𝑣, par.pro)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  35.965 μs … 214.001 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     36.177 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   36.276 μs ±   2.113 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
       ↓
     ▂█▇▃▁ 
  ▂▄▆█████▄▂▂▂▂▂▂▁▂▁▂▂▁▁▁▂▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▁▂▂▂▂▂ ▃
  36 μs           Histogram: frequency by time         38.4 μs <

 Memory estimate: 192 bytes, allocs estimate: 2.
```
"""
function taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                      pro::NamedTuple)::Nothing where 𝕋
    NX, NY = size(ρ)
    if NX == NY
        for j in axes(ρ, 2)
            for i in axes(ρ, 1)
                ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, NX, pro...)
            end
        end
    else
        for j in axes(ρ, 2)
            for i in axes(ρ, 1)
                ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, NX, NY, pro...)
            end
        end
    end
end

# Dictionary version
"""
```
> [...]
> @benchmark taylor_green(zero(𝕋), ρ, 𝑢, 𝑣, Dict(pairs(par.pro)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  66.482 μs … 107.583 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     66.794 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   66.925 μs ±   1.147 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
         ↓
      ▄██▅▂  
  ▂▂▃▇█████▇▆▄▃▃▃▃▃▃▂▃▂▂▂▂▂▂▂▂▁▂▁▂▂▂▂▁▁▁▁▁▁▁▂▁▁▂▁▁▁▁▁▁▂▂▂▂▂▂▂▂ ▃
  66.5 μs         Histogram: frequency by time         69.2 μs <

 Memory estimate: 800 bytes, allocs estimate: 7.
```
"""
function taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                      pro::Dict{Symbol,𝕋})::Nothing where 𝕋
    NX, NY = size(ρ)
    if NX == NY
        for j in axes(ρ, 2)
            for i in axes(ρ, 1)
                ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, NX, pro)
            end
        end
    else
        for j in axes(ρ, 2)
            for i in axes(ρ, 1)
                ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, NX, NY, pro)
            end
        end
    end
end

# At least, the above benchmark results make sense!


#----------------------------------------------------------------------------------------------#
#                    Equilibrium Mass Distribution Function Initialization                     #
#----------------------------------------------------------------------------------------------#

"""
```
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                          vec::NamedTuple)::Nothing where 𝕋
```
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.

```
> include("./11-ORIG-julia-BaseVector.jl");
> using BenchmarkTools
> par = init(Float64, 0);
> 𝕀, 𝕋 = par.typ;
> NX, NY = par.cas;
> f = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> g = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> ρ = Array{𝕋, 2}(undef, NX, NY);
> 𝑢 = Array{𝕋, 2}(undef, NX, NY);
> 𝑣 = Array{𝕋, 2}(undef, NX, NY);
> taylor_green(zero(𝕋), ρ, 𝑢, 𝑣, par.pro)
> 𝐃 = @benchmarkable init_equilibrium(f, ρ, 𝑢, 𝑣, par.lat.vec);
> tune!(𝐃);
> run(𝐃)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  22.803 μs … 92.733 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     23.194 μs (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   23.246 μs ±  1.011 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
                        ↓
                 ▂▃▅▆▇▆█▇▇█▅▅▂▂
  ▁▁▁▁▁▂▂▂▃▃▄▅▆▇█████████████████▆▆▄▄▃▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  22.8 μs         Histogram: frequency by time        23.8 μs <

 Memory estimate: 304 bytes, allocs estimate: 4.
```
"""
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                          vec::NamedTuple)::Nothing where 𝕋
    w, ξx, ξy = vec
    for 𝑦 in axes(𝑓, 2)
        for 𝑥 in axes(𝑓, 1)
            ϱ, 𝚞, 𝚟 = ρ[𝑥, 𝑦], 𝑢[𝑥, 𝑦], 𝑣[𝑥, 𝑦]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟                      # OP1
            for 𝑖 in axes(𝑓, 3)
                ξ𝘂 = ξx[𝑖] * 𝚞 + ξy[𝑖] * 𝚟
                𝑓[𝑥, 𝑦, 𝑖] = w[𝑖] * ϱ * (
                    + 𝕋(1.0)
                    + 𝕋(3.0) * ξ𝘂
                    + 𝕋(4.5) * ξ𝘂 * ξ𝘂
                    - 𝕋(1.5) * 𝘂𝘂
                )
            end
        end
    end
end

"""
```
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                          latvec::Dict{Symbol, Vector{𝕋}})::Nothing where 𝕋
```
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.

```
[...]
> 𝐄 = @benchmarkable init_equilibrium(f, ρ, 𝑢, 𝑣, Dict(pairs(par.lat.vec)));
> tune!(𝐄);
> run(𝐄)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  19.749 μs …  74.642 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     20.240 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   20.305 μs ± 899.474 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
                         ↓
                ▁▃▄▅▇▇███▇▆▆▅▃▃▂
  ▂▁▁▂▂▂▂▃▄▄▅▅▇███████████████████▇▇▅▅▄▄▄▄▄▄▃▃▃▃▃▃▃▃▃▃▃▂▂▂▂▂▂▂ ▅
  19.7 μs         Histogram: frequency by time           21 μs <

 Memory estimate: 880 bytes, allocs estimate: 9.
```
"""
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                          latvec::Dict{Symbol, Vector{𝕋}})::Nothing where 𝕋
    ξx = latvec[:ξx]
    ξy = latvec[:ξy]
    w  = latvec[:w]
    for 𝑦 in axes(𝑓, 2)
        for 𝑥 in axes(𝑓, 1)
            ϱ, 𝚞, 𝚟 = ρ[𝑥, 𝑦], 𝑢[𝑥, 𝑦], 𝑣[𝑥, 𝑦]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟                      # OP1
            for 𝑖 in axes(𝑓, 3)
                ξ𝘂 = ξx[𝑖] * 𝚞 + ξy[𝑖] * 𝚟
                𝑓[𝑥, 𝑦, 𝑖] = w[𝑖] * ϱ * (
                    + 𝕋(1.0)
                    + 𝕋(3.0) * ξ𝘂
                    + 𝕋(4.5) * ξ𝘂 * ξ𝘂
                    - 𝕋(1.5) * 𝘂𝘂
                )
            end
        end
    end
end

# Note: To me, it makes NO SENSE on why the second implementation of init_equilibrium(), with
# two extra function calls, namely, (i) Dict() and (ii) pairs(), and three extra variable
# assignments, namely, (a) ξx, (b) ξy, and (c) w; can be FASTER than the first implementation of
# it. Perhaps the passing of Vector{}'s to a function in a call is very costly.


#----------------------------------------------------------------------------------------------#
#                                          Streaming                                           #
#----------------------------------------------------------------------------------------------#

"""
```
function stream(𝑓::Array{𝕋, 3}, 𝑔::Array{𝕋, 3},
                vec::NamedTuple, _::Type{𝕀})::Nothing where {𝕋, 𝕀 <: Integer}
```
Function that performs streaming of the populations in a fully periodic domain, reading from 𝑓
and writing to 𝑔.

```
> include("./11-ORIG-julia-BaseVector.jl");
> using BenchmarkTools
> par = init(Float64, 0);
> 𝕀, 𝕋 = par.typ;
> NX, NY = par.cas;
> f = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> g = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> ρ = Array{𝕋, 2}(undef, NX, NY);
> 𝑢 = Array{𝕋, 2}(undef, NX, NY);
> 𝑣 = Array{𝕋, 2}(undef, NX, NY);
> taylor_green(zero(𝕋), ρ, 𝑢, 𝑣, par.pro)
> init_equilibrium(f, ρ, 𝑢, 𝑣, par.lat.vec);
> @benchmark stream(f, g, par.lat.vec, par.typ.i)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  57.821 μs … 131.787 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     57.881 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   58.005 μs ±   1.556 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
   ↓
  ██▇▅▃▂▁▁▁ ▁                                                  ▂
  ███████████████▇▇▆▅▁▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆▆▆▇▇█ █
  57.8 μs       Histogram: log(frequency) by time      60.3 μs <

 Memory estimate: 736 bytes, allocs estimate: 8.
```
"""
function stream(𝑓::Array{𝕋, 3}, 𝑔::Array{𝕋, 3},
                vec::NamedTuple, _::Type{𝕀})::Nothing where {𝕋, 𝕀 <: Integer}
    ξx = Vector{𝕀}(vec.ξx)
    ξy = Vector{𝕀}(vec.ξy)
    NX, NY, NL = size(𝑓)
    for 𝑦 in axes(𝑓, 2)
        for 𝑥 in axes(𝑓, 1)
            for 𝑖 in axes(𝑓, 3)
                # "from" indices, enforcing periodicity
                𝑝 = (NX + 𝑥 - ξx[𝑖]) % NX + 1       # NX is added as to guarantee positivity
                𝑞 = (NY + 𝑦 - ξy[𝑖]) % NY + 1       # NY is added as to guarantee positivity
                # Streaming from 𝑓 into 𝑔
                𝑔[𝑥, 𝑦, 𝑖] = 𝑓[𝑝, 𝑞, 𝑖]
            end
        end
    end
end

# Dictionary version
function stream(𝑓::Array{𝕋, 3}, 𝑔::Array{𝕋, 3},
                vec::Dict{Symbol, Vector{𝕋}}, _::Type{𝕀})::Nothing where {𝕋, 𝕀 <: Integer}
    ξx = Vector{𝕀}(vec[:ξx])
    ξy = Vector{𝕀}(vec[:ξy])
    NX, NY, NL = size(𝑓)
    for 𝑦 in axes(𝑓, 2)
        for 𝑥 in axes(𝑓, 1)
            for 𝑖 in axes(𝑓, 3)
                # "from" indices, enforcing periodicity
                𝑝 = (NX + 𝑥 - ξx[𝑖]) % NX + 1       # NX is added as to guarantee positivity
                𝑞 = (NY + 𝑦 - ξy[𝑖]) % NY + 1       # NY is added as to guarantee positivity
                # Streaming from 𝑓 into 𝑔
                𝑔[𝑥, 𝑦, 𝑖] = 𝑓[𝑝, 𝑞, 𝑖]
            end
        end
    end
end





