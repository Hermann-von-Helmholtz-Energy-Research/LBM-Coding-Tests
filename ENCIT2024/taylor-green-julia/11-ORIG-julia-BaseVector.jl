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
init(_::Type{𝕋} where 𝕋<:AbstractFloat, l::Int)::Dict{Symbol, Dict}
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
        cas=(sca=scale, NX=NX, NY=NY, IT=𝕀(round(maxIt / scale / scale))),
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
taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀, NY::𝕀;
             ν::𝕋, τ::𝕋, u_max::𝕋, ρ₀::𝕋)::NTuple{3, 𝕋} where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex decay.

```
> using BenchmarkTools;
> include("./11-ORIG-julia-BaseVector.jl");
> par = init(Float64, 0);
> t = par.typ.f(0.0); NX = par.cas.NX; NY = par.cas.NY; pro = par.pro;
> b = @benchmarkable taylor_green(t, 17, 17, NX, NY; \$pro...);
> b.params.evals = 1000;
> b.params.seconds = 25.0;
> run(b)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  67.514 ns …  3.383 μs  ┊ GC (min … max): 0.00% … 96.89%
 Time  (median):     69.805 ns (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   75.194 ns ± 89.412 ns  ┊ GC (mean ± σ):  4.74% ±  3.94%
     ↓
  ▆▄██▄▂                                                      ▁
  ███████▆▅▅▄▄▅▂▄▄▆▆▇███▇▇▆▅▄▄▅▅▅▅▄▃▄▅▅▇▇▇▇▅▅▅▅▆▄▄▅▅▆▂▅▆▇▅▆▆▆ █
  67.5 ns      Histogram: log(frequency) by time       114 ns <

 Memory estimate: 80 bytes, allocs estimate: 2.
```
"""
function taylor_green(t::𝕋, x::𝕀, y::𝕀, NX::𝕀, NY::𝕀;
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

"""
```
taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2}, NX::𝕀, NY::𝕀;
             pro::NamedTuple)::Nothing where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex in FIELDS (ρ, 𝑢, 𝑣).

```
> using BenchmarkTools;
> include("./11-ORIG-julia-BaseVector.jl");
> par = init(Float64, 0);
> t = par.typ.f(0.0); NX = par.cas.NX; NY = par.cas.NY; pro = par.pro;
> ρ = Array{Float64, 2}(undef, 32, 32);
> 𝑢 = Array{Float64, 2}(undef, 32, 32);
> 𝑣 = Array{Float64, 2}(undef, 32, 32);
> B = @benchmarkable taylor_green(t, ρ, 𝑢, 𝑣, NX, NY; pro=\$pro)
> tune!(B)
> run(B)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  35.970 μs … 78.185 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     36.195 μs (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   36.660 μs ±  3.068 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
   ↓
  █▇                                                          ▁
  ███▇▆▇▆▅▁█▇▆▅▅▄▃▄▁▃▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▅▆▇▇▅▅▅▄▅▆▇▇ █
  36 μs        Histogram: log(frequency) by time      51.4 μs <

 Memory estimate: 48 bytes, allocs estimate: 1.
```
"""
function taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2}, NX::𝕀, NY::𝕀;
                      pro::NamedTuple)::Nothing where {𝕋, 𝕀}
    for j in axes(ρ, 2)
        for i in axes(ρ, 1)
            ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, NX, NY; pro...)
        end
    end
end

"""
```
taylor_green_sq(t::𝕋, x::𝕀, y::𝕀;
                cas::Dict{Symbol, 𝕀},
                pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex decay in a square domain.

```
> using BenchmarkTools, Unitful
> include("./11-ORIG-julia-BaseVector.jl");
> par = init(Float64, 0);
> t = par.typ.f(0.0); NX = par.cas.NX; pro = par.pro;
> 𝑏 = @benchmarkable taylor_green_sq(t, 17, 17, NX; \$pro...);
> 𝑏.params.evals = 1000;
> 𝑏.params.seconds = 25.0;
> run(𝑏)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  64.539 ns …  3.647 μs  ┊ GC (min … max): 0.00% … 97.39%
 Time  (median):     66.546 ns (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   71.032 ns ± 85.556 ns  ┊ GC (mean ± σ):  4.69% ±  3.93%
     ↓
  ▅▇▇█▅▂                                                      ▂
  ████████▇▆▆▆▃▁▃▁▁▁▅▇█▇▇█▇▆▆▅▁▁▅▅▅▅▄▅▄▃▄▅▄▆▇▇▆▅▄▄▄▆▆▆▅▅▆▅▄▃▅ █
  64.5 ns      Histogram: log(frequency) by time       105 ns <

 Memory estimate: 80 bytes, allocs estimate: 2.
```
"""
function taylor_green_sq(t::𝕋, x::𝕀, y::𝕀, NX::𝕀;
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

"""
```
taylor_green_sq(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2}, NX::𝕀;
                pro::NamedTuple)::Nothing where {𝕋, 𝕀}
```
Function to compute the exact solution for Taylor-Green vortex in square FIELDS (ρ, 𝑢, 𝑣).

```
> using BenchmarkTools;
> include("./11-ORIG-julia-BaseVector.jl");
> par = init(Float64, 0);
> t = par.typ.f(0.0); NX = par.cas.NX; pro = par.pro;
> ρ = Array{Float64, 2}(undef, 32, 32);
> 𝑢 = Array{Float64, 2}(undef, 32, 32);
> 𝑣 = Array{Float64, 2}(undef, 32, 32);
> 𝐵 = @benchmarkable taylor_green_sq(t, ρ, 𝑢, 𝑣, NX; pro=\$pro)
> tune!(𝐵)
> run(𝐵)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  34.185 μs … 79.878 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     34.293 μs (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   34.738 μs ±  3.013 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
  ↓
  █▁                                                          ▁
  ██▇▆▇▆▆▄▆▇▆▅▄▄▃▄▃▃▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▄▅▆▆▆▅▅▃▅▆▆▆▆ █
  34.2 μs      Histogram: log(frequency) by time      49.6 μs <

 Memory estimate: 48 bytes, allocs estimate: 1.
```
"""
function taylor_green_sq(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2}, NX::𝕀;
                         pro::NamedTuple)::Nothing where {𝕋, 𝕀}
    for j in axes(ρ, 2)
        for i in axes(ρ, 1)
            ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green_sq(t, i, j, NX; pro...)
        end
    end
end


#----------------------------------------------------------------------------------------------#
#                    Equilibrium Mass Distribution Function Initialization                     #
#----------------------------------------------------------------------------------------------#

"""
```
init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                 w::Vector{𝕋}, ξx::Vector{𝕋}, ξy::Vector{𝕋})::Nothing where 𝕋
```
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.

```
> include("./11-ORIG-julia-BaseVector.jl");
> using BenchmarkTools
> par = init(Float64, 0);
> 𝕀, 𝕋 = par.typ; t = 𝕋(0.0);
> sca, NX, NY, IT = par.cas; pro = par.pro;
> f = Array{𝕋, 3}(undef, NX, NY, par.lat.int.vel);
> ρ = Array{𝕋, 2}(undef, NX, NY);
> 𝑢 = Array{𝕋, 2}(undef, NX, NY);
> 𝑣 = Array{𝕋, 2}(undef, NX, NY);
> taylor_green(t, ρ, 𝑢, 𝑣, NX, NY; pro=pro)
> vec = par.lat.vec
> 𝐃 = @benchmarkable init_equilibrium(f, ρ, 𝑢, 𝑣, par.lat.vec...);
> tune!(𝐃);
> run(𝐃)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  23.356 μs …  60.150 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     23.800 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   23.842 μs ± 705.056 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
                           ↓
                   ▁▃▄▅▆▆▇█▇▆▆▅▃▂▁
  ▂▁▁▂▂▂▂▂▃▃▃▄▄▅▆▇█████████████████▇▆▅▄▄▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▄
  23.4 μs         Histogram: frequency by time         24.4 μs <

 Memory estimate: 352 bytes, allocs estimate: 5.
```
"""
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2},
                          w::Vector{𝕋}, ξx::Vector{𝕋}, ξy::Vector{𝕋})::Nothing where 𝕋
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
init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2};
                 latvec::Dict{Symbol, Vector{𝕋}})::Nothing where 𝕋
```
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.

```
[...]
> 𝐄 = @benchmarkable init_equilibrium(f, ρ, 𝑢, 𝑣; latvec=Dict(pairs(par.lat.vec)));
> tune!(𝐄);
> run(𝐄)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  19.965 μs …  75.489 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     20.410 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   20.485 μs ± 997.724 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
                        ↓
               ▁▃▄▅▆▇▇▇█▇▆▆▄▃▂▁
  ▁▁▁▁▂▂▂▃▄▄▅▇██████████████████▆▆▅▅▅▄▄▃▃▃▃▃▃▃▃▃▂▂▂▂▂▂▂▁▁▁▁▁▁▁ ▄
  20 μs           Histogram: frequency by time         21.2 μs <

 Memory estimate: 912 bytes, allocs estimate: 11.
```
"""
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2};
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

# Note: To me, it makes absolutely NO SENSE on why the second impl of init_equilibrium(), with
# two extra function calls, namely, (i) Dict() and (ii) pairs(), and three extra variable
# assignments, namely, (a) ξx, (b) ξy, and (c) w; can be FASTER than the first implementation of
# it, which is based on NamedTuples, which have been measured 100's of times faster than Dict()
# to generate... really puzzled here!


