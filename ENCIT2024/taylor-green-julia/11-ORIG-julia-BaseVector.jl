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
`init(_::Type{𝕋} where 𝕋<:AbstractFloat, l::Int)::Dict{Symbol, Dict}`\n
Computes (i) types, (ii) case, (iii) lattice, and (iv) properties simulation parameters, and
returns as a `Dict{Symbol, Dict}`.

```julia-REPL
julia> using BenchmarkTools
julia> include("./11-ORIG-julia-BaseVector.jl");
julia> @benchmark init(Float64, 0)
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
taylor_green(t::𝕋, x::𝕀, y::𝕀;
             cas::NamedTuple,
             pro::NamedTuple)::NTuple{3, 𝕋} where {𝕋, 𝕀}```\n
Function to compute the exact solution for Taylor-Green vortex decay.

```julia-REPL
julia> using BenchmarkTools
julia> include("./11-ORIG-julia-BaseVector.jl");
julia> par = init(Float64, 0);
julia> b = @benchmarkable taylor_green(par.typ.f(0.0), 17, 17, cas=par.cas, pro=par.pro)
julia> b.params.evals = 1000;
julia> b.params.seconds = 25.0;
julia> run(b)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  402.291 ns …  5.000 μs  ┊ GC (min … max): 0.00% … 90.38%
 Time  (median):     406.224 ns (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   408.774 ns ± 61.312 ns  ┊ GC (mean ± σ):  0.42% ±  2.68%
        ↓
       █▇    
  ▂▂▃▅▆██▃▃▃▂▂▂▂▁▂▁▂▂▂▂▂▂▂▂▁▂▂▁▁▁▁▂▂▁▂▂▂▂▂▁▂▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▂
  402 ns          Histogram: frequency by time          442 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.
```
"""
function taylor_green(t::𝕋, x::𝕀, y::𝕀;
                      cas::NamedTuple{(:sca, :NX, :NY, :IT)},
                      pro::NamedTuple{(:ν, :τ, :u_max, :ρ₀)})::NTuple{3, 𝕋} where {𝕋, 𝕀}
    𝐍𝐱  = cas.NX
    𝐍𝐲  = cas.NY
    ϱ   = pro.ρ₀
    𝟐   = 𝕋(2.0)
    𝟐𝛑  = 𝟐 * π
    kx  = 𝟐𝛑 / 𝐍𝐱       # promote_type(UInt32, Float##) -> Float##
    ky  = 𝟐𝛑 / 𝐍𝐲
    td  = pro.ν * (kx*kx + ky*ky)
    𝐔𝐞  = pro.u_max * exp(-t * td)
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

function taylor_green(t::𝕋, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2};
                      cas::Dict{Symbol, 𝕀},
                      pro::Dict{Symbol, 𝕋})::Nothing where {𝕋, 𝕀}
    for j in axes(ρ, 2)
        for i in axes(ρ, 1)
            ρ[i, j], 𝑢[i, j], 𝑣[i, j] = taylor_green(t, i, j, cas=cas, pro=pro)
        end
    end
end

"""
```
taylor_green_sq(t::𝕋, x::𝕀, y::𝕀;
                cas::Dict{Symbol, 𝕀},
                pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕀}```\n
Function to compute the exact solution for Taylor-Green vortex decay in a square domain.

```julia-REPL
julia> using BenchmarkTools, Unitful
julia> include("./11-ORIG-julia-BaseVector.jl");
julia> par = init(Float64, 0);
julia> 𝑏 = @benchmarkable taylor_green_sq(par[:typ][:p](0.0), UInt64(17), UInt64(17), cas=par[:cas], pro=par[:pro]);
julia> 𝑏.params.evals = 1000;
julia> 𝑏.params.seconds = 25.0;
julia> run(𝑏)
BenchmarkTools.Trial: 10000 samples with 200 evaluations.
 Range (min … max):  399.807 ns …  2.559 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     405.780 ns (↓)          ┊ GC (median):    0.00%
 Time  (mean ± σ):   407.502 ns ± 45.823 ns  ┊ GC (mean ± σ):  0.39% ± 2.95%
             ↓
             █   
  ▂▃▃▃▄█▆▄▄▄██▄▃▂▂▂▂▂▂▁▂▂▂▁▁▂▂▂▂▂▂▂▁▁▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▂▂▂▂▂▂▂▂▂ ▂
  400 ns          Histogram: frequency by time          434 ns <

 Memory estimate: 96 bytes, allocs estimate: 3.
```
"""
function taylor_green_sq(t::𝕋, x::𝕀, y::𝕀;
                         cas::Dict{Symbol, 𝕀},
                         pro::Dict{Symbol, 𝕋})::NTuple{3, 𝕋} where {𝕋, 𝕀}
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


#----------------------------------------------------------------------------------------------#
#                    Equilibrium Mass Distribution Function Initialization                     #
#----------------------------------------------------------------------------------------------#

"""
```
init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2};
                 latvec::Dict{Symbol, Vector{𝕋}})::Nothing where 𝕋
```\n
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.

```julia-REPL
julia> using BenchmarkTools, Unitful
julia> include("./11-ORIG-julia-BaseVector.jl");
julia> par = init(Float64, 0);
julia> f = Array{par[:typ][:p], 3}(undef, (par[:cas][:NX], par[:cas][:NY], par[:lat][:int][:vel]));
julia> ρ = Array{par[:typ][:p], 2}(undef, (par[:cas][:NX], par[:cas][:NY]));
julia> 𝑢 = Array{par[:typ][:p], 2}(undef, (par[:cas][:NX], par[:cas][:NY]));
julia> 𝑣 = Array{par[:typ][:p], 2}(undef, (par[:cas][:NX], par[:cas][:NY]));
julia> 𝑏 = @benchmarkable init_equilibrium(f, ρ, 𝑢, 𝑣, latvec = par[:lat][:vec]);
julia> tune!(𝑏);
julia> run(𝑏)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  19.106 μs …  74.016 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     19.448 μs (↓)           ┊ GC (median):    0.00%
 Time  (mean ± σ):   19.489 μs ± 897.990 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
                                ↓
                     ▁▁▃▄▅▅▅▇▆▇▇▆▇▇█▆▄▃▅▂▁                      
  ▁▁▁▁▁▁▂▂▂▂▃▃▄▄▆▅▇▆████████████████████████▆▆▆▅▅▄▄▃▃▂▂▂▂▂▂▁▁▁ ▄
  19.1 μs         Histogram: frequency by time         19.8 μs <

 Memory estimate: 32 bytes, allocs estimate: 2.
```
"""
function init_equilibrium(𝑓::Array{𝕋, 3}, ρ::Array{𝕋, 2}, 𝑢::Array{𝕋, 2}, 𝑣::Array{𝕋, 2};
                          latvec::Dict{Symbol, Vector{𝕋}})::Nothing where 𝕋
    ξx  = latvec[:ξx]
    ξy  = latvec[:ξy]
    w   = latvec[:w]
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


