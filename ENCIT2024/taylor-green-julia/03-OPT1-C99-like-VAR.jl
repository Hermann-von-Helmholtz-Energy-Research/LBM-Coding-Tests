#!/usr/bin/env julia

#==============================================================================================#
#                             Improved, still C99-like julia code                              #
#==============================================================================================#

#----------------------------------------------------------------------------------------------#
#                                     Auxiliary Functions                                      #
#----------------------------------------------------------------------------------------------#

"""
`scalar_index(x::UInt, y::UInt; NX::UInt)::UInt`\n
Returns the linear index that corresponds to the 2D position [x, y] for SCALARS.
"""
scalar_index(x::UInt, y::UInt; NX::UInt)::UInt = NX * (y - 1) + x

"""
`field_index(x::UInt, y::UInt, d::UInt, NX::UInt, ndir::UInt)::UInt`\n
Returns the linear index that corresponds to the 2D position [x, y] for lattice FIELDS.
"""
field_index(x::UInt, y::UInt, d::UInt; NX::UInt, ndir::UInt)::UInt = begin
    ndir * (NX * (y - 1) + x - 1) + d
end

"""
`taylor_green`\n
Function to compute the exact solution for Taylor-Green vortex decay
"""
function taylor_green(t::𝕋, x::UInt, y::UInt;
    NX::UInt, NY::UInt, nu::𝕋, rho0::𝕋, u_max::𝕋
)::NTuple{3, 𝕋} where 𝕋
    kx = 𝕋(2.0 * π) / NX
    ky = 𝕋(2.0 * π) / NY
    td = 𝕋(1.0) / (nu * (kx*kx + ky*ky))
    X  = 𝕋(x - NX / 𝕋(2.0))     # Centered vortex
    Y  = 𝕋(y - NY / 𝕋(2.0))     # Centered vortex
    ux = - u_max * √(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-t / td)
    uy = + u_max * √(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-t / td)
    P  = - 𝕋(0.25) * rho0 * u_max * u_max * ( (ky / kx) * cos(2kx * X)
                                             +(kx / ky) * sin(2ky * Y) )
    rh = rho0 + 𝕋(3.0) * P
    return (rh, ux, uy)
end

function taylor_green(t::𝕋, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋};
    NX::UInt, NY::UInt, nu::𝕋, rho0::𝕋, u_max::𝕋
)::Nothing where 𝕋
    for j in UInt(1):NY
        for i in UInt(1):NX
            𝑖 = scalar_index(i, j, NX=NX)
            ρ[𝑖], 𝑢[𝑖], 𝑣[𝑖] = taylor_green(t, i, j, NX=NX, NY=NY, nu=nu, rho0=rho0, u_max=u_max)
        end
    end
end

"""
`init_equilibrium`\n
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.
"""
function init_equilibrium(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋};
    NX::UInt, NY::UInt, ndir::UInt, dirx::Vector{Int64}, diry::Vector{Int64}, wi::Vector{𝕋}
)::Nothing where 𝕋
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            i = scalar_index(𝑥, 𝑦, NX=NX)
            ϱ, 𝚞, 𝚟 = ρ[i], 𝑢[i], 𝑣[i]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟      # Optimization absent in the ref. C99 code
            for 𝑖 in UInt(1):ndir
                ξ𝘂 = 𝕋(dirx[𝑖] * 𝚞 + diry[𝑖] * 𝚟)
                𝑓[field_index(𝑥, 𝑦, 𝑖, NX=NX, ndir=ndir)] = wi[𝑖] * ϱ * (
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
`stream`\n
Function that performs streaming of the populations in a fully periodic domain, reading from 𝑓
and storing to 𝑔.
"""
function stream(𝑓::Vector{𝕋}, 𝑔::Vector{𝕋};
    NX::UInt, NY::UInt, ndir::UInt, dirx::Vector{Int64}, diry::Vector{Int64}
)::Nothing where 𝕋
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            for 𝑖 in UInt(1):ndir
                # "from" indices, enforcing periodicity
                𝑝 = (NX + 𝑥 - dirx[𝑖]) % NX + UInt(1)   # NX is added as to guarantee positivity
                𝑞 = (NY + 𝑦 - diry[𝑖]) % NY + UInt(1)   # NY is added as to guarantee positivity
                # Streaming from 𝑓 into 𝑔
                𝑔[field_index(𝑥, 𝑦, 𝑖, NX=NX, ndir=ndir)] =
                    𝑓[field_index(𝑝, 𝑞, 𝑖, NX=NX, ndir=ndir)]
            end
        end
    end
end

"""
`compute_rho_u`\n
Function that computes macroscopics from mesoscopics.
"""
function compute_rho_u(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋};
    NX::UInt, NY::UInt, ndir::UInt, dirx::Vector{Int64}, diry::Vector{Int64}
)::Nothing where 𝕋
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            # Initialize
            ϱ = zero(𝕋)
            𝚞 = zero(𝕋)
            𝚟 = zero(𝕋)
            𝑗 = scalar_index(𝑥, 𝑦, NX=NX)
            # Integrate
            for 𝑖 in UInt(1):ndir
                ϱ += 𝚏 = 𝑓[field_index(𝑥, 𝑦, 𝑖, NX=NX, ndir=ndir)]
                𝚞 += dirx[𝑖] * 𝚏
                𝚟 += diry[𝑖] * 𝚏
            end
            # Update
            ρ[𝑗], 𝑢[𝑗], 𝑣[𝑗] = ϱ, 𝚞, 𝚟
        end
    end
end

"""
`collide`\n
Function that performs the collision operation on the particle populations using pre-computed
density and velocity values.
"""
function collide(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋};
    NX::UInt, NY::UInt, ndir::UInt, dirx::Vector{Int64}, diry::Vector{Int64},
    wi::Vector{𝕋}, tau::𝕋
)::Nothing where 𝕋
    iτ = inv(tau)       # inverse:        1/τ
    cτ = one(𝕋) - iτ    # complement: 1 - 1/τ
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            # Initialize
            𝑗 = scalar_index(𝑥, 𝑦, NX=NX)
            ϱ, 𝚞, 𝚟 = ρ[𝑗], 𝑢[𝑗], 𝑣[𝑗]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟      # Optimization absent in the ref. C99 code
            for 𝑖 in UInt(1):ndir
                ξ𝘂 = 𝕋(dirx[𝑖] * 𝚞 + diry[𝑖] * 𝚟)
                # Equilibrium
                𝑓eq = wi[𝑖] * ϱ * (
                    + 𝕋(1.0)
                    + 𝕋(3.0) * ξ𝘂
                    + 𝕋(4.5) * ξ𝘂 * ξ𝘂
                    - 𝕋(1.5) * 𝘂𝘂
                )
                # Relax to equilibrium
                𝑓[field_index(𝑥, 𝑦, 𝑖, NX=NX, ndir=ndir)] =
                    cτ * 𝑓[field_index(𝑥, 𝑦, 𝑖, NX=NX, ndir=ndir)] + iτ * 𝑓eq
            end
        end
    end
end


#----------------------------------------------------------------------------------------------#
#                                             Main                                             #
#----------------------------------------------------------------------------------------------#

function main(𝕋::Type{T} where T<:AbstractFloat, log2scale::Unsigned)::Integer
    #--------------------------------------------------------------------------------------#
    #                            Constants (C99 code's globals)                            #
    #--------------------------------------------------------------------------------------#
    scale             = UInt(1) << log2scale      # 1 << n = 2^n
    chunk             = UInt(32)                  # Hardcoded in the ref. C99 code
    NX                = UInt(scale * chunk)
    NY                = NX
    ndir              = UInt(9)
    amountof_scalar   = UInt(NX * NY)
    amountof_vector   = UInt(NX * NY * ndir)
    w0                = 𝕋(4.0 /  9.0)             # zero velocity weight
    w1                = 𝕋(1.0 /  9.0)             # size velocity weight
    w2                = 𝕋(1.0 / 36.0)             # diag velocity weight
    wi                = [w0, w1, w1, w1, w1, w2, w2, w2, w2]
    dirx              = [+0, +1, +0, -1, +0, +1, -1, -1, +1]
    diry              = [+0, +0, +1, +0, -1, +1, +1, -1, -1]
    nu                = 𝕋(1.0 / 6.0)
    tau               = 𝕋(3.0 * nu + 0.5)
    u_max             = 𝕋(0.04 / scale)
    rho0              = 𝕋(1.0)
    NSTEPS            = UInt(round(204800 / scale / scale))
    #--------------------------------------------------------------------------------------#
    #                       Allocate memory, without initialization                        #
    #--------------------------------------------------------------------------------------#
    𝑓 = Vector{𝕋}(undef, amountof_vector)
    𝑔 = Vector{𝕋}(undef, amountof_vector)
    ρ = Vector{𝕋}(undef, amountof_scalar)
    𝑢 = Vector{𝕋}(undef, amountof_scalar)
    𝑣 = Vector{𝕋}(undef, amountof_scalar)
    # Initialize ρ, 𝑢, 𝑣 with macroscopic flow
    taylor_green(zero(𝕋), ρ, 𝑢, 𝑣, NX=NX, NY=NY, nu=nu, rho0=rho0, u_max=u_max)
    # Initialize 𝑓 at equilibrium
    init_equilibrium(𝑓, ρ, 𝑢, 𝑣, NX=NX, NY=NY, ndir=ndir, dirx=dirx, diry=diry, wi=wi)
    # Main loop
    for n in 1:NSTEPS
        # Stream
        stream(𝑓, 𝑔, NX=NX, NY=NY, ndir=ndir, dirx=dirx, diry=diry)
        # Macros
        compute_rho_u(𝑔, ρ, 𝑢, 𝑣, NX=NX, NY=NY, ndir=ndir, dirx=dirx, diry=diry)
        # Collide
        collide(𝑔, ρ, 𝑢, 𝑣, NX=NX, NY=NY, ndir=ndir, dirx=dirx, diry=diry, wi=wi, tau=tau)
        # (𝑓, 𝑔) swapping
        𝑓, 𝑔 = 𝑔, 𝑓
        # PROGRESS
        # if (n % 128 == 0) || (n == NSTEPS)
        #     if (n % 8192 == 0) || (n == NSTEPS)
        #         println(format(" ({1:6d}: {2:5.1f}%)", n, 𝕋(100n)/𝕋(NSTEPS)))
        #     else
        #         print(".")
        #     end
        # end 
    end
    #--------------------------------------------------------------------------#
    #    Memory de-allocation is automatically performed by julia's garbage    #
    #        collector when the 𝑓, 𝑔, ρ, 𝑢, 𝑣 Vectors are out of scope.        #
    #--------------------------------------------------------------------------#
    # Return
    return 0
end


