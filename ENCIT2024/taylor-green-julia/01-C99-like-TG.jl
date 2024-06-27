#!/usr/bin/env julia

#----------------------------------------------------------------------------------------------#
#                                        Case Constants                                        #
#----------------------------------------------------------------------------------------------#

# Precision parameters
const 𝕋                 = Float64                       # Absent in the ref. C99 code

# Lattice constants
const scale             = UInt(1) << 0                  # 1 << n = 2^n
const chunk             = UInt(32)                      # Hardcoded in the ref. C99 code
const NX                = UInt(scale * chunk)
const NY                = NX
const ndir              = UInt(9)
const amountof_scalar   = UInt(NX * NY)
const amountof_vector   = UInt(amountof_scalar * ndir)
const mem_size_scalar   = UInt(amountof_scalar * sizeof(𝕋))
const mem_size_vector   = UInt(mem_size_scalar * ndir)  # Optimized w/respect to ref. C99 code
const w0                = 𝕋(4.0 /  9.0)     # zero velocity weight
const ws                = 𝕋(1.0 /  9.0)     # size velocity weight
const wd                = 𝕋(1.0 / 36.0)     # diag velocity weight
const wi                = (w0, ws, ws, ws, ws, wd, wd, wd, wd)      # Tuples are immutable
const dirx              = (+0, +1, +0, -1, +0, +1, -1, -1, +1)
const diry              = (+0, +0, +1, +0, -1, +1, +1, -1, -1)

# Kinematic viscosity and parameter τ
const nu                = 𝕋(1.0 / 6.0)
const tau               = 𝕋(3.0 * nu + 0.5)

# Maximum macroscopic speed
const u_max             = 𝕋(0.04 / scale)

# Fluid density
const rho0              = 𝕋(1.0)

# Simulation time steps
const NSTEPS            = UInt(round(204800 / scale / scale))


#----------------------------------------------------------------------------------------------#
#                                     Auxiliary Functions                                      #
#----------------------------------------------------------------------------------------------#

"""
`scalar_index(x::UInt, y::UInt)::UInt`\n
Returns the linear index that corresponds to the 2D position [x, y] for SCALARS.
"""
scalar_index(x::UInt, y::UInt)::UInt = NX * (y - 1) + x

"""
`field_index(x::UInt, y::UInt, d::UInt)::UInt`\n
Returns the linear index that corresponds to the 2D position [x, y] for lattice FIELDS.
"""
field_index(x::UInt, y::UInt, d::UInt)::UInt = ndir * (NX * (y - 1) + x - 1) + d

"""
`taylor_green`\n
Function to compute the exact solution for Taylor-Green vortex decay
"""
function taylor_green(t::𝕋, x::UInt, y::UInt)::NTuple{3, 𝕋}
    kx = 𝕋(2.0 * π) / NX
    ky = 𝕋(2.0 * π) / NY
    td = 𝕋(1.0) / (nu * (kx*kx + ky*ky))
    X  = 𝕋(x + 0.5)
    Y  = 𝕋(y + 0.5)
    ux = - u_max * √(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-𝕋(t) / td)
    uy = + u_max * √(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-𝕋(t) / td)
    P  = - 𝕋(0.25) * rho0 * u_max * u_max * ( (ky / kx) * cos(2kx * X)
                                             +(kx / ky) * sin(2ky * Y) )
    rh = rho0 + 𝕋(3.0) * P
    return (rh, ux, uy)
end

function taylor_green(t::𝕋, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing
    for j in UInt(1):NY
        for i in UInt(1):NX
            𝑖 = scalar_index(i, j)
            ρ[𝑖], 𝑢[𝑖], 𝑣[𝑖] = taylor_green(t, i, j)
        end
    end
end

"""
`init_equilibrium(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing`\n
Function to initialise an equilibrium particle population `f` with provided `ρ, 𝑢, 𝑣`
macroscopic fields.
"""
function init_equilibrium(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            i = scalar_index(𝑥, 𝑦)
            ϱ, 𝚞, 𝚟 = ρ[i], 𝑢[i], 𝑣[i]
            𝘂𝘂 = 𝚞 * 𝚞 + 𝚟 * 𝚟      # Optimization absent in the ref. C99 code
            for 𝑖 in UInt(1):ndir
                ξ𝘂 = 𝕋(dirx[𝑖] * 𝚞 + diry[𝑖] * 𝚟)
                𝑓[field_index(𝑥, 𝑦, 𝑖)] = wi[𝑖] * ϱ * (
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
`stream(𝑓::Vector{𝕋}, 𝑔::Vector{𝕋})::Nothing`\n
Function that performs streaming of the populations in a fully periodic domain, reading from 𝑓
and storing to 𝑔.
"""
function stream(𝑓::Vector{𝕋}, 𝑔::Vector{𝕋})::Nothing
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            for 𝑖 in UInt(1):ndir
                # "from" indices, enforcing periodicity
                𝑝 = (NX + 𝑥 - dirx[𝑖]) % NX + UInt(1)   # NX is added as to guarantee positivity
                𝑞 = (NY + 𝑦 - diry[𝑖]) % NY + UInt(1)   # NY is added as to guarantee positivity
                # Streaming from 𝑓 into 𝑔
                𝑔[field_index(𝑥, 𝑦, 𝑖)] = 𝑓[field_index(𝑝, 𝑞, 𝑖)]
            end
        end
    end
end

"""
`compute_rho_u(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing`\n
Function that computes macroscopics from mesoscopics.
"""
function compute_rho_u(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            # Initialize
            ϱ = zero(𝕋)
            𝚞 = zero(𝕋)
            𝚟 = zero(𝕋)
            𝑗 = scalar_index(𝑥, 𝑦)
            # Integrate
            for 𝑖 in UInt(1):ndir
                ϱ += 𝚏 = 𝑓[field_index(𝑥, 𝑦, 𝑖)]
                𝚞 += dirx[𝑖] * 𝚏
                𝚟 += diry[𝑖] * 𝚏
            end
            # Update
            ρ[𝑗], 𝑢[𝑗], 𝑣[𝑗] = ϱ, 𝚞, 𝚟
        end
    end
end

"""
`collide(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing`\n
Function that performs the collision operation on the particle populations using pre-computed
density and velocity values.
"""
function collide(𝑓::Vector{𝕋}, ρ::Vector{𝕋}, 𝑢::Vector{𝕋}, 𝑣::Vector{𝕋})::Nothing
    iτ = inv(tau)       # inverse:        1/τ
    cτ = one(𝕋) - iτ    # complement: 1 - 1/τ
    for 𝑦 in UInt(1):NY
        for 𝑥 in UInt(1):NX
            # Initialize
            𝑗 = scalar_index(𝑥, 𝑦)
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
                𝑓[field_index(𝑥, 𝑦, 𝑖)] = cτ * 𝑓[field_index(𝑥, 𝑦, 𝑖)] + iτ * 𝑓eq
            end
        end
    end
end


#----------------------------------------------------------------------------------------------#
#                                             Main                                             #
#----------------------------------------------------------------------------------------------#

using Format

function main(argc::Integer = length(ARGS), argv::Vector{String} = ARGS)::Integer
    # Allocate memory, without initialization
    𝑓 = Vector{𝕋}(undef, amountof_vector)
    𝑔 = Vector{𝕋}(undef, amountof_vector)
    ρ = Vector{𝕋}(undef, amountof_scalar)
    𝑢 = Vector{𝕋}(undef, amountof_scalar)
    𝑣 = Vector{𝕋}(undef, amountof_scalar)
    # Initialize ρ, 𝑢, 𝑣 with macroscopic flow
    taylor_green(zero(𝕋), ρ, 𝑢, 𝑣)
    # Initialize 𝑓 at equilibrium
    init_equilibrium(𝑓, ρ, 𝑢, 𝑣)
    # Main loop
    for n in 1:NSTEPS
        # Stream
        stream(𝑓, 𝑔)
        # Macros
        compute_rho_u(𝑔, ρ, 𝑢, 𝑣)
        # Collide
        collide(𝑔, ρ, 𝑢, 𝑣)
        # (𝑓, 𝑔) swapping
        𝑓, 𝑔 = 𝑔, 𝑓
        # DEBUG
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

exit(main())

