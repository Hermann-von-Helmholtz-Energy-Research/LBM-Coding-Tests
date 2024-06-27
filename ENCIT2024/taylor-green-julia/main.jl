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
const NSTEPS            = UInt(204800 / scale / scale)


#----------------------------------------------------------------------------------------------#
#                                     Auxiliary Functions                                      #
#----------------------------------------------------------------------------------------------#

"""
`scalar_index(x::UInt, y::UInt)`\n
Returns the linear index that corresponds to the 2D position [x, y] for SCALARS.
"""
scalar_index(x::UInt, y::UInt) = NX * y + x

"""
`field_index(x::UInt, y::UInt, d::UInt = ndir)`\n
Returns the linear index that corresponds to the 2D position [x, y] for lattice FIELDS.
"""
field_index(x::UInt, y::UInt, d::UInt = ndir) = NX * (NY * d + y) + x

function taylor_green end
function init_equilibrium end
function stream end
function compute_rho_u end
function collide end


#----------------------------------------------------------------------------------------------#
#                                             Main                                             #
#----------------------------------------------------------------------------------------------#

using StaticArrays

function main()::Integer
    # Allocate (initialized) memory
    𝑓 = @SVector zeros(amountof_vector)
    𝑔 = @SVector zeros(amountof_vector)
    ρ = @SVector zeros(amountof_scalar)
    𝑢 = @SVector zeros(amountof_scalar)
    𝑣 = @SVector zeros(amountof_scalar)
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
    end
    # Return
    return 0
end

exit(main())

