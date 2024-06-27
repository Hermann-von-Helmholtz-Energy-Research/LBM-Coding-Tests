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
const mem_size_scalar   = UInt(NX * NY * sizeof(𝕋))
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



