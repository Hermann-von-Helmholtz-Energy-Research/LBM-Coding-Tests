# Precision parameters
const 𝕋                 = Float64                       # Absent in the ref. C99 code

# Lattice constants
const scale             = UInt(1)
const chunk             = UInt(32)                      # Hardcoded in the ref. C99 code
const NX                = UInt(scale * chunk)
const NY                = NX
const ndir              = UInt(9)
const mem_size_scalar   = UInt(NX * NY * sizeof(𝕋))
const mem_size_vector   = UInt(mem_size_scalar * ndir)  # Optimized w/respect to ref. C99 code


