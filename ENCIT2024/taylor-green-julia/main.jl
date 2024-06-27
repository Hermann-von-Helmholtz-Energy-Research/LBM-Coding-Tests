# Precision parameters
const 𝕋     = Float64

# Lattice constants
const scale             = UInt(1)
const chunk             = UInt(32)
const NX                = UInt(scale * chunk)
const NY                = NX
const ndir              = UInt(9)
const mem_size_scalar   = UInt(NX * NY * sizeof(𝕋))
const mem_size_vector   = UInt(mem_size_scalar * ndir)


