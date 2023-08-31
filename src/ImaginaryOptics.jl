"""
    ImaginaryOptics.jl

The code uses the negative sign convention for the propagation direction. Thus the gain and losses are defined using
Thus a value of ϵ″>0 indicates gain and a value of ϵ″<0 indicates losses. The code uses the following convention for the sign of the imaginary part of the dielectric function:
```math
ϵ=ϵ′+im*ϵ″
```
"""
module ImaginaryOptics
import Unitful: μ0,ϵ0,c0,μm,ns,μs,s,cm,V,rad
import Unitful.Length as UnitLength
import Unitful.Time as UnitTime
import Unitful: ustrip,unit, @u_str,𝐋,dimension
using SparseArrays
import Arpack: eigs # for eigenvalues
@enum SOURCE_TYPE TFSF SOFT_SOURCE
@enum MODE_TYPE TE TM
@enum PML_TYPE PML_CLASSIC PML_RUMPF 
@enum BOUNDARY_CONDITIONS begin
    DIRICHLET
    NEUMANN
    ROBIN
    PERIODIC
    TRANSPARENT
end

@enum SOLVER_TYPE BEAM_PROPAGATION_METHODS DOMAIN_DECOMPOSITION_METHODS FDFD_METHODS
abstract type AbstractField{T<:AbstractFloat} end
abstract type AbstractSource{T<:AbstractFloat} end
abstract type AbstractSolver{T<:AbstractFloat} end
abstract type AbstractMaterial{T<:AbstractFloat} end
include("utility.jl")
include("common_solver_functions.jl")


export μm,μ0,ϵ0,c0
export propagation_phase
# Write your package code here.

end