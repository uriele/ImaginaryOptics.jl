


""" 
    BPM_solver{T}(Nx::Int,Ny::Int,dx_norm::T,dy_norm::T,ERzz::Matrix{Complex{T}};nonlinear::Bool=true,Npmly=10,kwargs...) where T<:AbstractFloat
    BPM_solver{T}(NS::NTuple{2,Int},RES::NTuple{2,T},ERzz;kwargs...) where T<:AbstractFloat
    BPM_solver{T}(NS::NTuple{2,Int},dx_norm::T,ERzz;kwargs...)       where T<:AbstractFloat
    BPM_solver{T}(Nx::Int,Ny::Int,dx_norm::T,ERzz;kwargs...)         where T<:AbstractFloat

    BPM_solver is a subtype of AbstractDomain that contains the information to solve the numerical problem. As the name says it is used to solve the problem using the beam
    propagation method (BPM). The current version uses the nonlinear Crank Nicolson method to the direct problem and autodifferentiation to solve the inverse problem. The current
    version only supports the paraxial approximation, future version will support Pade approximants to solve large angles.

    # INPUTS
    - `Nx::Int`: number of points in the x direction
    - `Ny::Int`: number of points in the y direction
    - `dx_norm::T`: grid spacing in the x direction in normalized units (i.e. dx_norm=dx*k₀)
    - `dy_norm::T`: grid spacing in the y direction in normalized units (i.e. dy_norm=dy*k₀)
    - `ERzz ::Matrix{Complex{T}}`: matrix containing the real part of the dielectric function. It is a misnomer, since in the case of active devices it 
    can be complex and equal to ϵ_core-k_min^2.
    - `nonlinear::Bool=false`: if true, the nonlinear term is included in the BPM equation (BPM has been made to solve the nonlinear optimization in 1D and 2D, by default it is true)
    - `wpmly=10::Int=10`: number of PML points in the direction perpendicular to the propagation direction (y direction)
    See also: [DomainFDFD](@ref), [DomainDecomposition](@ref), [DomainDecompositionLowMemory](@ref),
"""
struct BPM_solver{T} <: AbstractDomain{T}
    nonlinear     :: Bool
    N             :: Tuple{Int,Int}
    laplacian     :: SparseMatrixCSC{Complex{T}}
    epsilon       :: Matrix{T}
    
    function BPM_solver(Nx::Int,Ny::Int,dx_norm::T,dy_norm::T,epsilon::Vector{T};
        nonlinear::Bool=true,wpmly=10,kwargs...) where T<:AbstractFloat   
        
        
        @info("Nx: $Nx, Ny: $Ny, dx_norm: $dx_norm, dy_norm: $dy_norm, Npmly: $Npmly")
        _,∇=Laplacian(1,Ny,dx_norm,dy_norm;wpmlx=0,wpmly=wpmly,check=2,kwargs...)
        @debug("Nx: $Nx, Ny: $Ny, dx_norm: $dx_norm, dy_norm: $dy_norm, Npmly: $Npmly")
        _,∇=Laplacian(1,Ny,dx_norm,dy_norm;wpmlx=0,wpmly=wpmly,check=2,kwargs...)
        # Assume μ=1 that is the case for most materials
        laplacian=∇[4]*∇[2]
       
        @debug("lapacialian created")
        N=(Nx,Ny)
        epsilon=reshape(epsilon,Nx,Ny)
        new{T}(nonlinear,(Nx,Ny),laplacian,epsilon)
    end
end

BPM_solver(NS::NTuple{2,Int},RES::NTuple{2,T},epsilon_real::Union{VecOrMat{T},VecOrMat{Complex{T}}};kwargs...) where T<:AbstractFloat =BPM_solver(NS...,RES...,epsilon_real;kwargs...)
BPM_solver(NS::NTuple{2,Int},dx_norm::T,epsilon_real::Union{VecOrMat{T},VecOrMat{Complex{T}}};kwargs...)       where T<:AbstractFloat =BPM_solver(NS...,dx_norm,dx_norm,epsilon_real;kwargs...)
BPM_solver(Nx::Int,Ny::Int,dx_norm::T,epsilon_real::Union{VecOrMat{T},VecOrMat{Complex{T}}};kwargs...)         where T<:AbstractFloat =BPM_solver(Nx,Ny,dx_norm,dx_norm,epsilon_real;kwargs...)

