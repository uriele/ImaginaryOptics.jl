
function meshgrid(x_range,y_range)
    Y=[ (j) for  i in x_range,j in y_range];
    X=[ (i) for i in x_range,j in y_range];
    return X,Y
end      

"""
propagation_phase(n::T,dx_norm::T,dy_norm::T,Nx::Int,Ny::Int;theta=0) where T<:AbstractFloat
propagation_phase(n::T,k0,dx_norm::M,dy_norm::M,Nx::Int,Ny::Int;;theta=0) where {T<:AbstractFloat,M}

Computes the propagation phase for a given propagation constant n, the wavelength k0 and the length L. The phase is computed as exp(-im*n*L*cos(theta)) where theta is the angle between the propagation direction and the normal to the surface. If theta=0, the phase is exp(-im*n*L). If theta=pi/2, the phase is exp(-im*n*L*cos(pi/2))=exp(im*n*L*cos(pi/2))=exp(0)=1. If theta=pi, the phase is exp(-im*n*L*cos(pi))=exp(im*n*L*cos(pi))=exp(-im*n*L)=-1. If theta=3pi/2, the phase is exp(-im*n*L*cos(3pi/2))=exp(im*n*L*cos(3pi/2))=exp(0)=1. Thus the phase is exp(-im*n*L*cos(theta)).
"""
function propagation_phase(n::T,dx::T,dy::T,Nx::Int,Ny::Int;theta=0) where T<:AbstractFloat    
    theta=oftype.(n,theta)
    Lx=(0:dx:(Nx-1)*dx).*cos(theta)
    Ly=0:dy:(Ny-1)*dy
    Ly=reshape(Ly,1,Ny).*sin(theta)
    return exp.(-im*n.*(Lx.+Ly))
end
function propagation_phase(n::T,k0,dx::M,dy::M,Nx::Int,Ny::Int;theta=0)  where {T<:AbstractFloat,M}
    dx_norm=oftype.(n,k0.*dx)
    dy_norm=oftype.(n,k0.*dy)
    return propagation_phase(n,dx,dy,Nx,Ny;theta=theta)
end


