function yeeder2d_grid(ER2,UR2)
    #Rumpf pp.132
    # To remember, look at the matrix sx_ez
    # E has 1,1 -> _ey has 1,2 -> _ex has 2,1 : the are the opposite of the z one in their comonent
    # H has 2,2 -> _ey has 2,1 -> _ex has 1,2 : the are the opposite of the z one in their comonent
   
    Nx2,Ny2=size(ER2);
    ERxx =ER2[2:2:Nx2,1:2:Ny2];
    ERyy =ER2[1:2:Nx2,2:2:Ny2];
    ERzz =ER2[1:2:Nx2,1:2:Ny2];
    URxx =UR2[1:2:Nx2,2:2:Ny2];
    URyy =UR2[2:2:Nx2,1:2:Ny2];
    URzz =UR2[2:2:Nx2,1:2:Ny2];
    return ERxx,ERyy,ERzz,URxx,URyy,URzz
end

"""
yeeder2d(NS,RES_unit;periodic=false,kinc,k0);
yeeder2d(NS,RES_norm;periodic=false,kinc);

Creates the derivative matrices for the finite difference frequency domain method for the staggered Yee grid.
    
 # INPUT
 
    - `NS::Tuple{Int,Int}`: number of points in the x and y direction
    - `RES::Tuple{T,T}`: grid spacing in the x and y direction in normalized units (i.e. dx_norm=dx*k₀)

=================
NS    (Nx, Ny) Grid Size
RES   (dx, dy) Grid Resolution
BC    [xbc ybc] Boundary Conditions
        0: Dirichlet boundary conditions
        1: Periodic boundary conditions
kinc  [kx ky] Incident Wave Vector
      This argument is only needed for PBCs.
%
Note: For normalized grids use k0*RES and kinc/k0
%
OUTPUT ARGUMENTS
=================
DEX   Derivative Matrix wrt x for Electric Fields
DEY   Derivative Matrix wrt to y for Electric Fields
DHX   Derivative Matrix wrt to x for Magnetic Fields
DHY   Derivative Matrix wrt to y for Magnetic Fields
"""
function yeeder2d(NS::Tuple{Int,Int},RES_norm::Tuple{T,T};boundary_conditions::BOUNDARY_CONDITIONS=DIRICHLET) where T<:AbstractFloat
    @assert(NS[1]>=1 && NS[2]>=1, "NS must be greater or equal to 1")
    @debug("use of yeeder2d unitless")
    RES2_norm=(RES_norm[1]/2,RES_norm[2]/2)
    # Build DEX
    M=prod(NS);
    # center diagonal
    d0=-ones(M);
    d1= ones(M-1);
    d1[NS[1]:NS[1]:M-1].=0.0;

    
    (NS[1]>1) ? DEX=spdiagm(0=>d0,1=>d1) : DEX=spzeros(M,M);
    d1= ones(M-NS[1]);
    
    (NS[2]>1) ? DEY=spdiagm(0=>d0,NS[1]=>d1) : DEY=spzeros(M,M);

    return (DEX./RES_norm[1],DEY./RES_norm[2], -1 .*transpose(copy(DEX))./RES_norm[1], -1 .*transpose(copy(DEY))./RES_norm[2]) 
end

function yeeder2d(NS::Tuple{Int,Int},RES_unit::Tuple{UnitLength,UnitLength};λ₀=1.55μm,theta=0)
    @debug("use of yeeder2d with units")
    λ₀=oftype(RES_unit[1],λ₀)
    k₀=2*π/λ₀;          # 1/𝐋
    
    RES_norm=(k₀*RES_unit[1],k₀*RES_unit[2])

    return yeeder2d(NS,RES_norm)

    # Build DEX
    M=prod(NS);
    # center diagonal
    d0=-ones(M);
    d1= ones(M-1);
    d1[NS[1]:NS[1]:M-1].=0.0;

    (NS[1]>1) ? DEX=spdiagm(0=>d0,1=>d1) : DEX=spzeros(M,M);

    # Build DEY
    d1= ones(M-NS[1]);
    
    (NS[2]>1) ? DEY=spdiagm(0=>d0,NS[1]=>d1) : DEY=spzeros(M,M);

    return (DEX./RES[1],DEY./RES[2], -1 .*transpose(copy(DEX))./RES[1], -1 .*transpose(copy(DEY))./RES[2]) 
end
yeeder2d(Nx,Ny,dx,dy;λ₀=1.55μm)=yeeder2d((Nx,Ny),(dx,dy);λ₀=1.55μm)
yeeder2d(Nx,Ny,dx,dy)=yeeder2d((Nx,Ny),(dx,dy))

"""
    Laplacian(Nx::Int,Ny::Int,dx::T,dy::T;μ_xx=nothing,μ_yy=nothing,ϵ_xx=nothing,ϵ_yy=nothing,wpmlx=0,wpmly=0,mode_type::MODE_TYPE=TE,use_pml=true,check=true, kwargs...) where T<:Union{AbstractFloat,Unitful.Length}

    Laplacian is a function that creates the laplacain operators for the finite difference frequency domain method. It is a wrapper for the yeeder2d function that creates the curl operators and the apply_scpml_TE function that applies the perfectly matched layer to the curl operators. It also creates the laplacian operator for the finite difference frequency domain method.
    The implementation is based on Rumpft work.
    
    # INPUTS
    - `Nx::Int`: number of points in the x direction
    - `Ny::Int`: number of points in the y direction
    - `dx::T`: grid spacing in the x direction in normalized units (i.e. dx_norm=dx*k₀)
    - `dy::T`: grid spacing in the y direction in normalized units (i.e. dy_norm=dy*k₀)
    - `μ_xx::Vector{T}`: vector containing the real part of the magnetic permeability in the x direction. It is a misnomer, since in the case of active devices it
    can be complex and equal to μ_core-k_min^2.
    - `μ_yy::Vector{T}`: vector containing the real part of the magnetic permeability in the y direction.
    - `ϵ_xx::Vector{T}`: vector containing the real part of the dielectric function in the x direction. 
    - `ϵ_yy::Vector{T}`: vector containing the real part of the dielectric function in the y direction. It is a misnomer, since in the case of active devices it
    can be complex and equal to ϵ_core-k_min^2.
    - `wpmlx::Int=0`: number of PML points in the direction perpendicular to the propagation direction (x direction)
    - `wpmly::Int=0`: number of PML points in the direction perpendicular to the propagation direction (y direction)
    - `mode_type::MODE_TYPE=TE`: type of mode to be solved, TE or TM
    - `use_pml::Bool=true`: if true, the perfectly matched layer is applied to the curl operators
    - `check::Bool=true`: if true, the function checks if the number of PML points is zero, if it is zero it will throw an error, if false it will ignore the error

    # OUTPUTS
    - `∇²::SparseMatrixCSC{Complex{T}}`: laplacian operator
    - `∇ᵉ_x::SparseMatrixCSC{Complex{T}}`: gradient operator in the x direction for the electric field
    - `∇ᵉ_y::SparseMatrixCSC{Complex{T}}`: gradient operator in the y direction for the electric field
    - `∇ʰ_x::SparseMatrixCSC{Complex{T}}`: gradient operator in the x direction for the magnetic field
    - `∇ʰ_y::SparseMatrixCSC{Complex{T}}`: gradient operator in the y direction for the magnetic field


"""
function Laplacian(Nx,Ny,dx::T,dy::T;μ_xx=nothing,μ_yy=nothing,ϵ_xx=nothing,ϵ_yy=nothing,wpmlx=0,wpmly=0,mode_type::MODE_TYPE=TE,use_pml=true,check=true, kwargs...) where T<:Union{AbstractFloat}
    if (wpmlx==0 || wpmly==0) 
        
        if (check==true) 
            @error "Npmlx or Npmly is zero, this will result in a non-physical PML
            set check=false to ignore this error" 
            return nothing
        else
            @warn "Npmlx or Npmly is zero, this will result in a non-physical PML
            Npmlx: $(wpmlx) Npmly: $(wpmly)"
        end
    end
    NS=(Nx,Ny)
    RES=(dx,dy) 
    @debug ("NS: $(NS), RES: $(RES)")
    @debug("NPML: $((wpmlx,wpmlx,wpmly,wpmly))")
    ∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y=yeeder2d(Nx,Ny,dx,dy)
    if use_pml
        (∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y)=apply_scpml_TE(∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y,Nx,Ny,(dx,dy);NPML=(wpmlx,wpmlx,wpmly,wpmly),kwargs...);
    end


    if (isnothing(μ_xx) && isnothing(μ_yy) && isnothing(ϵ_xx) && isnothing(ϵ_yy))
        @warn "No material properties specified, the laplacian will not be created only the curl operators"
        ∇²=nothing
        return (∇²,∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y)
    end
       
    if (mode_type==TE || mode_type==H) 
        @assert μ_xx!=nothing && μ_yy!=nothing "μ_xx and μ_yy must be specified for TE mode"
        ONE=one(typeof(μ_yy[1]))
        @debug "∇²: $(size(∇ᵉ_x)) $(size(∇ʰ_x)) $(size(∇ᵉ_y)) $(size(∇ʰ_y)) μ_xx: $(size(μ_xx)) μ_yy: $(size(μ_yy))"
        ∇²=∇ʰ_x*spdiagm(0=>ONE./μ_yy[:])*∇ᵉ_x.+∇ʰ_y*spdiagm(0=>ONE./μ_xx[:])*∇ᵉ_y
        return (∇²,(∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y))
    elseif (mode_type==TM || mode_type==E)
        @assert ϵ_xx!=nothing && ϵ_yy!=nothing "ϵ_xx and ϵ_yy must be specified for TM mode"
        ONE=one(typeof(ϵ_xx[1]))
        ∇²=∇ᵉ_x*spdiagm(0=>ONE./ϵ_yy[:])*∇ʰ_x.+DEY*spdiagm(0=>ONE./ϵ_xx[:])*∇ʰ_y
        return ∇²,∇ᵉ_x,∇ᵉ_y,∇ʰ_x,∇ʰ_y
    end

    error("Mode has to be TE(H), TM(E)")
end


"""
    SCPML_classic(a_max,σ_max,p,i,Nx)
Classic implementation of stretched coordinates perfectly matched layer (SCPML) for the finite difference frequency domain method (FDFD).
```math

σ=σ_max(i/Nx)^p

s=a_max(1+im*σ)
```
# Arguments
- `a_max::Float`: maximum value of the stretching factor
- `σ_max::Float`: maximum value of the conductivity
- `p::Float`: stretching factor
- `i::Int`: current position in the grid
"""
function SCPML_classic(a_max,σ_max,p,i,Nx)
    σ=σ_max*(i/Nx)^p;
    return a_max*(1-im*σ);
end

"""
    SCPML_rumpf(a_max,σ_max,p,i,Nx)
Implementation of stretched coordinates perfectly matched layer (SCPML) for the finite difference frequency domain method (FDFD) based on the work of Rumpft
```math

a=1+(a_max-1)*(i/Nx)^p
σ=σ_max*sin(π*i/(2*Nx))^2

s=a_max(1+im*σ*60*σ)
```
# Arguments
- `a_max::Float`: maximum value of the stretching factor
- `σ_max::Float`: maximum value of the conductivity
- `p::Float`: stretching factor
- `i::Int`: current position in the grid
"""
function SCPML_rumpf(a_max,σ_max,p,i,Nx)
    a=1+(a_max-1)*(i/Nx)^p;
    σ=σ_max*sin(π*i/(2*Nx))^2;

    return a*(1-im*60*σ);    
end

function SC_PML(Nx::Int,Ny::Int,RES_norm::NTuple{2};
    NPML::NTuple{4,Int}=(20,20,20,20),m=2.5,R=1e-7,λₘᵢₙ=1.55μm,
    pml_type::PML_TYPE=PML_CLASSIC,kwargs...)    
 
    sx2 = ones(Complex{Float32}, 2*Nx,2*Ny);
    sy2 = ones(Complex{Float32}, 2*Nx,2*Ny);
    (Nx2,Ny2)=2 .*(Nx,Ny);
    NPML2=[2 .*npml for npml in NPML]
    RES2_norm=(RES_norm[1]/2,RES_norm[2]/2)
    ###########################################
    d_PML=(NPML2[1]*RES2_norm[1],NPML2[2]*RES2_norm[1],NPML2[3]*RES2_norm[2],NPML2[4]*RES2_norm[2]);    
    if pml_type==PML_CLASSIC
        @debug("classic PML")
        σₘₐₓ=((-(m+1)*log(R)./(2 .*d_PML))); # unitless
        @info("σₘₐₓ: $(σₘₐₓ)")
        aₘₐₓ=1.0
        
    elseif pml_type==PML_RUMPF
        @debug("rumpf PML")
        σₘₐₓ=[1.0,1.0,1.0,1.0]
        aₘₐₓ=4.0
    else
        error("incorrect PML pml_type")
    end
    #return pml,σₘₐₓ
    @debug("σₘₐₓ: $(σₘₐₓ)")
    pml(σₘₐₓ,i,Nx)=(pml_type==PML_CLASSIC) ? SCPML_classic(aₘₐₓ,σₘₐₓ,m,i,Nx) :  SCPML_rumpf(aₘₐₓ,σₘₐₓ,m,i,Nx) 
    
    for nx=1:NPML2[1]
        sx2[NPML2[1]-nx+1,:].=pml(σₘₐₓ[1],float(nx),float(NPML2[1]));
    end
    
    for nx=1:NPML2[2]
        sx2[Nx2-NPML2[2]+nx,:].=pml(σₘₐₓ[2],float(nx),float(NPML2[2]));
    end
    
    for ny=1:NPML2[3]
        sy2[:,NPML2[3]-ny+1].=pml(σₘₐₓ[3],float(ny),float(NPML2[3]));
    end
    
    for ny=1:NPML2[4]
        sy2[:,Ny2-NPML2[4]+ny,:].=pml(σₘₐₓ[4],float(ny),float(NPML2[4]));
    end
   # To remember, look at the matrix sx_ez
   # E has 1,1 -> _ey has 1,2 -> _ex has 2,1 : the are the opposite of the z one in their comonent
   # H has 2,2 -> _ey has 2,1 -> _ex has 1,2 : the are the opposite of the z one in their comonent
   sx_ey = sx2[1:2:Nx2,2:2:Ny2]
   sx_ez = sx2[1:2:Nx2,1:2:Ny2]
   sy_ex = sy2[2:2:Nx2,1:2:Ny2]
   sy_ez = sy2[1:2:Nx2,1:2:Ny2]
   
   sx_hy = sx2[2:2:Nx2,1:2:Ny2]
   sx_hz = sx2[2:2:Nx2,2:2:Ny2]
   sy_hx = sy2[1:2:Nx2,2:2:Ny2]
   sy_hz = sy2[2:2:Nx2,2:2:Ny2]
    #returns the stretching to be used to update both the matrix Dx and Dy AND the Robin BC
    return  sx_ey,sx_ez,sy_ex,sy_ez,sx_hy,sx_hz,sy_hx,sy_hz #,sx,sy

end



function apply_scpml_TE(DEX,DEY,DHX,DHY,Nx::Int,Ny::Int,(dx,dy);kwargs...)
    @debug "kwargs: $(keys(kwargs))"
    (sx_ey,sx_ez,sy_ex,sy_ez,sx_hy,sx_hz,sy_hx,sy_hz)=SC_PML(Nx,Ny,(dx,dy);kwargs...)

    DEX=spdiagm(0=>1.0 ./(sx_hy[:]))*DEX;
    DEY=spdiagm(0=>1.0 ./(sy_hx[:]))*DEY;
    DHX=spdiagm(0=>1.0 ./(sx_ez[:]))*DHX;
    DHY=spdiagm(0=>1.0 ./(sy_ez[:]))*DHY;
    return (DEX,DEY,DHX,DHY,dx,dy)
end

function apply_scpml_TM(DEX,DEY,DHX,DHY,Nx::Int,Ny::Int,(dx,dy);kwargs...)
    
    (sx_ey,sx_ez,sy_ex,sy_ez,sx_hy,sx_hz,sy_hx,sy_hz)=SC_PML(Nx,Ny,(dx,dy);kwargs...)

    DEX=spdiagm(0=>1.0 ./(sx_hz[:]))*DEX;
    DEY=spdiagm(0=>1.0 ./(sy_hz[:]))*DEY;
    DHX=spdiagm(0=>1.0 ./(sx_ey[:]))*DHX;
    DHY=spdiagm(0=>1.0 ./(sy_ex[:]))*DHY;
    return (DEX,DEY,DHX,DHY,dx,dy)
end

#=
function apply_scpml_TE(DEX,DEY,DHX,DHY,Nx::Int,Ny::Int,(dx,dy);kwargs...)
    (_,sx_ez,_,sy_ez,sx_hy,_,sy_hx,_
    )=SC_PML(Nx,Ny,(dx,dy);kwargs...)

    DEX=spdiagm(0=>1.0 ./(sx_hy[:]))*DEX;
    DEY=spdiagm(0=>1.0 ./(sy_hx[:]))*DEY;
    DHX=spdiagm(0=>1.0 ./(sx_ez[:]))*DHX;
    DHY=spdiagm(0=>1.0 ./(sy_ez[:]))*DHY;
    return (DEX,DEY,DHX,DHY,dx,dy)
end
=#


function slabmode(W::T,Ny,dy_norm::V; mode= 0, k₀=2π/1.55μm,n_core,n_cladding,n_substrate=nothing,shift=0.0μm) where {T,V<:AbstractFloat}
    if isnothing(n_substrate);
        n_substrate=n_cladding
    end
    
    println("k₀   : $k₀");
    NA=sqrt(n_core^2-n_cladding^2);  
    @debug("k₀: $k₀ NA: $NA W: $W")
    R=k₀*W*NA;  
    M_allowed=Int(floor(2*R/pi));
    println("Numerical Aperture   : $NA")
    println("Normalized Frequency : $R")
    println("Mode Allowed : $M_allowed")
    N_used=Int(min(mode,M_allowed)+1);
    
    Ny_buffer=Int(floor(Ny/2));
    Ny=Ny+2*Ny_buffer;
    Ny2=Int(Ny*2);
    dy2_norm=dy_norm/2;
    dy2_unit=dy2_norm/k₀;
    ER2=ones(1,Ny2).*n_cladding^2;
    UR2=ones(1,Ny2);
    Y2_wg=range(1,Ny2,Ny2).*dy2_unit;
    Y2_wg= Y2_wg.-mean(Y2_wg);
    
    ER2= [abs(y.+shift).<= W/2 ? n_core.^2 : (((y.+shift).<-W/2) ? n_substrate.^2 : n_cladding.^2) for y in Y2_wg];
    ER2=reshape(ER2,1,Ny2);

    # EXTRACT YEE GRID ARRAYS
    ERxx = ER2[1,1:2:Ny2]; ERxx=spdiagm(0=> ERxx[:]);
    ERyy = ER2[1,2:2:Ny2]; ERyy=spdiagm(0=> ERyy[:]);
    ERzz = ER2[1,1:2:Ny2]; ERzz=spdiagm(0=> ERzz[:]);
    URxx = UR2[1,2:2:Ny2]; URxx=spdiagm(0=> URxx[:]);
    URyy = UR2[1,1:2:Ny2]; URyy=spdiagm(0=> URyy[:]);
    URzz = UR2[1,2:2:Ny2]; URzz=spdiagm(0=> URzz[:]);
    

    NS  = (1, Ny);
    RES = (ones(V),dy_norm);
    (DEX,DEY,DHX,DHY) = yeeder2d(NS,RES);
    A = -(DHY/URzz*DEY + ERxx);
    B = inv(Matrix(URyy));
    ev      = -n_core.^2;
    eigenvalues, eigenvectors = eigs(Matrix(A), B,which=:SR, check=2) #N_used+1,  which=:SR);
    
    D       = [sqrt(Complex{V}(x)) for x in eigenvalues];
    NEFF    = -im*diagm(D);
    
    return real(NEFF)[:],eigenvectors[Ny_buffer+1:end-Ny_buffer,:],ERzz[Ny_buffer+1:end-Ny_buffer,:]

end




mutable struct Tensor
    xx::Array;
    yy::Array;
    zz::Array;
    xy::Array;
    xz::Array;
    yz::Array;
end

function create_FDFD2D(device,RES::NTuple{2};
    n_core,n_cladding,λ=1.55μm,source_position,NPMLx=20,NPMLy=20,padding_x=0.0μm,padding_y=0.0μm)
    #Nx=100;
    #Ny=100;
    #Δλ=0.01;
    #λ0=1.55;
    #device=Device(λ
    dx=RES[1];
    dy=RES[2];

    NPML=(NPMLx,NPMLx,NPMLy,NPMLy);

    xlim = (minimum(x for (_, (x, _)) in device.ports),maximum(x for (_, (x, _)) in  device.ports));
    ylim = (minimum(y for ((_, y)) in device.polygon),maximum(y for ((_, y)) in  device.polygon));
    
    xlim=(xlim[1]-padding_x,xlim[2].+padding_x)
    ylim=(ylim[1]-padding_y,ylim[2].+padding_y)
    Nx=Int(Unitful.ustrip(ceil((xlim[2]-xlim[1])/dx)));
    Ny=Int(Unitful.ustrip(ceil((ylim[2]-ylim[1])/dx)));

    Nx2=2*Nx;
    Ny2=2*Ny;
    dx2=dx/2;
    dy2=dy/2;

    ER2=ones(Nx2,Ny2);
    UR2=copy(ER2);

    x_range=range(0.5,Nx-0.5)*dx;
    y_range=range(0.5,Ny-0.5)*dy;
    x_range=x_range.-mean(x_range);
    y_range=y_range.-mean(y_range);

    x2_range=range(0.5,Nx2-0.5)*dx2;
    y2_range=range(0.5,Ny2-0.5)*dy2;
    x2_range=x2_range.-mean(x2_range);
    y2_range=y2_range.-mean(y2_range);

    Y=[ (j) for  i in x_range,j in y_range];
    X=[ (i) for i in x_range,j in y_range];


    grid_points=[(x,y) for x in x2_range,y in y2_range];


    # Define permittivities for inside and outside the polygon
    epsilon_inside = n_core.^2;
    epsilon_outside = n_cladding.^2;
    pl=[ustrip.((x,y)) for (x,y) in device.polygon];

# Assign permittivities based on whether the rotated points are inside or outside the polygon
    for k in CartesianIndices(ER2)
        (i,j)=Tuple(k)
        point=grid_points[i,j];
        ER2[i, j] = (inpolygon(ustrip.(point),pl)==1) ? epsilon_inside : epsilon_outside;
    end

    result=yeeder2d_grid(ER2,UR2)
    ERzz=result[3];
    URxx=result[4];
    URyy=result[5];


    ω=uconvert(u"rad/s",2*pi*c0/λ);
    k₀=2π/λ;
    println("create the FDFD matrix")
    (DEX,DEY,DHX,DHY)=yeeder2d((Nx,Ny),(dx,dx).*k₀);
    
    println("Apply SCPML")
    (DEX,DEY,DHX,DHY)=apply_scpml_TE(DEX,DEY,DHX,DHY,Nx,Ny,(dx,dx);NPML=NPML);
    Dxx=DHX*spdiagm(0=>URyy[:])*DEX;
    Dyy=DHY*spdiagm(0=>URxx[:])*DEY;

    Ae=Dxx+Dyy+spdiagm(0=>real.(ERzz[:]));
    println("Create source vector for TF/SF")
    QQ=zeros(Nx,Ny);

    QQ=X.<=source_position[1][1][1];

    b0=spdiagm(0=>QQ[:])*Ae-Ae*spdiagm(0=>QQ[:]);

    ER=Tensor(ERxx,ERyy,ERzz,[],[],[]);
    UR=Tensor(URxx,URyy,URzz,[],[],[]);
    RES=(dx*k₀,dy*k₀)
    return Ae,b0, ER,UR,QQ, DEX,DEY,DHX,DHY,Nx,Ny,RES,X,Y,x_range,y_range
end

function create_optimization_elements(device)    
    optimization_area=[];
    if !isempty(device.monitor[OPTIMIZATION])
        for (i,optim) in enumerate(device.monitor[OPTIMIZATION])
            push!(optimization_area,
            [(minimum(x for (x,_) in optim),(minimum(y for (_,y) in optim))),
            (maximum(x for (x,_) in optim),(maximum(y for (_,y) in optim)))]);
        end
    end

    source_position=[];
    for source in device.monitor[SOURCE]
        push!(source_position,[
            (source[1][1],source[1][2]),(source[2][1],source[2][2])
        ]);
    end

    measurment_position=[];
    for measurment in device.monitor[MEASURMENT]
        push!(measurment_position,[
            (measurment[1][1],measurment[1][2]),(measurment[2][1],measurment[2][2])
        ]);
    end
    return optimization_area,source_position,measurment_position    
end

"""
    create_sources(W_wg_in,Nx,Ny,dy,position)

    Create the sources for the TF/SF method and the modes for the optimization
    # Arguments
    * `W_wg_in::Float64`: width of the waveguide
    * `X`: discretized x coordinates
    * `Y`: discretized y coordinates
    * `dx::Unitful.Length`: grid spacing in the y direction
    * `dy::Unitful.Length`: grid spacing in the y direction
    * `position`           : position of the source

    # Returns
    * `sources`            : source vector for TF/SF
    * `mode`               : mode vector for optimization
    * `phase`              : phase vector for optimization
    * `eigenvalue`         : eigenvalue of the mode
"""
function create_sources(W_wg,X,Y,dx,dy,position;λ, noshift=false,kwargs...)
    Nx= size(X,1)
    Ny= size(X,2)
    eigenvalue,e2,ycoor=slabmode(W_wg,Ny,dy;λ,kwargs...);
    k₀=2π/λ;
    sources = zeros(Complex{Float32},Nx, Ny, length(position));
    phase   = zeros(Complex{Float32},Nx);
    for (k, source)  in enumerate(position)
        shifted_e2 = zeros(Complex{Float32}, size(e2));
        if noshift
            shifted_e2 = e2;
        else
            center = (source[2] .+ source[1]) ./ 2;
            idx = findmin(abs.(Y[1, :] .- center[2]))[2];
        
            # Calculate the shift required for e2
            local tmp_shift = idx - findmin(abs.(Y[1, :] .- 0μm))[2];
            shifted_e2[max(1, 1+tmp_shift):min(end, end+tmp_shift)] = e2[max(1, 1-tmp_shift):min(end, end-tmp_shift)];
        end
        for i in axes(sources,1)
            sources[i, :, k] = shifted_e2;
        end
    end
    N=size(sources);#
    modes=deepcopy(sources);
    modes=reshape(modes,prod(N[1:2]),N[3]);

    ddx=dx*k₀;
    for i in axes(sources,1)
        phase[i]=exp(-im*eigenvalue*(i-1)*ddx);
    end
    for k in CartesianIndices(sources)
        i, j, l = Tuple(k)
        sources[i, j, l] = sources[i, j, l] .* phase[i];
    end

    sources=reshape(sources,prod(N[1:2]),N[3]);

    return sources,modes,phase,eigenvalue,center
end

function slabmode(W::T,Ny,dy_norm::V; mode= 0, k₀=2π/1.55μm,n_core,n_cladding,n_substrate=nothing,shift=0.0μm) where {T,V<:AbstractFloat}
    if isnothing(n_substrate);
        n_substrate=n_cladding
    end
    
    println("k₀   : $k₀");
    NA=sqrt(n_core^2-n_cladding^2);  
    @debug("k₀: $k₀ NA: $NA W: $W")
    R=k₀*W*NA;  
    M_allowed=Int(floor(2*R/pi));
    println("Numerical Aperture   : $NA")
    println("Normalized Frequency : $R")
    println("Mode Allowed : $M_allowed")
    N_used=Int(min(mode,M_allowed)+1);
    
    Ny_buffer=Int(floor(Ny/2));
    Ny=Ny+2*Ny_buffer;
    Ny2=Int(Ny*2);
    dy2_norm=dy_norm/2;
    dy2_unit=dy2_norm/k₀;
    #ER2=ones(1,Ny2).*n_cladding^2;
    UR2=ones(Ny2);
    Y2_wg=range(1,Ny2,Ny2).*dy2_unit;
    Y2_wg= Y2_wg.-mean(Y2_wg);
    
    ER2= [abs(y.+shift).<= W/2 ? n_core.^2 : (((y.+shift).<-W/2) ? n_substrate.^2 : n_cladding.^2) for y in Y2_wg];
    #ER2=reshape(ER2,1,Ny2);

    # EXTRACT YEE GRID ARRAYS
    ERxx = ER2[1:2:Ny2]; # ERxx=spdiagm(0=> ERxx[:]);
    ERyy = ER2[2:2:Ny2]; #ERyy=spdiagm(0=> ERyy[:]);
    ERzz = ER2[1:2:Ny2]; #ERzz=spdiagm(0=> ERzz[:]);
    URxx = UR2[2:2:Ny2]; #URxx=spdiagm(0=> URxx[:]);
    URyy = UR2[1:2:Ny2]; #URyy=spdiagm(0=> URyy[:]);
    URzz = UR2[2:2:Ny2]; #URzz=spdiagm(0=> URzz[:]);
    
    NEFF,eigenvectors=slabmode(ERzz,URxx,URyy,buffer=0,dy_norm=dy_norm;mode=mode,n_core=n_core,n_cladding=n_cladding,n_substrate=n_substrate)
    return NEFF,eigenvectors[Ny_buffer+1:end-Ny_buffer,:],ERzz[Ny_buffer+1:end-Ny_buffer,:]
    
end

function slabmode(epsilon_zz::Vector{T},mu_xx::Vector{T},mu_yy::Vector{T},buffer::Int,dy_norm::T;
    mode=0, n_core,n_cladding,n_substrate=nothing,k0=1.0) where T<:AbstractFloat
    
    ERzz=zeros(length(epsilon_zz)+2*buffer);
    URxx=zeros(length(mu_xx)+2*buffer);
    URyy=zeros(length(mu_yy)+2*buffer);
    
    @info("ERxx: $(size(ERyy)) URzz: $(size(URzz)) URyy: $(size(URxx))")

    @info("epsilon_zz: $(size(epsilon_zz)) mu_yy: $(size(mu_yy)) mu_xx: $(size(mu_xx))")


    ERzz[buffer+1:end-buffer]=epsilon_zz; ERzz[1:buffer].=ERzz[buffer+1]; ERzz[end-buffer+1:end].=ERzz[end-buffer];
    
    URyy[buffer+1:end-buffer]=mu_yy; URyy[1:buffer].=URyy[buffer+1]; URyy[end-buffer+1:end].=URyy[end-buffer];
    URxx[buffer+1:end-buffer]=mu_xx; URxx[1:buffer].=URxx[buffer+1]; URxx[end-buffer+1:end].=URxx[end-buffer];

    @info("ERxx: $(size(ERzz)) URzz: $(size(URyy)) URyy: $(size(URxx))")

    Ny=length(ERzz);
    NS  = (1, Ny);
    RES = (oftype(dy_norm,1.0),dy_norm);
    (DEX,DEY,DHX,DHY) = yeeder2d(NS,RES);
    @info("DEX: $(size(DEX)) DEY: $(size(DEY)) DHX: $(size(DHX)) DHY: $(size(DHY)) URzz: $(size(URzz)) ERzz")
    A = -(DHY/spdiagm(0=>URxx)*DEY + spdiagm(0=>ERzz.*k0^2));
    B = inv(Matrix(spdiagm(0=>URyy)));
    ev      = -n_core.^2;
    eigenvalues, eigenvectors = eigs(Matrix(A), B,which=:SR, check=2) #N_used+1,  which=:SR);
    
    D       = [sqrt(Complex{T}(x)) for x in eigenvalues];
    NEFF    = -im*diagm(D);
    
    return real(NEFF)[:],eigenvectors[buffer+1:end-buffer,mode+1] 

end