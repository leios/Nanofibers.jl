using SpecialFunctions, MAT

# Set of parameters to generate fields
mutable struct Param
    resx::Int64
    resy::Int64
    resz::Int64
    mode_a::Int64
    mode_b::Int64

    xmax::Float64
    ymax::Float64
    zmax::Float64
    lambda::Float64
    n1::Float64
    n2::Float64
    N1::Float64
    N2::Float64
    beta::Float64
    h::Float64
    q::Float64
    a::Float64    
    s::Float64
    A::Float64
end

# struct to hold all fields
struct Fields
    Er::Array{1,Float64}
    Ephi::Array{1,Float64}
    Ez::Array{1,Float64}
    Ex::Array{1,Float64}
    Ey::Array{1,Float64}
end

# This function relies on "packeddata.mat" and will create realiztic parameters
# based on the provided file for experimentally realistic results
function read_mat(value)

    return params
end

function derive(f, l, arg, dx)
    return (f(l, arg) - f(l, arg+dx))/dx
end

function calc_s!(params)
    h = params.h
    q = params.q
    a = params.a
    l = params.mode_a
    dx = a/100
    s = (1/(h*h*a*a))+(1/(q*q*a*a))
        /((derive(besselj, l, h*a, dx)/(h*a*besselj(l,h*a)))
         +(derive(besselk, l, q*a, dx)/(q*a*besselk(l,q*a))))
    params.s = s
end

function calc_A!(params)
    beta = params.beta
    q = params.q
    h = params.h
    n1 = params.n1
    n2 = params.n2
    N1 = params.N1
    N2 = params.N2
    l = params.mode_a
    a = params.a

    A = (beta/(2*q))
        *((besselj(l,h*a)/besselk(l,q*a))
          / sqrt(2*pi*a*a(n1*n1*N1 + n2*n2*N2)))

    params.A = A
end

function calc_N1!(params)
    beta = params.beta
    h = params.h
    s = params.s
    l = params.mode_a
    a = params.a
    N1 = (beta*beta / (4*h*h))
         *(((1 - s)^2)*(besselj(l-1,h*a)^2 + besselj(l, h*a)^2)
            +((1 + s)^2)*(besselj(l+1,h*a)^2-besselj(l, h*a)*besselj(l+2,h*a)))
         + 0.5*(besselj(l, h*a)^2-besselj(l-1,h*a)*besselj(l+1,h*a))

    params.N1 = N1
end

function calc_N2!(params)
    beta = params.beta
    h = params.h
    s = params.s
    q = params.q
    l = params.mode_a
    a = params.a

    N2 = (besselj(l,h*a)^2/(2*besselk(l,q*a))^2)
         *((beta*beta / (4*q*q))
          *(((1 - s)^2)*(besselk(l-1,q*a)^2 - besselk(l, q*a)^2)
             +((1+s)^2)*(besselk(l+1,q*a)^2-besselk(l, q*a)*besselk(l+2,q*a)))
          +0.5*(besselk(l, q*a)^2+besselk(l-1,q*a)*besselk(l+1,q*a)))
    params.N2 = N2
end

function calc_E_circ(params, x, y, z)
    q = params.q
    A = params.A
    l = params.mode_a
    s = params.s

    beta = params.beta
    lambda = params.lambda

    r = sqrt(x*x + y*y)

    Er = im*A*((1-s)*besselk(l-1,q*r)+(1+s)*besselk(l+1,q*r))
    Ephi = -A*((1-s)*besselk(l-1,q*r)-(1+s)*besselk(l+1,q*r))
    Ez = 2*A*(q/beta)*besselk(l,q*r)
end

function calc_E_lin(params, x, y, z)
    q = params.q
    A = params.A
    l = params.mode_a
    s = params.s

    phi = atan2(y,x)

    beta = params.beta
    lambda = params.lambda

    r = sqrt(x*x + y*y)

    Ex = sqrt(2)*A*((1-s)besselk(l-1,q*r)cos(angle)
                    +(1+s)besselk(l+1,q*r)cos(2*phi-angle))
    Ey = sqrt(2)*A*((1-s)besselk(l-1,q*r)sin(angle)
                    +(1+s)besselk(l+1,q*r)sin(2*phi-angle))
    Ez = 2*sqrt(2)*im*A*(q/beta)besselk(l,q*r)cos(phi-angle)

end

function generate_fields(params, polarization)
    calc_s!(params)
    calc_N1!(params)
    calc_N2!(params)
    calc_A!(params)

    for i = 1:params.resx
        x = -params.xmax/2 + i*params.xmax/params.resx
        for j = 1:params.resy
            y = -params.ymax/2 + i*params.ymax/params.resy
            for k = 1:params.resz
                z = -params.zmax/2 + i*params.zmax/params.resz
            end
        end
    end
    return Fields(Er, Ephi, Ez, zeros(res), zeros(res), zeros(res))
end

function output_field(field)
end
