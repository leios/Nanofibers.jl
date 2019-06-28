using SpecialFunctions, MAT

# Set of parameters to generate
struct Param
    resx::Int64
    resy::Int64
    resz::Int64
    mode_a::Int64
    mode_b::Int64
    lambda::Float64
    n1::Float64
    n2::Float64
    beta::Float64
    h::Float64
    q::Float64
    
end

# This function relies on "packeddata.mat" and will create realiztic parameters
# based on the provided file for experimentally realistic results
function read_mat(value)

    return params
end

function calc_s(x, y, z, params)
end

function calc_A(x, y, z, params)
end

function calc_N1(x, y, z, params)
    
end

function calc_N2(x, y, z, params)
end

function generate_field()

    return field
end

function output_field(field)
end
