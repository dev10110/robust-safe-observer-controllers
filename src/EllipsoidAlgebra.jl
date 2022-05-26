module EllipsoidAlgebra

using Plots, LinearAlgebra

"""
E(x, P) = { z : (z-x) P ( z-x) <= 1
"""
struct Ellipsoid{Tx, TP} 
    x0::Tx
    P::TP
end

Ellipsoid(x0, P, V) = Ellipsoid(x0, P/V)

function plotEllipse!(e::Ellipsoid; kwargs...)

    @assert length(e.x0) == 2
    
    Pinv_half = real.(sqrt(inv(e.P)))
    
    x(θ) = e.x0[1] + first(Pinv_half * [cos(θ), sin(θ)])
    y(θ) = e.x0[2] + last(Pinv_half * [cos(θ), sin(θ)])
    
    plot!(x, y, 0, 2π; kwargs...)
    
end

function project(e::Ellipsoid, dims=[1,2])

    N = length(e.x0)
    Pinv_half = sqrt(inv(e.P))
    
    # construct projection matrix
    T = zeros(2, N)
    T[1,dims[1]] = 1
    T[2, dims[2]] = 1
    
    # get new shape matrix
    TPinvhalf = T * Pinv_half
    PPinv = TPinvhalf * TPinvhalf'
    Pp = real.(inv(PPinv))

    # get new center
    xp0 = T * e.x0
    

    return Ellipsoid(xp0, Pp)

end



function sampleRand(e::Ellipsoid, N=1)
    
    @assert isposdef(e.P)

    Pinv_half = sqrt(inv(e.P))
    
    # sample random from unit sphere, and map to ellipsoid
     y = 2*rand(length(e.x0), N) .- 1
     for i=1:N
        while norm(y[:,i]) > 1.0
            y[:,i] = 2*rand(length(e.x0), 1) .- 1
        end

    end
    
    xs = similar(y)
    for i=1:N
        xs[:,i] = e.x0 + Pinv_half * y[:,i]
    end
    
    return xs
end


function creategrid(d::Integer, n::Integer)

    @assert d >= 1 ("d (number of dimensions) must be a positive integer")
    @assert n >= 2 ("n (number of points) must be a at least 2")

    r = range(-1, 1, length = n)

    iter = Iterators.product((r for _ in 1:d)...)

    return vec([collect(i) for i in iter if dot(i, i) <= 1])
end

const sphere_6_grid = creategrid(6, 5)

function sampleUniformGrid(e::Ellipsoid, n=5)
    N = length(e.x0)
    Pinv_half = sqrt(inv(e.P))
    if N == 6
        pts = [e.x0 + Pinv_half * i for i in sphere_6_grid]
        return pts
    else
        pts = [e.x0 + Pinv_half * i for i in creategrid(N, n)]
        return pts
    end
end

end