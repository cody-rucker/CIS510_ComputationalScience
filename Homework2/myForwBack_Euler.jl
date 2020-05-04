using LinearAlgebra
using Plots


function exact(t)
    x = zeros(2)
    x[1] = -2*exp(-2t)*(cos(8t) + 8sin(8t)) + (5/8)*exp(-2t)*(8*cos(8t) - sin(8t))
    x[2] = -10*exp(-2t)*cos(8t) - (5/8)*exp(-2t)*5*sin(8t)

    return x
end


function my_forward_Euler!(Y, Δt, t1, tf, A, y1)

    y .= y1      # initial guess
    Y[:,1] = y[:]

    Exact = Matrix{Float64}(undef,2,N+1)
    Exact[:,1] = y[:]

    for n = 2:N+1
        y[:] = y[:] + Δt*A*y[:]
        Y[:,n] = y[:]
        Exact[:,n] = exact(t[n])
    end

    return (t, Exact)

end


function my_backward_Euler!(Y, Δt, t1, tf, A, y1)

    y .= y1      # initial guess
    Y[:,1] = y[:]

    Exact = Matrix{Float64}(undef,2,N+1)
    Exact[:,1] = y[:]

    Id = Matrix{Float64}(I,N+1,N+1)

    for n = 2:N+1
        y[:] = (I - Δt*A)\y[:]
        Y[:,n] = y[:]
        Exact[:,n] = exact(t[n])
    end

    return (t, Exact)

end

A = [-3 13;-5 -1]
y1 = Array{Float64}(undef,2,1)
y1 = [3;-10]
Δt = 0.01
t1 = 0
tf = 3
t = t1:Δt:tf
N = Integer((tf - t1)/Δt)  # total number of points is N + 1

Y = Matrix{Float64}(undef,2,N+1)
(t, Exact) = my_backward_Euler!(Y, Δt, t1, tf, A, y1)
#(t, Exact) = my_forward_Euler!(Y, Δt, t1, tf, A, y1)

labels = ["approx_y1" "exact_y1" "approx_y2" "exact_y2"]
markershapes = [:circle :star :circle :star]
markercolors = [:green :green :purple :purple]
strd = 1


#plot(t[1:strd:end],[Y[1,1:strd:end], Exact[1,1:strd:end], Y[2,1:strd:end], Exact[2,1:strd:end]], label = labels, shape = markershapes, color = markercolors)
plot(Y[1,1:strd:end],Y[2,1:strd:end])
