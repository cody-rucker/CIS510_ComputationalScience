using SparseArrays
using Plots
include("myForwBack_Euler.jl")

#theme(:vibrant)

# TODO: this code is in development and has not been verified!

# Solve steady state and time-dependent heat equation in 1D
# u_t  = k u_xx + F(x,t) on 0 ≦ x ≦ 1
# u(x, 0) = f(x)
# u(0, t) = u(1, t) = 0

# Steady state: 0 = k u_xx + F(x), u(0) = u(1) = 0


function F(x,t)
    return 2*(pi^2 -1)*exp(-2t)*sin.(pi*x)
end

function G(x,t)
    return k*exp(-t)*pi^2*sin.(pi*x)
end

function f(x)
    return sin.(pi*x)
end

function exact(x,t)
    return exp(-2t)*sin.(pi*x)
    #return exp(-pi^2 * t)*sin.(pi*x)
end

function steady_state_solve(k, Δx, my_source, my_exact)
    N  = Integer((1-0)/Δx) # N+1 total nodes, N-1 interior nodes
    x = 0:Δx:1
    # A is N-1 by N-1
    A = (k/Δx^2)*(sparse(1:N-1,1:N-1,-2*ones(N-1), N-1, N-1) + sparse(2:N-1,1:N-2,ones(N-2),N-1,N-1) +
         sparse(1:N-2,2:N-1,ones(N-2),N-1,N-1))

    b = -my_source(x[2:N])

    u = A\b # vector solution at the interior nodes

    U = [0;u;0] # vector solution at all nodes
    E = my_exact(x) # exact solution evaluated at all nodes

    return (U, E)
end

function time_dependent_heat(k, Δx, Δt, T, my_source, my_initial)
    N  = Integer(ceil((1-0)/Δx)) # N+1 total nodes, N-1 interior nodes
    x = 0:Δx:1
    t = 0:Δt:T
    M = Integer(ceil((T-0)/Δt)) # M+1 total temporal nodes

    λ = Δt/Δx^2

    # A is N-1 by N-1

    A = (1-2*λ*k)*sparse(1:N-1,1:N-1,ones(N-1), N-1, N-1) + (λ*k)*sparse(2:N-1,1:N-2,ones(N-2),N-1,N-1) +
         (λ*k)*sparse(1:N-2,2:N-1,ones(N-2),N-1,N-1)




    u = Array{Float64}(undef,N-1)
    u .= my_initial(x[2:N])  # setting the initial condition for interior nodes
    Y = zeros(N+1,M+1)
    #Y[2:N,1] .= u[:]

#=
    for n = 1:M
        b = Δt * my_source(x[2:N],t[n])
        u[:] = A * u[:] + b
        Y[2:N,n] .= u[:]

    end
=#


    my_forward_Euler!(Y, Δt, 0, T, A, F, u, x, N)




    #U = [0;u;0] # vector solution at all nodes
    #E = my_exact(x) # exact solution evaluated at all nodes

    return Y#(t,Y)
end

k = 2
Δx = 0.1

T = 1
λ = 0.5
Δt = λ*Δx^2 / k

x = 0:Δx:1
t = 0:Δt:T
xfine = 0:0.01:1
#=
(U, E) = steady_state_solve(k, Δx, F, exact)
labels = ["approx" "exact"]
markershapes= [:circle :star5];
markercolors= [:green :red];
plot(x,[U E], label = labels, shape = markershapes, color = markercolors)
=#

υ = time_dependent_heat(k, Δx, Δt, T, F, f)

#=
pyplot(legend=false, xlim=(0,1), ylim=(0,1.25))

anim = @animate for i = 1:length(t)
    #plot(x,exact(x[:],t[i]))
    #plot(x[:], [υ[:,i], exact(x[:], t[i])])
    p1 = plot(xfine,exact(xfine[:],t[i]), title="exact")
    p2 = plot(x,υ[:,i], title="forward Euler", shape=:star)
    plot(p1,p2, layout=(1,2))


end

gif(anim, "anim_fps15.gif", fps = 30)
=#
