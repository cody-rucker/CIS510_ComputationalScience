using SparseArrays
using Plots
include("myForwBack_Euler.jl")
include("/home/cody/github/Finite_Element_Methods/convergence_rates.jl")

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

function time_dependent_heat(k, Δx, Δt, T, my_source, my_initial, my_exact)
    N  = Integer(ceil((1-0)/Δx)) # N+1 total nodes, N-1 interior nodes
    x = 0:Δx:1
    t = 0:Δt:T
    M = Integer(ceil((T-0)/Δt)) # M+1 total temporal nodes

    λ = Δt/Δx^2

    # A is N-1 by N-1
    A = (1-2*λ*k)*sparse(1:N-1,1:N-1,ones(N-1), N-1, N-1) +
            (λ*k)*sparse(2:N-1,1:N-2,ones(N-2),N-1,N-1) +
            (λ*k)*sparse(1:N-2,2:N-1,ones(N-2),N-1,N-1)

    A2 = (1+2*λ*k)*sparse(1:N-1,1:N-1,ones(N-1), N-1, N-1) -
             (λ*k)*sparse(2:N-1,1:N-2,ones(N-2),N-1,N-1) -
             (λ*k)*sparse(1:N-2,2:N-1,ones(N-2),N-1,N-1)

    u = Array{Float64}(undef,N-1)
    u .= my_initial(x[2:N])  # setting the initial condition for interior nodes
    u2 = copy(u)
    Y = zeros(N+1,M)         # initializing solution matrix
    Z = copy(Y)

    my_forward_Euler!(Y, Δt, t, A, F, u, x, N, M)
    my_backward_Euler!(Z, Δt, t, A2, F, u2, x, N, M)


    uₑ₁ = my_exact(x[:], t[M]) - Y[:, M]
    uₑ₂ = my_exact(x[:], t[M-1]) - Z[:, M-1]

    e1 = sqrt(Δx * uₑ₁' * uₑ₁)
    e2 = sqrt(Δx * uₑ₂' * uₑ₂)

    return Y, Z, e1, e2
end

# repackage the solution as a function of grid spacing which returns the discrete
# L²-error
function ℜ(Δx)
    k = 2
    #Δx = 0.1

    T = 2
    λ = 0.2
    Δt = λ*Δx^2 / k

    x = 0:Δx:1
    t = 0:Δt:T

    U1, U2, e1, e2 = time_dependent_heat(k, Δx, Δt, T, F, f, exact)

#=
    xfine = 0:0.01:1
    plot(legend=true, size=(500,500), xlim=(0,1), ylim=(0,1.25))
    labels = ["forward" "backward"]
    markershapes = [:star :none]
    anim = @animate for i = 1:length(t)
        #plot(x,exact(x[:],t[i]))
        #plot(x[:], [U[:,i], exact(x[:], t[i])])
        p1 = plot(xfine,exact(xfine[:],t[i]), title="exact")
        p2 = plot(x, [U1[:,i],U2[:,i]], label=labels, shape=markershapes)
        plot(p1,p2, layout=(1,2))

    end

    gif(anim, "anim_fps15.gif", fps = 30)
=#
#=
    xfine = 0:Δx/20:1

    stride_time = 10

    plot(legend=true)#, size=(500,500))#, xlim=(0,1), ylim=(0,1.25))

    for i = 1:stride_time:length(t)
        #IJulia.clear_output(true)

        p = plot(x,U1[:,i],label = "forward", shape = :circle, color = :blue)
        ylims!((0,1))
        display(p)

        p2 = plot!(x,U2[:,i],label = "backward", shape = :star, color = :green)
        ylims!((0,1))
        display(p2)

        Efine = exact(xfine, t[i])
        pexact = plot!(xfine,Efine, label = "exact", color = :red)
        display(pexact)

        sleep(0.25)

    end=#
    return e1#  U1, U2, e1, e2
end
convergence_rates(ℜ, 0.1, 4)
#U1, U2, e1, e2 = ℜ(0.05);
