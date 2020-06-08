using LinearAlgebra
using SparseArrays
using Plots
include("../../Finite_Element_Methods/convergence_rates.jl")

function my_ode_solve(Y, t, Δt, N, A)

    u = Y[1:N]
    v = Y[N+1 : 2*N]

    U = zeros(2*N, length(t))

    unew = zeros(N)
    #vnew = zeros(N-2)
    vnew = zeros(N)

    for i = 1:length(t)
        unew .= u .+ (Δt .* v)

        #vnew .= v[2:N-1] .+(Δt .* A*u[2:N-1])
        vnew .= v .+(Δt .* A*u)

        u .= unew
        #v[2:N-1] .= vnew
        v .= vnew

        v[1] = -2*pi*cos(2*pi*t[i])
        v[N] = -2*pi*cos(2*pi*(1-t[i]))

        U[:, i] = [u;v]
    end

    return U

end
function ℜ(Δx)
    #Δx = 0.1
    c = 1
    λ = 0.2
    Δt = λ*(Δx)^2/ c

    x = 0:Δx:1
    t = 0:Δt:1

    N = length(x)

    u = zeros(N)
    v = zeros(N)

    uexact(x,t) = sin(2*pi*(x-t))
    vexact(x,t) = -2*pi*cos(2*pi*(x-t))

    # initial displacement
    u .= uexact.(x,0)

    # initial velocity
    v .= vexact.(x,0)

    Y0 = [u;v]

    #A = (1/Δx^2)*(sparse(1:N-2,1:N-2,-2*ones(N-2), N-2, N-2) + sparse(2:N-2,1:N-3,ones(N-3),N-2,N-2) +
    #     sparse(1:N-3,2:N-2,ones(N-3),N-2,N-2))

    A = (1/Δx^2)*(sparse(1:N,1:N,-2*ones(N), N, N) + sparse(2:N,1:N-1,ones(N-1),N,N) +
              sparse(1:N-1,2:N,ones(N-1),N,N))

    y = my_ode_solve(Y0,t, Δt, N, A)

    uₑ = y[1:N, length(t)] .- uexact.(x, t[length(t)])

    error = sqrt(Δx * uₑ' * uₑ)

    #xfine = 0:Δx/20:1
    return error
end

convergence_rates(ℜ, 0.1, 4)
#stride_time = 10

#pyplot(legend=true)#, size=(500,500))#, xlim=(0,1), ylim=(0,1.25))

#=
for i = 1:stride_time:length(t)
    #IJulia.clear_output(true)

    p = plot(x,y[1:N,i],label = "approx. u", shape = :circle, color = :blue)
    ylims!((-1,1))
    display(p)

    Efine = uexact.(xfine, t[i])
    pexact = plot!(xfine,Efine, label = "exact", color = :red)
    display(pexact)

    sleep(0.25)

end
=#
#=
xfine = 0:0.01:1
pyplot(legend=true, size=(1000,750), xlim=(0,1), ylim=(-1.25,1.25))
#labels = ["forward" "backward"]
#markershapes = [:star :none]
anim = @animate for i = 1: 500#length(t)
    #plot(x,exact(x[:],t[i]))
    #plot(x[:], [U[:,i], exact(x[:], t[i])])
    p1 = plot(xfine,uexact.(xfine,t[i]), label="exact")
    p2 = plot!(x, y[1:N,i], label="approximate")
    #plot(p1,p2, layout=(1,2))

end

gif(anim, "anim_fps30.gif", fps = 30)
=#
