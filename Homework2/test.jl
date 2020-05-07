using Plots


x = 0:0.05:1
t = 0:0.05:1

u = zeros(length(x), length(t))

for j = 1:length(t)
    for i = 1:length(x)
        u[i,j] = sin(2*pi*(x[i] - t[j]))
    end
end

anim = @animate for i = 1:length(t)
    plot(u[:,i])
end

gif(anim, "anim_fps15.gif", fps = 10)
