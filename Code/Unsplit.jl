using Plots

u = 1
beta(x) = 1

mu = 0.25
sigma2 = 0.003
q0(x) = exp.(-(x .- mu) .^ 2 / (2 * sigma2)) / sqrt(2 * pi * sigma2)

q(x, t) = exp.(-beta(x) .* t) .* q0(x .- u * t)

x = range(0, 1, 100)
dx = 0.01
x_ = 0:dx:1
Q0 = q0(x_)
b = similar(Q0)
fill!(b, 1)
t_end = 0.3
dt = 0.0005

Q = similar(Q0)
Q[1] = Q0[1]

for t = 0:dt:t_end
    for i = 1:(length(Q0)-1)
        Q[i+1] = Q0[i+1] - u * dt / dx * (Q0[i+1] - Q0[i]) - dt * b[i] * Q0[i+1]
    end
    global Q0 = Q
end

plot(x, q(x, t_end), label="Exact")
title!("β = 1")
scatter!(x_, Q, label="Unsplit")
xaxis!("x")
yaxis!("q(x,t)")
ylims!(0, 8)
plot!(legend=:outerbottom, legendcolumns=3)
png("Unsplit")
