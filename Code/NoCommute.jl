using Plots

u = 1
beta(x) = 1 .- x

mu = 0.25
sigma2 = 0.003
q0(x) = exp.(-(x .- mu) .^ 2 / (2 * sigma2)) / sqrt(2 * pi * sigma2)

q(x, t) = exp.(-beta(x) .* t) .* q0(x .- u * t)

x = range(0, 1, 100)
dx = 0.01
x_ = 0:dx:1
Q0 = q0(x_)
b = beta(x_)
t_end = 0.3
dt = dx / 10

Q = similar(Q0)
Qstar = similar(Q0)
Qr = similar(Q0)

Q[1] = Q0[1]
Qstar[1] = Q0[1]
Qr[1] = Q0[1]

for t = 0:dt:t_end
    for i = 1:(length(Q0)-1)
        Qstar[i+1] = Q0[i+1] - u * dt * (Q0[i+1] - Q0[i]) / dx
        Qr[i+1] = Qstar[i+1] - b[i+1] * dt * Qstar[i+1]
    end
    global Q0 = Qr
end

plot(x, q(x, t_end), label="Exact")
title!("β = 1-x")
scatter!(x_, Qr, label="AB")

Q0 = q0(x_)
Q = similar(Q0)
Qstar = similar(Q0)
Qr = similar(Q0)

Q[1] = Q0[1]
Qstar[1] = Q0[1]
Qr[1] = Q0[1]

for t = 0:dt:t_end
    for i = 1:(length(Q0)-1)
        Qstar[i+1] = Q0[i+1] - b[i+1] * dt * Q0[i+1]
        Qr[i+1] = Qstar[i+1] - u * dt * (Qstar[i+1] - Qstar[i]) / dx
    end
    global Q0 = Qr
end

scatter!(x_, Qr, label="BA")

Q0 = q0(x_)
Q = similar(Q0)
Q[1] = Q0[1]

for t = 0:dt:t_end
    for i = 1:(length(Q0)-1)
        Q[i+1] = Q0[i+1] - u * dt / dx * (Q0[i+1] - Q0[i]) - dt * b[i] * Q0[i+1]
    end
    global Q0 = Q
end

scatter!(x_, Q, label="Unsplit")

xaxis!("x")
yaxis!("q(x,t)")
ylims!(0, 8)
plot!(legend=:outerbottom, legendcolumns=4)
png("NoCommute")