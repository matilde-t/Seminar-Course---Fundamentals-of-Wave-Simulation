using Plots


u = 1.0
mu = 0.25
sigma2 = 0.003

q0(x) = @. exp(-(x - mu)^2 / (2 * sigma2)) / sqrt(2 * π * sigma2)

beta(x) = @. 1 - x

q(x, t) = @. exp(-beta(x) * t) * q0(x - u * t)

x = range(0, stop=1, length=100)
dx = 0.02
x_ = 0:dx:1

Q0 = q0.(x_)
b = beta.(x_)

t_end = 0.6
dt = dx

Q = copy(Q0)
Qstar = copy(Q0)
Qr = copy(Q0)

for t in 0:dt:t_end
    for i in 2:length(Q0)
        Qstar[i] = Q0[i] - u * dt * (Q0[i] - Q0[i-1]) / dx
        Qr[i] = Qstar[i] - b[i] * dt * Qstar[i]
    end
    Q0 .= Qr
end

plot(x, q.(x, t_end), label="Exact")
scatter!(x_, Qr, label="AB")

QAB = copy(Qr)

Q0 .= q0.(x_)
Q .= copy(Q0)
Qstar .= copy(Q0)
Qr .= copy(Q0)

for t in 0:dt:t_end
    for i in 2:length(Q0)
        Qstar[i] = Q0[i] - b[i] * dt * Q0[i]
        Qr[i] = Qstar[i] - u * dt * (Qstar[i] - Qstar[i-1]) / dx
    end
    Q0 .= Qr
end

scatter!(x_, Qr, label="BA")

QBA = copy(Qr)

xlabel!("x")
ylabel!("q(x,t)")
xlims!(0.7, Inf)
ylims!(0, 7)
title!("Solution, β = 1-x")
plot!(legend=:outerbottom, legendcolumns=3)
png("NoCommute")

error = @. (QAB - QBA) / QAB
ab_error = @. QAB - QBA
scatter(x_, [error ab_error], label=["Relative" "Absolute"])
title!("Error, β = 1-x")
xlabel!("x")
ylabel!("Error")
xlims!(0.7, Inf)
ylims!(-0.01, 0.08)
png("NoCommuteErr")