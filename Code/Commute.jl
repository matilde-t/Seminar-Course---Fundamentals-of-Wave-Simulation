using Plots


u = 1.0
mu = 0.25
sigma2 = 0.003

q0(x) = @. exp(-(x - mu)^2 / (2 * sigma2)) / sqrt(2 * π * sigma2)

beta(x) = @. 1 - 0 * x

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

xlims!(0.7, Inf)
ylims!(0, 7)
title!("β = 1")
plot!(legend=:outerbottom, legendcolumns=3)
png("Commute")

error = @. (QAB - QBA) / QAB
scatter(x_, error, legend=false)
title!("Relative Error β = 1")
xlabel!("x")
ylabel!("(Q_AB - Q_BA)/Q_AB")
xlims!(0.7, Inf)
ylims!(-0.05, 0.05)
png("CommuteErr")