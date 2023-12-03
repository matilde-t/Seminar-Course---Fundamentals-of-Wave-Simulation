using Plots

u = 1

mu = 0.25
sigma2 = 0.003
q0(x) = @. exp.(-(x .- mu) .^ 2 / (2 * sigma2)) / sqrt(2 * pi * sigma2)

q(x, t) = @. exp.(-t) .* q0(x .- u * t)

x = range(0, 1, 100)

plot(x, [q(x, 0) q(x, 0.25) q(x, 0.5)], label=["t=0" "t=0.25" "t=0.5"])
title!("Advection-Reaction Equation - Exact Solution")
xaxis!("x")
yaxis!("q(x,t)")
plot!(legend=:outerbottom, legendcolumns=3)
png("Advection")
