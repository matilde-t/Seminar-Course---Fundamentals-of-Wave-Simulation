using Plots


beta(x) = 1 .- x
q0(x) = exp.(-0.5 * (x .- 0.25) .^ 2 / (2 * 0.003)) ./ sqrt(2 * pi * 0.003)
u = 1

q(x, t) = exp.(-beta(x) .* t) .* q0(x .- u * t)


x = range(0, 1, 100)

plot(x, [q(x, 0) q(x, 0.25) q(x, 0.5)], label=["t=0" "t=0.25" "t=0.5"])
title!("Advection-Reaction Equation - Exact Solution")
xaxis!("x")
yaxis!("q(x,t)")
plot!(legend=:outerbottom, legendcolumns=3)
png("Advection")
