using Plots
using LinearAlgebra

function Astep(qi, qi_, dt, dx, u)

    return @. qi - u * dt / dx * (qi - qi_)

end

function Bstep(qi, beta, dt)

    return @. qi - beta * dt * qi

end

u = 1.0
mu = 0.25
sigma2 = 0.003

q0(x) = @. exp(-(x - mu)^2 / (2 * sigma2)) / sqrt(2 * π * sigma2)

beta(x) = @. 1 - x

q(x, t) = @. exp(-beta(x) * t) * q0(x - u * t)

x = range(0, stop=1, length=100)
t_end = 0.6

errors = Any[]

for ex in range(-3, -1, 6)

    dx = 10^ex
    global dt = dx

    local x_ = 0:dx:1

    local Q0 = q0.(x_)
    local b = beta.(x_)


    local Qref = q(x_, t_end)

    for t in 0:dt:t_end
        Qstar = Astep(Q0[2:end], Q0[1:end-1], dt, dx, u)
        Q = Bstep(Qstar, b[2:end], dt)
        global Q0[2:end] = Q
    end

    QgodunovAB = copy(Q0)
    errAB = norm(Qref .- QgodunovAB, Inf)

    plot(x, q.(x, t_end), label="Exact")
    scatter!(x_, QgodunovAB, label="Godunov AB")

    Q0 .= q0.(x_)

    for t in 0:dt:t_end
        Qstar = Bstep(Q0, b, dt)
        Q = Astep(Qstar[2:end], Qstar[1:end-1], dt, dx, u)
        global Q0[2:end] = Q
    end

    QgodunovBA = copy(Q0)
    errBA = norm(Qref .- QgodunovBA, Inf)

    scatter!(x_, QgodunovBA, label="Godunov BA")

    Q0 .= q0.(x_)

    for t in 0:2*dt:t_end

        Qstar = Astep(Q0[2:end], Q0[1:end-1], dt, dx, u)
        global Q0[2:end] = Bstep(Qstar, b[2:end], 2 * dt)
        Q = Astep(Q0[2:end], Q0[1:end-1], dt, dx, u)
        global Q0[2:end] = Q

    end

    Qstrang = copy(Q0)
    errS = norm(Qref .- Qstrang, Inf)

    scatter!(x_, Qstrang, label="Strang")
    xlabel!("x")
    ylabel!("q(x,t)")
    xlims!(0.7, Inf)
    # ylims!(0, 7)
    title!("δt = δx = 10^" * string(ex))
    plot!(legend=:outerbottom, legendcolumns=4)
    png("Accuracy_advection_1e" * string(ex))

    global errors = cat(errors, [dx errAB errBA errS], dims=1)

end

errors = cat(errors, [NaN NaN NaN NaN], dims=1)

plot(errors[:, 1], errors[:, 2], xaxis=:log, yaxis=:log, label="Godunov AB", marker=:circle)
plot!(errors[:, 1], errors[:, 3], xaxis=:log, yaxis=:log, label="Godunov BA", marker=:circle)
plot!(errors[:, 1], errors[:, 4], xaxis=:log, yaxis=:log, label="Strang", marker=:circle)
xlabel!("δt = δx")
ylabel!("max norm of error")
title!("Advection-Reaction Equation")
png("ErrorAdvection")