using Plots
using LinearAlgebra

function Astep(q, dt)

    return @. q + dt * q

end

function Bstep(q, dt)

    return @. q - dt * q * q

end

q0 = 0.1

q(t) = @. 1 / (9 * exp(-t) + 1)

t_end = 8

errors = Any[]

for dt in [1e-2 1e-1 1]

    t_ = 0:dt:t_end
    local Qref = q(t_)

    local Q0 = q0

    QAB = copy(Qref)

    i = 1
    for t in dt:dt:t_end
        Qstar = Astep(Q0, dt)
        Q_ = Bstep(Qstar, dt)
        QAB[i+1] = Q_
        Q0 = Q_
        i = i + 1
    end

    errAB = norm(Qref .- QAB, Inf)

    plot(t_, Qref, label="Exact")
    scatter!(t_, QAB, label="Godunov AB", markerstrokewidth=.1)

    Q0 = q0

    QBA = copy(Qref)

    i = 1
    for t in dt:dt:t_end
        Qstar = Bstep(Q0, dt)
        Q_ = Astep(Qstar, dt)
        QBA[i+1] = Q_
        Q0 = Q_
        i = i + 1
    end

    errBA = norm(Qref .- QBA, Inf)

    scatter!(t_, QBA, label="Godunov BA", markerstrokewidth=.1)

    Q0 = q0

    QS = copy(Qref)

    i = 1
    for t in dt:dt:t_end
        Qstar = Astep(Q0, 0.5 * dt)
        Q0 = Bstep(Qstar, dt)
        Q_ = Astep(Q0, 0.5 * dt)
        QS[i+1] = Q_
        Q0 = Q_
        i = i + 1
    end

    errS = norm(Qref .- QS, Inf)

    scatter!(t_, QS, label="Strang", markerstrokewidth=.1)
    xlabel!("t")
    ylabel!("u(t)")
    title!("δt = " * string(dt))
    plot!(legend=:outerbottom, legendcolumns=4)
    png("Accuracy_logistic_" * string(dt))

    global errors = cat(errors, [dt errAB errBA errS], dims=1)

end

errors = cat(errors, [NaN NaN NaN NaN], dims=1)

plot(errors[:, 1], errors[:, 2], xaxis=:log, yaxis=:log, label="Godunov AB", marker=:circle)
plot!(errors[:, 1], errors[:, 3], xaxis=:log, yaxis=:log, label="Godunov BA", marker=:circle)
plot!(errors[:, 1], errors[:, 4], xaxis=:log, yaxis=:log, label="Strang", marker=:circle)
xlabel!("δt")
ylabel!("max norm of error")
title!("Logistic Equation")
png("ErrorLogistic")