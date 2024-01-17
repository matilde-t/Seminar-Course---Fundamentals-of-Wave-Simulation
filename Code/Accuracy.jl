using Plots

u = 1.0
mu = 0.25
sigma2 = 0.003

q0(x) = @. exp(-(x - mu)^2 / (2 * sigma2)) / sqrt(2 * π * sigma2)

β(x) = @. 1 - x

q(x, t) = @. exp(-β(x) * t) * q0(x - u * t)

x = range(0, stop=1, length=100)
t_end = 0.5

dx_ = [0.1, 0.01, 0.001]
errAB = zeros(0)
errBA = zeros(0)
errStrang = zeros(0)

for dx in dx_
    x_ = 0:dx:1

    Q0 = q0(x_)
    b = β(x_)

    dt = dx

    Q = copy(Q0)
    Qstar = copy(Q0)
    Qr = copy(Q0)

    for t in 0:dt:t_end
        for i in 2:length(Q0)
            Qstar[i] = Q0[i] - u * dt * (Q0[i] - Q0[i-1]) / dx
            Qr[i] = Qstar[i] - b[i] * dt * Qstar[i]
        end
        Q0 = copy(Qr)
    end

    plot(x, q(x, t_end), label="Exact")
    scatter!(x_, Qr, label="Godunov AB")

    QAB = copy(Qr)

    Q0 .= q0(x_)
    Q = copy(Q0)
    Qstar = copy(Q0)
    Qr = copy(Q0)

    for t in 0:dt:t_end
        for i in 2:length(Q0)
            Qstar[i] = Q0[i] - b[i] * dt * Q0[i]
            Qr[i] = Qstar[i] - u * dt * (Qstar[i] - Qstar[i-1]) / dx
        end
        Q0 = copy(Qr)
    end

    scatter!(x_, Qr, label="Godunov BA")

    QBA = copy(Qr)

    Q0 .= q0(x_)
    Q = copy(Q0)
    Qstar = copy(Q0)
    Qr = copy(Q0)

    i = 0
    for t in 0:dt:t_end
        if i % 2 == 0
            for i in 2:length(Q0)
                Qstar[i] = Q0[i] - b[i] * dt * Q0[i]
                Qr[i] = Qstar[i] - u * dt * (Qstar[i] - Qstar[i-1]) / dx
            end
        else
            for i in 2:length(Q0)
                Qstar[i] = Q0[i] - b[i] * dt * Q0[i]
                Qr[i] = Qstar[i] - u * dt * (Qstar[i] - Qstar[i-1]) / dx
            end
        end
        i += 1
        Q0 = copy(Qr)
    end

    scatter!(x_, Qr, label="Strang")

    QStrang = copy(Qr)

    xlabel!("x")
    ylabel!("q(x,t)")
    xlims!(0.5, 1)
    ylims!(0, 7)
    title!("β = 1-x")
    plot!(legend=:outerbottom, legendcolumns=4)
    png("Plot"*string(dx))

    Qreal = q(x_, t_end)

    push!(errAB, copy(norm(Qreal .- QAB, 2)))
    push!(errBA, copy(norm(Qreal .- QBA, 2)))
    push!(errStrang, norm(Qreal .- QStrang, 2))
end

scatter(log10.(dx_), log10.([errAB'; errBA'; errStrang']), label=["AB" "BA" "Strang"])
xlabel!("log_{10}(dx)")
ylabel!("log_{10}(Error)")
plot!(legend=:outerbottom, legendcolumns=3)