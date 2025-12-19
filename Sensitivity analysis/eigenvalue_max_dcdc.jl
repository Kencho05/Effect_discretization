include("../main.jl")
using Plots

function get_system(type)
    # Parametros
    l = 15.91e-3    # H
    r = 25          # ohm
    c = 50e-6       # F
    v0 = 24         # V
    rl = 0.05       # ohm

    if type == "buck"
        Q  = [1/l 0; 0 1/c]
        Js = (x,u) -> [0 -1; 1 0]
        Rs = [rl 0; 0 1/r]
        Gs = [v0; 0]

    elseif type == "boost"
        Q  = [1/l 0; 0 1/c]
        Js = (x,u) -> [0 -(1-u); (1-u) 0]
        Rs = [rl 0; 0 1/r]
        Gs = [0; 0]

    elseif type == "buck_boost"
        Q  = [1/l 0; 0 1/c]
        Js = (x,u) -> [0 (1-u); -(1-u) 0]
        Rs = [rl 0; 0 1/r]
        Gs = [v0; 0]

    elseif type == "ni_buck_boost"
        Q  = [1/l 0; 0 1/c]
        Js = (x,u) -> [0 -(1-u); (1-u) 0]
        Rs = [rl 0; 0 1/r]
        Gs = [v0; 0]

    else
        error("Tipo de sistema no reconocido")
    end

    return Q, Js, Rs, Gs
end

# Definiciones 
converters = ["buck", "boost", "buck_boost", "ni_buck_boost"]
u_all = 0:0.01:1
plt = plot()
dt = 1/45e3

# Ciclo para cada convertidor
for conv in converters
    Q, Jfrk4, R, G = get_system(conv)
    λmin = []

    for u in u_all
        J = Jfrk4(1,u)
        M = (J - R)*Q
        P = I(2) + dt*M + (dt^2/2)*M*M + (dt^3/6)*M*M*M + (dt^4/24)*M*M*M*M
        Chek = P'*Q*P - Q
        λs = eigvals(Chek)
        Lmin = maximum(λs)
        push!(λmin, Lmin)
    end 

    if conv == "buck"
        global plt = plot(u_all, λmin, label="buck Converter")
    else
        plot!(u_all, λmin, label="$(conv) Converter", xlabel="u", ylabel="Maximum λ", title="Maximum λ vs u",legend=:right)
    end
end

display(plt)