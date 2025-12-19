include("../main.jl")
using JuMP
using Plots
import SCS

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

# Ciclo para cada convertidor
for conv in converters
    Q, Jdc, R, G = get_system(conv)
    δ_values = Float64[]
    #println("Analizando convertidor: ", conv)

    for u in u_all
        #println("  u = ", u)
        J = Jdc(1,u)
        model = Model(SCS.Optimizer)

        @variable(model, δ >= 0)

        M = 2*R - δ*(J - R)'*Q*(J - R)

        @constraint(model, M in PSDCone())

        @objective(model, Max, δ)
        optimize!(model)

        if termination_status(model) == MOI.OPTIMAL
            push!(δ_values, value(δ))
        else
            push!(δ_values, 0.0)
        end
    end

    if conv == "buck"
        global plt = plot(u_all, δ_values, label="buck Converter")
    else 
        plot!(plt, u_all, δ_values, label="$(conv) Converter", xlabel="u", ylabel="Maximum δ", title="Maximum δ vs u", legend=:topleft)
    end
end

display(plt)
