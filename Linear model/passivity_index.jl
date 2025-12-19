include("../main.jl")
using Plots
using StatsPlots
using LaTeXStrings

# Solo aplica para modelo lineal
function get_y(dt,x,u,method)
    # Pasar de función a matriz
    G = Gs(1)
    R = Rs(1)
    J = Js(1,1)
    M = (J - R) * Q
    if method == "explicit_euler"
        x = x[:,1:end-1]
        #y = dt*G'*Q*x - dt^2*G'*Q*R*Q*x + dt^2/2*G'*Q*G*u[2:end]'
        Fe = dt*G
        Pe = I(2) + dt*M
        y = Fe' * Q * Pe * x + 1/2 * Fe' * Q * Fe * u[2:end]'
        return y
    end
    if method == "implicit_euler"
        y = dt*G'*Q * x[:,2:end] # x_{k+1}
        return y
    end
    if method == "midpoint"
        xm = (x[:,1:end-1] .+ x[:,2:end])/2 # (x_k + x_{k+1})/2
        y = dt*G'*Q * xm
        return y
    end
    if method == "rk4"
        P = I(2) + dt*M + (dt^2/2)*M*M + (dt^3/6)*M*M*M + (dt^4/24)*M*M*M*M
        y = (inv(M)*(P - I(2)) * G)' * Q * (P * x[:,1:end-1] + 1/2*inv(M)*(P - I(2)) * G * u[2:end]') # x_{k+1}
        return y
    end
    if method == "exact"
        x = x[:,1:end-1]
        F = M \ (exp(dt*M) - I(2)) * G
        y = F'*Q* exp(dt*M) * x + 1/2 * F'*Q*F * u[2:end]'
        return y
    end
    return y
end

# Simulación
u_input(x,t) = [0.25 * (sign(t - 0.01) + 1)]

# Parametos
nt = 2000
dt = collect(range(7.06e-4, 7.12e-4, step=1e-7))
meths = [explicit_euler, implicit_euler, midpoint, rk4, exact]

# Crear el DataFrame vacío
data = DataFrame(
    dt = Float64[],
    method = String[],
    value = Float64[]
)

for i in dt
    for (k, solver) in enumerate(meths)
        df = solver(f, u_input, [0;0], i, nt)
        x = [df.x1 df.x2]'
        u = df.u1
        # Hallar H 
        H = [0.5 * x[:,j]' * Q * x[:,j] for j in 1:size(x,2)]
        # Calcular el cambio de energía
        ΔH = H[2:end] .- H[1:end-1]
        # Calcular la entrada y 
        y = get_y(i,x,u,string(solver))
        # Restricción pasiva
        s = y' .* u[2:end]
        # y'u - (H+ - H) >= 0  (passvity check)
        balance = s .- ΔH   
        # tomar el infimo
        infimo = minimum(balance)

        # Añadir filas al dataframe
        append!(data, DataFrame(
            dt = i,
            method = string(solver),
            value = infimo
        ))
    end
end 

plt = plot()
# Recorra subs dataframes
for g in groupby(data, :method)
    plot!(plt,g.dt,g.value, label=g.method[1], xlabel="dt", ylabel="Passivity Index ", title="Passivity check", legend=:left)
end
display(plt)