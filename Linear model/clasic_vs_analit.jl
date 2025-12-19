include("../main.jl")
using Plots

# Simulación
u_input(x,t) = [0.25 * (sign(t - 0.01) + 1)]
dt = 1/5e3

# Ejecución del método
df5 = exact(f, u_input, [0;0], dt, 200)

# Datos de los df1
df = df5
u = df.u1 
x = [df.x1 df.x2]'
# Hallar H 
H = [0.5 * x[:,i]' * Q * x[:,i] for i in 1:size(x,2)]
# Calcular el cambio de energía
ΔH = H[2:end] .- H[1:end-1]

# Pasar de función a matriz
G = Gs(1)
R = Rs(1)
x = x[:,1:end-1]

J = Js(1,1)
M = (J - R) * Q
F = M \ (exp(dt*M) - I(2)) * G

# Cálculo y analítica
y = F'*Q* exp(dt*M) * x + 1/2 * F'*Q*F * u[2:end]'         # Exacta
# Restricción pasiva
s = y' .* u[2:end]

# Calculo y clasica 
y2 = dt*G'*Q*x
s2 = y2' .* u[2:end]

scatter(s, ΔH, label="y analitic", legend=:topleft, grid=true)
scatter!(s2, ΔH, label="y classic", legend=:topleft, grid=true)
plot!([0,0.0019],[0,0.0019], linestyle=:dash, color=:orange,label="Passivity Bound") 