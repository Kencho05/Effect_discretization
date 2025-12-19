include("../main.jl")
using Plots
using CSV
using StatsPlots
using LaTeXStrings

# # Parametros
l = 15.91e-3    # H
r = 25          # ohm
c = 50e-6       # F
v0 = 24         # V
rl = 0.05       # ohm

# Matrices del sistema
global Q = [1/l 0; 0 1/c]
global Js(x,u) = [0 -1; 1 0]
global Rs(x) = [rl 0; 0 1/r]
global Gs(x) = [v0; 0]

dHs(x) = Q*x
dt = 1/45e3

# Construcción del sistema PCHS
f = (x,u) -> (Js(x,u)-Rs(x))*dHs(x) + Gs(x).*ones(2,1)*u

# Punto de equilibrio
xref = [0.48; 12]

# Simulación
u_input(x,t) = [0.25 * (sign(t - 0.01) + 1)]

# Ejecución del método
df5 = exact(f, u_input, [0;0], dt, 2000)
# Ploteo de resultados
pl_5 = plot(df5.t, (df5.x1 .* 1/l .- xref[1]), label="Δx1", title="State Variables")
plot!(df5.t, (df5.x2 .* 1/c .- xref[2]), label="Δx2")

# Mostrar el Hamiltoniano en el tiempo (fig A. y fig B.)
u = df5.u1 
x = [df5.x1 df5.x2]' .- (xref .* [l; c])
# Hallar H 
H = [0.5 * x[:,i]' * Q * x[:,i] for i in 1:size(x,2)]
pl_H = plot(df5.t, H, label="H(x)", title="Hamiltonian function", xlabel="Time (s)")

# Subplots 
pl = plot(pl_5, pl_H, layout=(2,1), size=(700,500))
display(pl)