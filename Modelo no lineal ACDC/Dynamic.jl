include("../main.jl")
using Plots


function abc2dq_mod(mag, ang)
    md = mag * sin(ang * π/180)
    mq = mag * cos(ang * π/180)
    return [md; mq]
end

function u_input(x, t)
    if t < 0.1
        mag = 0.2
    else 
        mag = 0.6
    end 
    ang = 10
    return abc2dq_mod(mag, ang)
end


# # Parametros
f = 50
ω = 2*π*f
l = 383e-6
c = 160e-6
r = 24e-3
ldc = 9.5

# Parametros externos
gd = 0
gq = 380 * sqrt(2/3)

dt = 50e-5

# Matrices del sistema
J0 = [0 -l*ω 1 0 0; l*ω 0 0 1 0; -1 0 0 -c*ω 0; 0 -1 c*ω 0 0; 0 0 0 0 0]
JL(u) = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 u[1]; 0 0 0 0 u[2]; 0 0 -u[1] -u[2] 0]

global Gs(x) = zeros(5,2)
global Q = [1/l 0 0 0 0; 0 1/l 0 0 0; 0 0 1/c 0 0;0 0 0 1/c 0;0 0 0 0 1/ldc]
global Js(x,u) = J0 + JL(u)
global Rs(x) = [r 0 0 0 0; 0 r 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0] 
global d = [-gd; -gq; 0; 0; 0]

# Construcción del sistema PCHS
f = (x,u) -> (Js(x,u)-Rs(x))*(Q*x) + d # Revisar exacta 

# Modelo incremental 
xref = [-3.35;70.5;8.3832;312.35;120]

# Ejecución del método
x0 = [-0.019*l; 15.6907*l; 312.1562*c; 0.3789*c; 120*ldc]
df1 = rk4(f, u_input, x0, dt, 400)
# Ploteo de resultados
pl_1 = plot(df1.t, (df1.x1 * 1/l .- xref[1]), label="Δid")
plot!(df1.t, (df1.x2 * 1/l .- xref[2]), label="Δiq")

pl_2 = plot(df1.t, (df1.x3 * 1/c .- xref[3]), label="Δvd")
plot!(df1.t, (df1.x4 * 1/c .- xref[4]), label="Δvq")

pl_3 = plot(df1.t, (df1.x5 * 1/ldc .- xref[5]), label="Δidc")

display(plot(pl_1, pl_2, pl_3, layout = (3,1), size=(800,600)))