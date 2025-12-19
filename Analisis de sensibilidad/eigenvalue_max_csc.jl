include("../main.jl")
using Plots
using LinearAlgebra

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

dt = 50e-6

# Matrices del sistema
J0 = [0 -l*ω 1 0 0; l*ω 0 0 1 0; -1 0 0 -c*ω 0; 0 -1 c*ω 0 0; 0 0 0 0 0]
JL(u) = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 u[1]; 0 0 0 0 u[2]; 0 0 -u[1] -u[2] 0]

Q = [1/l 0 0 0 0; 0 1/l 0 0 0; 0 0 1/c 0 0;0 0 0 1/c 0;0 0 0 0 1/ldc]
Jf(x,u) = J0 + JL(u)
R = [r 0 0 0 0; 0 r 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0] 

# Rango de u1 y u2
u1 = range(-1, 1, length=15)
u2 = range(-1, 1, length=15)

# Crear matriz 
Z = zeros(length(u2), length(u1))

for i in 1:length(u2)
    for j in 1:length(u1)
        u = [u1[j]; u2[i]]
        J = Jf(1,u)
        M = (J - R)*Q
        P = I(5) + dt*M + (dt^2/2)*M*M + (dt^3/6)*M*M*M + (dt^4/24)*M*M*M*M
        Chek = P'*Q*P - Q
        λs = eigvals(Chek)
        Lmax = maximum(λs)
        Z[i, j] = Lmax
    end
end

display(surface(u1, u2, Z,
    xlabel="md",
    ylabel="mq",
    zlabel="λmax",
    title="CSC Converter Eigenvalue Analysis",
    colorbar = false
))