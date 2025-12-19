using LinearAlgebra
using DataFrames
using Plots


function get_exact(x,u,dt)
    a = 1
    if length(u) == 1
        a = I(1)
    end
    b = u*a
    xa = (I(length(x)) - dt*(Js(x,u) - Rs(x))*Q) \ (x + dt*Gs(x)*(u*a))
    return xa
end

function explicit_euler(f, uf, x0, dt, nt)
    # Iniciliazaciones
    discreto = DataFrame()
    discreto.t = (1:nt) * dt

    x = x0
    u = uf(x,0)

    xode = zeros(length(x), nt)
    uode = zeros(length(u), nt)

    # Condiciones iniciales
    xode[:, 1] = x
    uode[:, 1] = u
    
    # Simulación
    for k = 2:nt
        # Dato de la entrada en el tiempo
        u = uf(x, (k-1)*dt) # u(x,t) 
        # Discretización 
        x = x + dt * f(x,u)
        # Almacenamiento de resultados
        xode[:, k] = x
        uode[:, k] = u
    end 

    for i in 1:length(x0)
        discreto[!, "x$(i)"] = vec(xode[i, :])
    end
    for j in 1:length(u)
        discreto[!, "u$(j)"] = vec(uode[j, :])
    end

    return discreto
end 

function implicit_euler(f,uf,x0,dt,nt)
    # Iniciliazaciones
    discreto = DataFrame()
    discreto.t = (1:nt) * dt

    x = x0
    u = uf(x,0)

    xode = zeros(length(x), nt)
    uode = zeros(length(u), nt)

    # Condiciones iniciales
    xode[:, 1] = x
    uode[:, 1] = u
    
    # Simulación
    for k = 2:nt
        # Dato predictivo
        xa = get_exact(x,u,dt) 
        # Dato de la entrada en el tiempo
        u = uf(x, (k-1)*dt) # u(x,t) 
        # Discretización 
        xn = x + dt * f(xa,u)  # siguiente paso        
        x = xn
        # Almacenamiento de resultados
        xode[:, k] = x
        uode[:, k] = u
    end 

    for i in 1:length(x0)
        discreto[!, "x$(i)"] = vec(xode[i, :])
    end
    for j in 1:length(u)
        discreto[!, "u$(j)"] = vec(uode[j, :])
    end

    return discreto
end 

function midpoint(f,uf,x0,dt,nt)
    # Iniciliazaciones
    discreto = DataFrame()
    discreto.t = (1:nt) * dt

    x = x0
    u = uf(x,0)

    xode = zeros(length(x), nt)
    uode = zeros(length(u), nt)

    # Condiciones iniciales
    xode[:, 1] = x
    uode[:, 1] = u
    
    # Simulación
    for k = 2:nt
        # Dato predictivo
        xa = get_exact(x,u,dt) 
        # Dato de la entrada en el tiempo
        u = uf(x, (k-1)*dt) # u(x,t) 
        # Discretización 
        xm = x + dt/2 * f(xa,u)  # siguiente paso 
        # Punto medio
        x = 2*xm - x
        # Almacenamiento de resultados
        xode[:, k] = x
        uode[:, k] = u
    end 

    for i in 1:length(x0)
        discreto[!, "x$(i)"] = vec(xode[i, :])
    end
    for j in 1:length(u)
        discreto[!, "u$(j)"] = vec(uode[j, :])
    end

    return discreto
end 

function rk4(f,uf,x0,dt,nt)
    # Iniciliazaciones
    discreto = DataFrame()
    discreto.t = (1:nt) * dt

    x = x0
    u = uf(x,0)

    xode = zeros(length(x), nt)
    uode = zeros(length(u), nt)

    # Condiciones iniciales
    xode[:, 1] = x
    uode[:, 1] = u
    
    # Simulación
    for k = 2:nt
        # Dato de la entrada en el tiempo
        u = uf(x, (k-1)*dt) # u(x,t) 

        # RK4
        k1 = f(x,u)
        k2 = f(x + dt/2 * k1, u)
        k3 = f(x + dt/2 * k2, u)
        k4 = f(x + dt * k3, u)

        x = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

        # Almacenamiento de resultados
        xode[:, k] = x
        uode[:, k] = u
    end 

    for i in 1:length(x0)
        discreto[!, "x$(i)"] = vec(xode[i, :])
    end
    for j in 1:length(u)
        discreto[!, "u$(j)"] = vec(uode[j, :])
    end

    return discreto
end

function exact(f,uf,x0,dt,nt)
    # Iniciliazaciones
    discreto = DataFrame()
    discreto.t = (1:nt) * dt

    x = x0
    u = uf(x,0)

    xode = zeros(length(x), nt)
    uode = zeros(length(u), nt)

    # Condiciones iniciales
    xode[:, 1] = x
    uode[:, 1] = u
    
    # Simulación
    for k = 2:nt
        # Dato de la entrada en el tiempo
        u = uf(x, (k-1)*dt) # u(x,t) 
        # Discretización 
        M(x,u) = (Js(x,u) - Rs(x))*Q
        x = exp(dt * M(x,u)) * x + (M(x,u) \ (exp(dt * M(x,u)) - I(length(x)))*(Gs(x).*ones(2,1)*u))
        # Almacenamiento de resultados
        xode[:, k] = x
        uode[:, k] = u
    end 

    for i in 1:length(x0)
        discreto[!, "x$(i)"] = vec(xode[i, :])
    end
    for j in 1:length(u)
        discreto[!, "u$(j)"] = vec(uode[j, :])
    end

    return discreto
end

