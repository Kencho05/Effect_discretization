# Effect of discretization

## Preliminaries 

### Port-Controlled Hamiltonian Systems

Port-controlled Hamiltonian systems (pH-system) are a specialized class of dynamical systems that integrate concepts from Hamiltonian mechanics with port-based modeling, allowing for the analysis and control of complex physical systems. A large class of power-electronic converters can be represented by the following state-space model:

$$\dot{x} = \left(J(u) - R\right)\nabla H(x) + Gu , $$
$$ y = G^\top \nabla H(x) , $$

where, $x\in\mathbb{R}^n$ is the state vector, $u\in\mathbb{R}^m$ the input vector, and $y\in\mathbb{R}^m$ the output vector. The matrix $J(u)$ represents interconnection while $R$ represents dissipation.  The Hamiltonian function $H:\mathbb{R}^n\rightarrow\mathbb{R}$, for the applications studied here, is given by a quadratic form:

$$ H(x) = \frac{1}{2}x^\top Q x, $$

where $Q\succ 0$

#### Passivity of PCHS

Passivity is a key property shared by many physical sistems (see [[Van-der-Schaft & Jeltsema 2014]](https://people.math.ethz.ch/~hiptmair/Seminars/PHS_24/VSJ14.pdf) for more details):  

**Definition:** A system $\dot{x} = f(x, u)$, $y = h(x, u)$, where $x \in \mathcal{X}\subseteq \mathbb{R}^n$ and $u, y \in \mathbb{R}^m$, is called *passive* if there exists a differentiable storage function $S : \mathcal{R} \rightarrow \mathbb{R}$ with $S(x) \geq 0, x \in \mathcal{X}$ , satisfying the differential dissipation inequality:

$$\dot{S} \leq u^\top y$$

## Motivation 

Time discretization techniques are required to simulate the dynamics modeled as port-Hamiltonian systems; however, the discretization process may lead to the loss of passivity properties. Therefore, this repository contains multiple examples of power electronic converters to illustrate how passivity can be analytically preserved in the discrete-time domain for different discretization techniques, including Forward Euler, Backward Euler, Midpoint, Exact, and RK4 methods. The repository presents sufficient conditions for passivity preservation and validates them through passivity indices and sensitivity analysis. Additionally, it explains how to obtain optimal parameter values when the discretization method depends on specific tuning parameters.
