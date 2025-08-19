# Discrete Markov Chain Approximation of an SDE on a 2D Lattice

This project implements a **discrete-space approximation** of a conservative differential equation (SPDE) on a 2D periodic lattice (it becomes an SDE there).
We work on a 2 dimensional regular grid with nearest-neighbor edges and periodic boundary conditions.   

Our aim is to compare three levels of description:

1. The **limiting PDE** (conservative equation) on the grid,  
2. A **particle system** of many independent Markov chains on the grid,  
3. The **scaling limit** of the particle system.  

The PDE can describes the evolution of the probability of finding a particle moving according to a stochastic dynamic, the precise definition of whihch will be given later. 
A key interest is the **entropy-dissipation inequality (EDI)**: how quickly the relative entropy decays under the dynamics and, if the system fluctuates, how does the entropy change relate with the kinetic and potential cost.
In particular, we focus on **non-reversible Markov chains**, where the entropy dissipation inequlaities need to be modified. The practical reason to study non-reversible systems is that are know to often converge faster than reversible ones, see ([Rey-Bellet & Spiliopoulos, 2015](https://iopscience.iop.org/article/10.1088/0951-7715/28/7/2081)).  
Indeed, Chen & Hwang (2013) showed that adding an antisymmetric “circulation” to a reversible chain produces a new chain with uniformly smaller asymptotic variance [[Chen & Hwang, 2013]](https://ideas.repec.org/p/nus/nuswps/1375.html).  

This motivates our use of non-reversible rates to potentially accelerate convergence (e.g. in simulated annealing) and to explore the EDI in this setting.

---

## Method Overview

To study these questions we follow a three-step approach:

- **Grid Discretization**  
  Formulate the limiting conservative PDE on a 2D periodic grid (e.g., a Fokker–Planck equation with drift $-\nabla V(x)$).

- **Particle Approximation**  
  Replace the continuum equation by a system of many independent particles evolving as a continuous-time Markov chain on the grid.  
  Initially we take **no interactions** between particles; future work will add interactions (zero-range, Stein gradient, etc.).

- **Scaling Limit**  
  Analyze the behavior as the number of particles grows, comparing the empirical measure and currents to the PDE solution.  
  In particular, we monitor the entropy dissipation over time.

This pipeline lets us verify numerically that the discrete Markov model converges to the intended limit and satisfies the expected entropy decay inequality, even when detailed balance is broken.

---

## Model and Dynamics

We consider a potential $V(x)$ at each lattice site $x \in \text{Grid}$.  
The transition rate from site $x$ to a neighboring site $y$ is defined by

\[
r(x,y) = \bigl(1 + c(x,y)\bigr)\exp\!\Bigl(-\tfrac{\beta}{2}[V(y) - V(x)]\Bigr),
\]

where:

- $c(x,y)$ is **antisymmetric** ($c(x,y) = -c(y,x)$),
- $c$ is chosen so that the invariant measure remains $\pi(x) \propto e^{-\beta V(x)}$,
- Namely, we impose $\sum_{x,y} e^{-\beta[V(x)+V(y)]/2}\,c(x,y) = 0$.

👉 An easy illustrative case is when $c(\cdot,\cdot)$ is *orthogonal* to $\nabla V$, so that the non-reversible perturbation circulates along level sets of $V$.  
In this setting the stationary law is still $e^{-\beta V(x)}$, but the dynamics **breaks detailed balance**.

---

### Probability Currents

Under these rates, the probability current on each edge $x \to y$ is:

\[
J_{xy} = r(x,y)\,\pi(x) - r(y,x)\,\pi(y).
\]

We decompose:

- Symmetric part: $J^S_{xy} = \tfrac{1}{2}(J_{xy} + J_{yx})$  
- Skew part: $J^A_{xy} = \tfrac{1}{2}(J_{xy} - J_{yx})$

The skew part $J^A$ represents the **net non-reversible flow** and drives entropy production.  

Our code computes these currents and the evolving empirical density.  
We then check the **entropy dissipation inequality (EDI)**, i.e. whether the time-derivative of relative entropy $H(\rho_t|\pi)$ matches the computed entropy production.

---

## Implementation

The repository contains code to:

- **Set up the lattice** (2D periodic grid, nearest-neighbor structure),
- **Define transition rates** $r(x,y)$,
- **Simulate particle trajectories** as continuous-time Markov chains,
- **Visualize results** (densities, currents, entropy dissipation).

### Code Structure

- A **graph structure** represents the lattice and its edges.  
- A **simulation routine** runs many independent particles on the graph.  
- **Python scripts** process the output and produce plots, such as:
  - Empirical density $\rho_t(x)$ over time,
  - Symmetric and antisymmetric edge currents,
  - Relative entropy $H(\rho_t|\pi)$ decay,
  - Entropy dissipation gap: observed $dH/dt$ vs. theoretical bound.

The code is modular:  
each model instance (grid, potential $V$, perturbation $c$, temperature $\beta$) is a reusable object.  
You can easily change $c$ or add interactions in future.

---

## Future Work

Planned extensions:

- **Interacting Particles**  
  Add interactions (e.g. zero-range processes, Stein variational dynamics) and study their effect on convergence and EDI.

- **Rigorous Scaling Limit**  
  Formal proof that the empirical measure converges to the PDE solution and satisfies EDI.

- **Unknown Invariant Measure**  
  Apply the methods to cases where $\pi$ is not known.  
  The framework already estimates $\pi$ empirically from long runs.

- **Non-Reversible Optimization**  
  Test applications to simulated annealing and sampling.  
  Non-reversible dynamics are known to accelerate mixing [[Chen & Hwang, 2013]](https://ideas.repec.org/p/nus/nuswps/1375.html).

---

## References

- Stroock, D. W., & Zheng, W. (1997). *Markov chain approximations to symmetric diffusions*.  
  Annales de l’Institut Henri Poincaré (B) Probability and Statistics, 33(5), 619–649.  
  [Link](https://www.numdam.org/article/AIHPB_1997__33_5_619_0.pdf)

- Chen, Y., & Hwang, C.-R. (2013). *Accelerating reversible Markov chains*.  
  National University of Singapore, Working Paper No. 1375.  
  [Link](https://ideas.repec.org/p/nus/nuswps/1375.html)

- Rey-Bellet, Laurent, & Spiliopoulos, Konstantinos (2015). *Irreversible Langevin samplers and variance reduction: a large deviations approach*. **Nonlinearity**, **28**(7), 2081–2103.  
  [Link](https://iopscience.iop.org/article/10.1088/0951-7715/28/7/2081)

- Liu, Qiang, & Wang, Dilin (2016). *Stein variational gradient descent: A general purpose Bayesian inference algorithm*. In **Advances in Neural Information Processing Systems** (NeurIPS), **29**.  
  [Link](https://proceedings.neurips.cc/paper/2016/hash/45f0b2bf851c9a435b783f0d86fbb4c8-Abstract.html)

---
