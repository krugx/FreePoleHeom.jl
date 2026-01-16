# Theory: Free-Pole HEOM

This section outlines the theoretical foundation of the Free-Pole Hierarchy of Equations of Motion (FP-HEOM), specifically focusing on its formulation in the rotating frame.

## Open System Formulation

The total Hamiltonian is given by:

$$
\hat{H}_{tot}(t) = \hat{H}_S(t) + \hat{H}_I(t) + \hat{H}_B
$$

With the reservoir (bath) operator:

$$
\hat{H}_B = \sum_k \left( \frac{\hat{p}_k^2}{2m_k} + \frac{1}{2}m_k\omega_k\hat{x}_k^2 \right)
$$

The System + Control Hamiltonian $\hat{H}_S(t)$ is:

$$ 
\hat{H}_S(t) = \frac{\omega_q}{2}\hat{\sigma}_z + \epsilon(t)\hat{\sigma}_x
$$ 

The interaction Hamiltonian is assumed to be:

$$ 
\hat{H}_I(t) = - \hat{V} \sum_k c_k \hat{x}_k = -\hat{V} \hat{X}_B
$$ 

where in the lab frame $\hat{V} = \hat{\sigma}_x$.

### Rotating Frame Transformation

We transform to the rotating frame using the unitary operator:

$$ 
\hat{U}(t) = \exp\left( i\frac{\omega_d}{2}t\hat{\sigma}_z \right)
$$ 

The transformed Hamiltonian $\hat{H}_{rot}(t)$ is:

$$ 
\begin{aligned}
\hat{H}_{rot}(t) &= \hat{U}(t)\hat{H}_{tot}(t)\hat{U}^\dagger(t) + i \dot{\hat{U}}(t)\hat{U}^\dagger(t) \\
&= \hat{U}(t) \hat{H}_S(t) \hat{U}^\dagger(t) - \frac{\omega_d}{2}\hat{\sigma}_z + \hat{H}_I^{rot}(t) + \hat{H}_B
\end{aligned}
$$ 

**System Part:**

$$ 
\begin{aligned}
\hat{H}_{S}^{rot}(t) &= \hat{U}(t) \left( \frac{\omega_q}{2}\hat{\sigma}_z + \epsilon(t)\hat{\sigma}_x \right) \hat{U}^\dagger(t) - \frac{\omega_d}{2}\hat{\sigma}_z \\
&= \frac{\Delta \omega}{2}\hat{\sigma}_z + \epsilon(t) \left( \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t} \right)
\end{aligned}
$$ 

where $\Delta \omega = \omega_q - \omega_d$.
Assuming the drive $\epsilon(t)$ has the form:
$$ 
\epsilon(t) = 2 \epsilon_x(t) \cos(\omega_d t) + 2 \epsilon_y(t) \sin(\omega_d t)
$$ 
Using the Rotating Wave Approximation (RWA) by neglecting terms oscillating at $\pm 2\omega_d$, we obtain:
$$ 
\hat{H}_{S}^{RWA}(t) = \frac{\Delta \omega}{2}\hat{\sigma}_z + \frac{1}{2}\epsilon_x(t)\hat{\sigma}_x + \frac{1}{2}\epsilon_y(t)\hat{\sigma}_y
$$ 

**Interaction Part:**

The interaction operator $\hat{V} = \hat{\sigma}_x$ transforms as:

$$ 
\begin{aligned}
\hat{V}^{rot}(t) &= \hat{U}(t) \hat{\sigma}_x \hat{U}^\dagger(t) \\
&= \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t}
\end{aligned}
$$ 

Thus, the interaction Hamiltonian in the rotating frame is:

$$ 
\hat{H}_I^{rot}(t) = - \left( \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t} \right) \hat{X}_B
$$ 

## HEOM Derivation with RWA

The Feynman-Vernon influence functional is determined by the bath correlation function $C(t) = \langle \hat{X}_B(t) \hat{X}_B(0) \rangle_B$.
In the rotating frame, the effective interaction involves phases $e^{\pm i \omega_d t}$.

The phase of the influence functional is (schematically):
$$ 
\Phi = \int_0^t ds \int_0^s d\tau \; \hat{V}^{rot}(s)^\times \left( C(s-\tau) \hat{V}^{rot}(\tau)^\to - C^*(s-\tau) \hat{V}^{rot}(\tau)^\leftarrow \right)
$$ 
where $\times$ denotes commutator ($A^\times B = [A, B]$), $\to$ denotes left multiplication, and $\leftarrow$ denotes right multiplication.

Substituting $\hat{V}^{rot}(t) = \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t}$:

We expand the products and apply the RWA (neglecting terms like $e^{\pm i \omega_d (s+\tau)}$ which oscillate fast, keeping terms like $e^{\pm i \omega_d (s-\tau)}$).

This splits the correlation function into two effective branches:

1.  **Positive branch ($+$):** Couples $\hat{\sigma}_-$ at $s$ with $\hat{\sigma}_+$ at $\tau$.
    $$ C_+(t) = C(t) e^{-i\omega_d t} $$
    Associated operator: $\hat{q}_+ = \hat{\sigma}_+$
2.  **Negative branch ($-$):** Couples $\hat{\sigma}_+$ at $s$ with $\hat{\sigma}_-$ at $\tau$.
    $$ C_-(t) = C(t) e^{+i\omega_d t} $$
    Associated operator: $\hat{q}_- = \hat{\sigma}_-$

### Correlation Function Expansion

We expand the original bath correlation function as a sum of exponentials (e.g., using FP-HEOM decomposition):

$$ 
C(t) = \sum_{k=1}^{M} d_k e^{-\gamma_k t}
$$ 

Under RWA, we obtain a set of $2M$ effective modes. Let us index them by $j=1,\dots,2M$:

*   For $j = 1, \dots, M$ (derived from $C_+$):
    $$ \tilde{d}_j = d_j, \quad \tilde{\gamma}_j = \gamma_j + i\omega_d, \quad \hat{q}_j = \hat{\sigma}_+ $$
*   For $j = M+1, \dots, 2M$ (derived from $C_-$):
    $$ \tilde{d}_j = d_{j-M}, \quad \tilde{\gamma}_j = \gamma_{j-M} - i\omega_d, \quad \hat{q}_j = \hat{\sigma}_- $$

### Hierarchy of Equations of Motion (HEOM)

We define the Auxiliary Density Operators (ADOs) $\hat{\rho}_{\mathbf{m}, \mathbf{n}}$ where $\mathbf{m} = \{m_1, \dots, m_{2M}\}$ and $\mathbf{n} = \{n_1, \dots, n_{2M}\}$.
Here, $m_j$ tracks the excitation from the "forward" path (Left multiplication) and $n_j$ from the "backward" path (Right multiplication).

The equation of motion is:

$$
\begin{aligned}
\frac{d}{dt}\hat{\rho}_{\mathbf{m,n}} &= -i \left[ \hat{H}_S^{RWA}, \hat{\rho}_{\mathbf{m,n}} \right] 
- \sum_{j=1}^{2M} (m_j \tilde{\gamma}_j + n_j \tilde{\gamma}_j^*) \hat{\rho}_{\mathbf{m,n}} \\
&\quad \text{// Interaction Terms (Creation)} \\
&\quad -i \sum_{j=1}^{2M} \sqrt{(m_j+1)\tilde{d}_j} \; \hat{q}_j \hat{\rho}_{\mathbf{m}_j^+, \mathbf{n}}  \quad (\text{Left coupling})\\
&\quad +i \sum_{j=1}^{2M} \sqrt{(n_j+1)\tilde{d}_j^*} \; \hat{\rho}_{\mathbf{m}, \mathbf{n}_j^+} \hat{q}_j  \quad (\text{Right coupling})\\
&\quad \text{// Interaction Terms (Annihilation)} \\
&\quad -i \sum_{j=1}^{2M} \left( \sqrt{m_j \tilde{d}_j} \left[ \hat{q}_j, \hat{\rho}_{\mathbf{m}_j^-, \mathbf{n}} \right] + \sqrt{n_j \tilde{d}_j^*} \left[ \hat{q}_j, \hat{\rho}_{\mathbf{m}, \mathbf{n}_j^-} \right] \right)
\end{aligned}
$$
