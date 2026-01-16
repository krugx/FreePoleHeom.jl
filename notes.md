### Open System Formulation

$$
\hat{H}_{tot}(t) = \hat{H}(t) + \hat{H}_I(t) + \hat{H}_R
$$

With the reservoir operator

$$
\hat{H}_R = \sum_j \left( \frac{\hat{p}_j^2}{2m_j} + \frac{1}{2}m_j\omega_j\hat{x}_j^2 \right)
$$

Where $\hat{H}(t)$ is the System + Control Hamiltonian

$$
\hat{H}(t) = \frac{\omega_q}{2}\hat{\sigma}_z + \epsilon(t)\hat{\sigma}_x
$$

The interaction Hamiltonian

$$
\hat{H}_I(t) = - \hat{V}(t)\sum_j c_j\hat{x}_j = -\hat{V}(t)\hat{X}
$$

Transfer this to the rotating frame using the unitary

$$
\hat{R}(t) = \exp\left( i\frac{\omega_d}{2}t\hat{\sigma}_z \right) \otimes \hat{\mathbb{1}}_R
$$

Resulting in 

$$
\begin{aligned}
\hat{H}_{tot}^{rot}(t) &= \hat{R}(t)\hat{H}_{tot}(t)\hat{R}^\dagger(t) + i \dot{\hat{R}}(t)\hat{R}^\dagger(t)\\
&= \frac{\Delta \omega}{2}\hat{\sigma}_z + \frac{1}{2} \left( \epsilon_x(t)\hat{\sigma}_x - \epsilon_y(t)\hat{\sigma}_y \right)\\
&- \exp\left( i\frac{\omega_d}{2}t\hat{\sigma}_z \right) \hat{V}(t) \exp\left( -\frac{\omega_d}{2}t\hat{\sigma}_z \right)\\
&+ \sum_j \left(\frac{\hat{p}_j^2}{2m_j} + \frac{1}{2}m_j\omega_j \hat{x}_j^2 \right)
\end{aligned}
$$

In the lab frame

$$
\hat{V}(t) = \hat{\sigma}_x
$$

In the rotating frame

$$
\begin{aligned}
\hat{F}(t) &= \exp\left( i\frac{\omega_d}{2}t\hat{\sigma}_z \right) \hat{\sigma}_x \exp\left( -i\frac{\omega_d}{2}t\hat{\sigma}_z \right)\\
&= \cos(\omega_d t) \hat{\sigma}_x - \sin(\omega_d t)\hat{\sigma}_y\\
&=\hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t}
\end{aligned}
$$

Where $\hat{\sigma}_+ = (\hat{\sigma}_x + i\hat{\sigma}_y)/2$ and $\hat{\sigma}_- = (\hat{\sigma}_x - i \hat{\sigma}_y)/2$

And

$$
\begin{aligned}
\epsilon(t) &= \cos(\omega_d t)\epsilon_x(t) + \sin(\omega_d t)\epsilon_y(t)\\
&= \frac{1}{2}\left( e^{i\omega_d t} + e^{-i\omega_d t} \right)\epsilon_x(t) - \frac{i}{2}\left( e^{i\omega_d t} - e^{-i\omega_d t} \right) \epsilon_y(t)
\end{aligned}
$$

Then neglect fast rotating terms $\exp\left( \pm 2i \omega_d t \right)$

Thus we find

$$
\epsilon(t)e^{i\omega_d t} = \frac{1}{2}\left( \epsilon_x(t) - i \epsilon_y(t) \right)
$$

As a consequence

$$
\epsilon(t)\left( \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-i\omega_d t} \right) = \frac{1}{2}\epsilon_x(t)\hat{\sigma}_x - \frac{1}{2}\epsilon_y(t)\hat{\sigma}_y(t)
$$

### HEOM derivation

Feynman-Vernon influence functional

$$
\mathcal{F}_{FV} = \exp\left( -\int_0^t ds (F(s) - F(s'))\int_0^s d\tau \left(C(s-\tau)F(\tau) - C^*(s-\tau)F'(\tau) \right) \right)
$$

Where the matrix representation of $F(t)$:

$$
\hat{F}(t) = \hat{\sigma}_+ e^{i\omega_d t} + \hat{\sigma}_- e^{-\omega_d t}
$$

Represent as $F(s) = q_+(s)e^{i\omega_d s} + q_-(s)e^{-i\omega_d s}$

Neglect $e^{\pm i \omega_d (s+\tau)}$ and we find

$$
\mathcal{F}_{FV} = \mathcal{F}_\downarrow \cdot \mathcal{F}_{\uparrow}
$$

Define as well

$$
\begin{aligned}
C_+(s-\tau) &= C(s-\tau)e^{-i\omega_d (s-\tau)}\\
C_-(s-\tau) &= C(s-\tau)e^{+i\omega_d (s-\tau)}\\
\end{aligned}
$$

Where

$$
\mathcal{F}_\uparrow = \exp\left( -\int_0^t ds (q_-(s) - q_-'(s)) \right)\int_0^t d\tau \left( C_+(s-\tau)q_+(\tau) - C_+^*(s-\tau)q_+'(\tau) \right)
$$

and

$$
\mathcal{F}_\downarrow = \exp\left( -\int_0^t ds (q_+(s) - q_+'(s)) \right)\int_0^t d\tau \left( C_-(s-\tau)q_-(\tau) - C_-^*(s-\tau)q_-'(\tau) \right)
$$

Take the derivative

$$
\begin{aligned}
\frac{d\mathcal{F}_{FV}}{dt} &= \frac{d\mathcal{F}_\uparrow}{dt}\mathcal{F}_\downarrow + \mathcal{F}_\uparrow \frac{d \mathcal{F}_\downarrow}{dt}\\
&= -\left( q_-(t) - q_-'(t) \right) \int_0^t d\tau \left( C_+(t-\tau)q_+(\tau) - C_+^*(t-\tau)q_+'(\tau) \right)\\
&= -\left( q_+(t) - q_+'(t) \right) \int_0^t d\tau \left( C_-(t-\tau)q_-(\tau) - C_-^*(t-\tau)q_-'(\tau) \right)\\
\end{aligned}
$$

Make the expansion

$$
C(t) = \sum_{j=1}^{M} d_j e^{-\gamma_j t}
$$

Where $d_j,\gamma_j \in \mathbb{C}$.

Define the ADOs

$$
\begin{aligned}
\rho_{\bf m,n}(q_f,q_f',t) &= \int dq_i\int dq_i' \int \mathcal{D}[q]\int \mathcal{D}[q']\, \exp\left( i \mathcal{S}[q] - i \mathcal{S}[q'] \right)\\
&\cdot \prod_{j=1}^{2K} \left[ \int_0^t ds \tilde{d}_j e^{-\tilde{\gamma}_j (t-s)} q_j(s) \right]^{m_j}\\
&\cdot \prod_{j=1}^{2K} \left[ \int_0^t ds \tilde{d}_j^* e^{-\tilde{\gamma}_j^* (t-s)} q_j'(s) \right]^{n_j}\\
&\mathcal{F}_{FV} \, \rho(q_i,q_i',0)
\end{aligned}
$$

Where we define the new fit parameters for the multi exponential

- $j=1,\dots,K$: 
  $$
  \begin{aligned}
  \tilde{d}_j &= d_j\\
  \tilde{\gamma}_j &= \gamma_j - i\omega_d\\
  q_j &= q_+
  \end{aligned}
  $$

- $j=K+1,\dots,2K$: 
  $$
  \begin{aligned}
  \tilde{d}_j &= d_j\\
  \tilde{\gamma}_j &= \gamma_j + i\omega_d\\
  q_j &= q_-
  \end{aligned}
  $$

The new ADOs have extended multi-indeces compared to the standard FP-HEOM ${\bf (m,n)} = (m_1,\dots, m_{2K},n_1,\dots,n_{2K})$

Resulting in a HEOM for the rotating frame

$$
\begin{aligned}
\frac{d}{dt}\hat{\rho}_{\bf m,n} &= -\left[\hat{H}_{tot}^{tot}, \hat{\rho}_{\bf m,n}\right] - \sum_{j=1}^{2K} \left( m_j \tilde{\gamma}_j + n_j \tilde{\gamma}_j^* \right) \hat{\rho}_{\bf m,n}\\
&-i \sum_{j=1}^{2K} \sqrt{(m_j +1)\tilde{d}_j} \left[ \hat{q}_j, \hat{\rho}_{\bf m_j^+,n} \right]\\
&-i \sum_{j=1}^{2K} \sqrt{(n_j +1)\tilde{d}_j^*} \left[ \hat{q}_j, \hat{\rho}_{\bf m,n_j^+} \right]\\
&-i \sum_{j=1}^{2K}\left( \sqrt{m_j \tilde{d}_j} \hat{q}_j \hat{\rho}_{\bf m_j^-,n} - \sqrt{n_j \tilde{d}_j^*} \hat{\rho}_{\bf m,n_j^-}\hat{q}_j \right)
\end{aligned}
$$
