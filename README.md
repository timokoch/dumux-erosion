Erosion toy model
------------------

Similar to Derr et al. (2020) Phys. Rev. Lett., [10.1103/PhysRevLett.125.158002](https://doi.org/10.1103/PhysRevLett.125.158002).

Governing equations are

```math
\begin{align}
\kappa \frac{\partial p}{\partial t} - \nabla \cdot \left( K(\phi) \nabla p \right) &= 0, \\
\frac{\partial \phi}{\partial t} - \nabla \cdot \left( D \nabla \phi \right) &= -\phi \mathrm{max}\left\lbrace 0, \nabla p : \nabla p - \sigma^2(g) \right\rbrace, \\
R \frac{\partial g}{\partial t} - \nabla \cdot \left( B \nabla g \right) &= \phi - g,
\end{align}
```

where $\phi$ is the solid volume fraction, $K(\phi) = \frac{(1-\phi)^3}{\phi^2}$,
$\sigma^2(g) = \frac{H(g)-H(0)}{H(1)-H(0)}$ with $H(g) = \frac{1}{2}\left[1 -  \mathrm{tanh}(\omega (g^* - g))\right]$.
We erode when $\nabla p : \nabla p > \sigma^2(g)$, where $\sigma^2(g)$ is a yield stress that depends on the solidity $g$.
