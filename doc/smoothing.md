# Elliptic Smoothing

The implemented smoothing is based on [^1], Chapter 4: Elliptic Generation
Systems. The discretization and solution method is described in 4.2.2 and
applied for a mesh global smoothing. In the
following the used assembling of the linear solution system is derived and documented.

[^1]: J. F. Thompson, B. K. Soni, and N. P. Weatherill, Eds., Handbook of grid generation. Boca Raton, Fla: CRC Press, 1999.

## Finite Differences

See https://en.wikipedia.org/wiki/Finite_difference

$$\vec{x}_{\xi} \approx \frac{1}{2} \left( \vec{x}_{i+1,j} - \vec{x}_{i-1,j}
\right)$$

$$\vec{x}_{\eta} \approx \frac{1}{2} \left( \vec{x}_{i,j+1} - \vec{x}_{i,j-1}
\right)$$

$$\vec{x}_{\xi\xi} \approx \left( \vec{x}_{i+1,j} - 2 \vec{x}_{i,j} +
\vec{x}_{i-1,j} \right)$$

$$\vec{x}_{\eta\eta} \approx \left( \vec{x}_{i,j+1} - 2 \vec{x}_{i,j} +
\vec{x}_{i,j-1} \right)$$

$$\vec{x}_{\xi\eta} \approx \frac{1}{4} \left( \vec{x}_{i+1,j+1} -
\vec{x}_{i+1,j-1} - \vec{x}_{i-1,j+1} + \vec{x}_{i-1,j-1} \right)$$

## Linear system of equations

$$0 = P \vec{x}_{\xi\xi} - 2Q \vec{x}_{\xi\eta} + R \vec{x}_{\eta\eta} +
S \vec{x}_{\xi} + T \vec{x}_{\eta}$$
where
$$P = \vec{x}_{\eta} \cdot \vec{x}_{\eta}$$
$$Q = \vec{x}_{\xi} \cdot \vec{x}_{\eta}$$
$$R = \vec{x}_{\xi} \cdot \vec{x}_{\xi}$$
$$S=PP_{11}^1 - 2QP_{12}^1 + RP_{22}^1$$
$$T=PP_{11}^2 - 2QP_{12}^2 + RP_{22}^2$$

For Laplace smoothing: $S=0$ and $T=0$.

Introducing the FD approximations and rearranging yields the 9 point stencil for
every point $(i,j)$:
$$
\begin{align}
0 = & && \\
  & \left( -2P-2R \right)           && \vec{x}_{i,j} \\
  & \left( P + \frac{1}{2}S \right) && \vec{x}_{i+1,j} \\
  & \left( P - \frac{1}{2}S \right) && \vec{x}_{i-1,j} \\
  & \left( R + \frac{1}{2}T \right) && \vec{x}_{i,j+1} \\
  & \left( R - \frac{1}{2}T \right) && \vec{x}_{i,j-1} \\
  & \left( -\frac{1}{2}Q \right)    && \vec{x}_{i+1,j+1} \\
  & \left( \frac{1}{2}Q \right)     && \vec{x}_{i+1,j-1} \\
  & \left( \frac{1}{2}Q \right)     && \vec{x}_{i-1,j+1} \\
  & \left( -\frac{1}{2}Q \right)    && \vec{x}_{i-1,j-1} \\
\end{align}
$$