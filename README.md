# SKOOPI: Spectral KOopman Approach to Optimal Control Using Path-Integral Formulation

This paper presents a novel framework for computing optimal control by exploiting the spectral properties of the Koopman operator associated with the uncontrolled dynamical system. The innovative aspect lies in the ability to compute optimal control using information extracted from the open-loop system's behavior. A key advantage of this approach is that it does not require global knowledge of Koopman eigenfunctions. Instead, the optimal control at a given state is determined using the gradient (local information) of the Koopman eigenfunction. The eigenfunction values are computed using a path-integral formula, which is particularly appealing as it only requires open-loop trajectory information starting from the state of interest. By relying solely on local trajectory data, the approach offers a scalable and efficient solution for optimal control computation without explicit knowledge of the system dynamics. 

## Example 1: Stabilization

The figure below illustrates the stabilization performance of the proposed eigenfunction-based controller (SKOOPI). The stabilization performance is evaluated over 100 initial conditions, uniformly sampled within a unit disk centered at $\bx_0 = (-2,\, 3)$, and simulated for $5\,\mathrm{s}$ at $20\,\mathrm{Hz}$. Note that all trajectories converged to the equilibrium point at the origin.

<p align="center">
<img src="Figures/ACC/example1_SKOOPI.png" width="500">
</p>

For comparison, we consider a finite-time LQR controller (FT-LQR) based on linearization at the origin. The same cost matrices are used for both the baseline and the proposed method. The figure below shows that the FT-LQR failed to stabilize the system at the origin. This can be attributed to the fact that the FT-LQR controller relies on the linearization of the system at the origin, which leads to a smaller domain of attraction compared to the eigenfunction-based controller (SKOOPI).

<p align="center">
<img src="Figures/ACC/example1_LQR.png" width="500">
</p>

## Example 2: Trajectory Tracking

Case 1: We initialize the system close to the stable equilibrium point with the end effector position at $\bx_0 = (0,\;0)$. We sample 10 random configurations close to $x_0$. The figure below shows the snapshots of the system tracking the reference trajectory. The eigenfunction-based controller (SKOOPI) is able to track the reference with an RMSE error of $0.0848 \pm 0.0177$ while the finite-time LQR controller (FT-LQR) leads to an RMSE of $0.0801 \pm 0.0204$. Note that both methods are able to track the reference trajectory when initialized close to the stable equilibrium point. In the next case, we explore the performance of each controller when the system is initialized farther away from the stable equilibrium point.

<p align="center">
<img src="Figures/ACC/example2_comparison.png" width="500">
</p>

Case 2: We initialize the system with the position of the end effector at $x_0 = (1.5,\;0)$. The Figure below shows the snapshots of the system tracking the reference trajectory (green) with the proposed eigenfunction-based controller (SKOOPI) and the finite-time LQR controller (FT-LQR). Note that the proposed SKOOPI controller tracks the reference more closely compared to the FT-LQR. The total state cost $\sum\left(\left(x-x^\star)^\top \Q (x-x^\star) \right)\right)$ for SKOOPI is $757.8$ which is much lower than $72463.3$ for FT-LQR (an improvement of $98.95\%$ for this case).

<p align="center">
<img src="Figures/ACC/example2_comparison2.png" width="500">
</p>
