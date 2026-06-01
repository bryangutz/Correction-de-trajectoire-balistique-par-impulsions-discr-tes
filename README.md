\# Discrete Pulse Ballistic Trajectory Correction



This repository contains MATLAB code developed during a research project at the Centre Automatique et Systèmes (CAS), Mines Paris PSL.



The project studies the correction of a nonlinear ballistic trajectory using discrete impulsive velocity corrections. The approach combines a 6-DOF projectile simulation, trajectory linearization, state transition matrices, sensitivity analysis and regularized least-squares optimization.



\## Project context



The objective was to analyze how small impulsive corrections applied during flight can modify the final impact point of a ballistic projectile.



The work focuses on:



\- nonlinear 6-DOF projectile dynamics;

\- aerodynamic forces and moments;

\- quaternion-based attitude representation;

\- numerical trajectory integration until impact;

\- linearization around a nominal trajectory;

\- computation of the State Transition Matrix (STM);

\- sensitivity analysis of the impact point;

\- impulse-based trajectory correction;

\- reachable set analysis under energy constraints.



\## Main idea



A nominal trajectory is first simulated using a nonlinear dynamical model.



The system is then linearized around this trajectory through the numerical computation of the Jacobian matrix.



The State Transition Matrix is used to propagate small perturbations and estimate how a velocity impulse applied during flight affects the final impact point.



This leads to a correction problem formulated as a regularized least-squares problem.



\## Repository structure



\### `sim\_missile\_6dof.m`



Main script for simulating the nonlinear 6-DOF projectile trajectory.



\### `proj\_dynamics.m`



Nonlinear projectile dynamics: translational motion, attitude dynamics, aerodynamic forces and moments.



\### `impact\_event.m`



Event function used to stop the numerical integration when the projectile reaches the ground.



\### `compute\_A\_numeric.m`



Numerical computation of the Jacobian matrix using finite differences.



\### `augmented\_dynamics.m`



Augmented system used to integrate both the nonlinear state and the State Transition Matrix.



\### `forward\_sensitivity.m`



Forward sensitivity propagation based on the linearized system.



\### `A\_of\_t.m`



Helper function used to evaluate/interpolate the time-varying Jacobian matrix.

