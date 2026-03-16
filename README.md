# Aeroservoelasticity Research Projects ✈️🚁

[![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Academic](https://img.shields.io/badge/University-Politecnico%20di%20Milano-blue.svg)](https://www.polimi.it/)

A comprehensive study of aeroelastic phenomena in aircraft and rotorcraft systems, developed during the **Aeroservoelasticity** course at Politecnico di Milano (A.A. 2025-2026). This repository explores structural dynamics, stability boundaries, and nonlinear oscillations.

---

## 🔬 Core Projects

### 1. Fixed-Wing Stability (Goland Wing)
Systematic aeroelastic study focusing on structural dynamics and stability boundaries.
* **Structural FEM:** Bending (cubic Hermite) and torsion (linear Lagrange) discretization for a cantilever wing.
* **Flutter Analysis:** Boundary determination using the **p-k method** and **Theodorsen’s unsteady theory**.
* **Control Surfaces:** Evaluation of aileron effects on flutter speed, control reversal, and static divergence.
* **Optimization:** Mass balancing strategies using tip weights to expand flight envelopes.

### 2. Rotary-Wing Aeroelasticity (Puma Blade)
Stability analysis of a helicopter blade in hover, focusing on the pitch-flap flutter phenomenon.
* **Dynamics:** Accounting for **centrifugal stiffening** and inertial coupling between bending and torsion.
* **Aerodynamic Refinement:** Comparative analysis between quasi-steady models and exact spanwise-varying unsteady aerodynamics.
* **Stability Maps:** Parametric sweeps of pitch link stiffness and chordwise CG positions to define safe operational regimes.

### 3. Nonlinear LCO Analysis
Investigation of **Limit Cycle Oscillations (LCO)** induced by actuator free-play nonlinearities in ailerons.
* **Quasi-linearization:** Implementation of the **Sinusoidal Input Describing Function (SIDF)**.
* **Numerical Continuation:** Newton-Raphson iterations to track subcritical LCO branches and frequency hardening.
* **Stability:** Static stability analysis of the LCO based on damping-amplitude variations.

---

## 🛠 Tech Stack & Methods
* **Platform:** MATLAB.
* **Numerical Techniques:** Richardson Extrapolation for modal convergence, Rayleigh-Ritz modal synthesis, and Galerkin projection.
* **Aerodynamics:** Theodorsen theory with Prandtl-Glauert compressibility corrections.

---

## 📂 Repository Structure
* `/01_Goland_wing_analysis`: Fixed-wing FEM and flutter scripts.
* `/02_Puma_blade_stability`: Rotor blade hover stability and maps.
* `/03_Nonlinear_lco_wing`: Nonlinear actuator dynamics and LCO solvers.

## 👨‍💻 Author
**Alejandro Rivera Míguez** *M.Sc. Aeronautical Engineering, Politecnico di Milano* *Professor: Giuseppe Quaranta*

---
*Disclaimer: This repository is intended for academic and research purposes in the field of aeroservoelasticity.*
