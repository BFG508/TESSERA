# T.E.S.S.E.R.A. ☄️
**T**hermal **E**valuation & **S**pherical **S**urface **T**essellation for **R**eentry **A**nalysis

A comprehensive MATLAB-based toolkit designed for the geometric modeling and transient thermal analysis of spacecraft heat shields during atmospheric reentry. This personal project is divided into two main engineering modules: a highly uniform spherical dome tessellation engine and a robust 1D transient heat conduction solver capable of simulating both ablative single-layer materials and complex multi-layer Thermal Protection Systems (TPS).

## 🚀 Features
* **Spherical Dome Tessellation:** Generates quasi-uniform geometric meshes on spherical caps utilizing Fibonacci spirals and golden angle increments.
* **Geodesic Poisson-Disk Filtering:** Implements a rigorous greedy acceptance filter using latitude/longitude grids to prevent polygon overlap and guarantee minimum geodesic separation between tessellation nodes.
* **Transient 1D Heat Conduction (Crank–Nicolson):** Features a semi-implicit, unconditionally stable finite difference solver to compute temperature profiles $T(x,t)$ across the shield thickness over time.
* **Multi-Layer TPS Modeling:** Supports the definition of stacked material layers (e.g., outer ceramic tiles, insulation) with temperature-dependent thermal properties (conductivity $k$, specific heat $c_p$, density $\rho$) using constant, linear, or polynomial models.
* **Ablation Mechanics:** Incorporates a dynamic moving-boundary ablation model for single-layer configurations, automatically calculating material removal (depth and mass loss) and applying the specific enthalpy of ablation once a critical temperature threshold is exceeded. 
* **Advanced Boundary Conditions:** Handles complex environmental interactions at the hot face, including imposed time-varying convective heat fluxes and linearized radiative cooling (&epsilon;&sigma;(T<sup>4</sup> - T<sub>&infin;</sub><sup>4</sup>)), alongside customizable cold-face boundary conditions (Robin convection or insulated).
* **Data-Driven Heatflux Integration:** Capable of reading external reentry heating profiles (e.g., Apollo 4 trajectory data) from CSV/Excel datasets and safely interpolating them (PCHIP) into the thermal solver.

## 🛠️ Technology Stack
Developed entirely in **MATLAB** (R2024b), relying on core matrix operations, implicit linear system solvers, and custom spherical geometry mathematics. No external toolboxes are strictly required.

## 📂 Repository Structure
* `/Functions` - Core mathematical utilities for geometry generation (`fibonacciPoints`, `poissonGeodesicFilter`, `polygonOnSphere`) and exact spherical polygon area evaluation.
* `/Heatflux` - Trajectory and heat flux datasets (e.g., `Apollo4.xlsx`) along with the `loadHeatflux` parsing scripts.
* `/ThermalProtectionSystem` - The primary thermal analysis models:
  * `transientConduction_onelayer` - Single-layer solver featuring the moving-boundary ablation model.
  * `transientConduction_multilayer` & `multilayerTPS` - Multi-layer Crank-Nicolson solvers with temperature-dependent property evaluation.
  * Config generators like `alfaTPS.m` and `betaTPS.m` for easy scenario setup.

## ⚙️ Installation & Usage
1. **Clone the repository:**
   ```bash
   git clone https://github.com/BFG508/TESSERA.git
2. **Open MATLAB** and navigate to the cloned `TESSERA` directory.
3. **Add to Path**: Ensure that all subdirectories (`/Functions`, `/Heatflux`, `/ThermalProtectionSystem`) are added to your MATLAB path to allow the main scripts to access necessary functions and data files. You can do this by right-clicking the `TESSERA` folder in the Current Folder browser and selecting Add to Path > Selected Folders and Subfolders.
4. **Run Geometric Tessellation**: Execute `sphericalDomeTesselation.m` to generate and visualize the 3D polygon distribution over the spacecraft's forward heat shield.
5. **Run Thermal Analysis**: 
* Execute `transientConduction_onelayer.mlx` to evaluate the ablative performance and thickness loss for different initial shield configurations.
* Execute `transientConduction_multilayer.mlx` (which calls `betaTPS` or `alfaTPS`) to evaluate energy balances and interface temperatures in a stacked TPS.