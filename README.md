# Epidemic Control with MPC and Agent-Based Modelling

Bachelor thesis MATLAB repository for modelling and control of a two-city epidemic using:

- a SIRVHDB macro model in MATLAB/Simulink,
- a stochastic agent-based micro model,
- model predictive control for hospital-load regulation,
- ABM-to-ODE calibration for control-oriented plant design.

## Main Entry Points

- `main.m` runs the closed-loop MPC + ABM simulation and compares it with an open-loop baseline.
- `validate_abm.m` compares the ABM ensemble response with the SIRVHDB ODE model.
- `run_visual_abm.m` visualises the ABM dynamics using the same micro-model functions as the main simulation.
- `integration/build_plant_ga.m` calibrates and linearises the ODE plant used by MPC.

## Repository Structure

- `macro/` contains the parameter file, Simulink model, calibrated plant files, and Simulink block images.
- `micro/` contains the ABM initialisation, update, and aggregation functions.
- `integration/` contains ODE and plant-calibration code.
- `figures/` contains selected thesis/result figures.

## Running the Main Simulation

Open MATLAB in this repository root and run:

```matlab
main
```

Required MATLAB toolboxes include Simulink and Model Predictive Control Toolbox. The calibrated plant files `macro/plant.mat` and `macro/plant_ga.mat` are included because they are required by the scripts.

Generated results are written to `integration/results_closed_loop.mat` and are intentionally ignored by Git.
