# A Predictive Load Flow Control in Contingency-Affected Power Systems Using FRT Enhanced GUPFC and MPC

This repository contains the complete MATLAB source code for a Final Year Project investigating the performance of a Fault-Ride-Through (FRT) enhanced Generalized Unified Power Flow Controller (GUPFC) with Model Predictive Control (MPC) under predictive (time-varying) load conditions.

## How to Run the Simulation
The simulation pipeline consists of 11 integrated scripts. They must be run in the following order:

1. `load_data.m` - Loads IEEE 30-bus network data.
2. `build_Ybus.m` - Constructs the bus admittance matrix.
3. `run_NR.m` - Runs the base-case Newton-Raphson load flow.
4. `contingency.m` - Models the critical N-1 line outage (Branch 5).
5. `add_gupfc.m` - Inserts the GUPFC at buses (3, 7, 17).
6. `mpc_gupfc.m` - Runs the 24-hour load-aware MPC simulation.
7. `frt_test.m` - Evaluates FRT performance across varying load levels.
8. `frt_ksweep.m` - Performs the FRT gain (K-factor) sensitivity sweep.
9. `generate_plots.m` - Generates the final result figures.

*(Note: `run_gupfc_lf.m` and `mpc_cost_fn.m` are helper functions called automatically by the main scripts).*

## Requirements
- MATLAB (Tested on R2018a and R2024b)
- Optimization Toolbox (for `fmincon`)

## Author
**Nwokemodo Martins 210403505** 
