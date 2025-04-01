# CO2 Storage Numerical Modeling - SPE11 Implementation

This repository contains the implementation of the 11th Society of Petroleum Engineers Comparative Solution Project (SPE11 CSP) for CO2 storage simulation using DuMuX. The implementation focuses on Versions B and C of the benchmark, which represent field-scale simulations of CO2 injection and long-term storage in heterogeneous geological formations.

## Project Overview

This code was developed as part of a master thesis titled "Influence of Relative Permeability Functions on Numerical Modelling of CO2 Storage Performance" at Clausthal University of Technology. The implementation analyzes how different relative permeability models affect CO2 migration, phase behavior, and trapping mechanisms in geological storage scenarios.

## What's Implemented

- **Version B**: 2D field-scale model with dimensions 8.4km x 1.2km
- **Version C**: 3D field-scale model (8.4km x 5km x 1.2km) with parabolic structure
- **Boundary Volume Multipliers**: Implementation of artificial boundary volumes to mitigate pressure buildup
- **Master of Time Manager**: Dynamic timestep adaptation during injection and post-injection phases
- **Multiple Relative Permeability Models**:
  - Brooks-Corey (Base SPE11 model)
  - Van Genuchten
  - Brooks-Corey with saturated CO2 parameters
  - Brooks-Corey with undersaturated CO2 parameters
- **Shifted Grid Geometry**: Implementation of vertical shift to represent realistic depths
- **Non-isothermal Effects**: Thermal gradients and cooling effects from CO2 injection
- **7-Facies Heterogeneous Model**: Implementation of the SPE11 geology with varying permeability and porosity

## Divergence from SPE11 Specifications

Several adaptations were made to the original SPE11 specifications to ensure numerical stability and computational feasibility:

1. **Pre-Injection Period**: Added a 50-year pre-simulation period before injection starts to allow the system to stabilize
2. **Simulation Duration**: Reduced from 1000 years to 400 years post-injection while maintaining key behavioral insights
3. **Facies 1 Capillary Entry Pressure**: Modified from 193,531 Pa to 200 Pa to avoid numerical convergence issues
4. **Facies 7 Properties**: Implemented with non-zero but extremely low permeability (1e-20 m²) and porosity (1e-4) to avoid solver failures
5. **Maximum Capillary Pressure**: Set to 3e7 Pa as defined in the thesis
6. **No Leakage Modeling**: The implementation does not consider CO2 leakage through faults or caprock

## What's Missing

The implementation has several limitations compared to the full SPE11 specification:

1. **No Version A**: Laboratory-scale model at atmospheric conditions is not implemented
2. **Simplified Well Geometry**: 
   - Wells are implemented as point sources rather than horizontal wells
   - Well 2 is modeled as a straight well instead of following the curved path specified in Version C
3. **Dispersivity and Diffusion Constants**: Not fully incorporated due to time constraints
4. **Salinity Effects**: Implementation uses pure water rather than brine as specified in SPE11

## Code Structure

The repository contains the following key files:

- `main.cc`: Entry point for simulation execution
- `problem.hh`: Definition of the problem setup, boundary conditions, and well injection 
- `properties.hh`: Type tag definitions and property settings
- `spatialparams.hh`: Implementation of geological properties for the seven facies
- `masteroftime.hh`: Dynamic timestep management for different simulation phases
- `boundarycellmanager.cc` and `boundarycellmanager.hh`: Management of artificial boundary volumes
- `materiallaw.hh`: Implementation of different relative permeability models
- `shiftedgridgeometry.hh`: Vertical shifting of the grid geometry to represent correct depths
- `params.input`: Configuration parameters for the simulation

## How to Use

### Prerequisites

- DuMuX (tested with version 3.5)
- DUNE modules
- C++ compiler with C++17 support

### Building and Running

1. Clone the repository
2. Configure using CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   ```
3. Build the executable:
   ```bash
   make
   ```
4. Run the simulation:
   ```bash
   ./spe11 params.input
   ```

## Results and Visualization

The simulation outputs VTK files that can be visualized using ParaView or other VTK-compatible viewers. The output includes:

- Phase distributions (gas saturation, water saturation)
- Component distributions (CO2 mole fractions in gas and liquid phases)
- Pressure fields
- Temperature distributions
- Velocity fields

## License

This code is released under the GPL-3.0 License, consistent with the DuMuX framework licensing.

## Acknowledgments

This implementation is based on the SPE11 CSP specification developed by Jan M. Nordbotten, Martin A. Fernø, Bernd Flemisch, Anthony R. Kovscek, and Knut-Andreas Lie. Special thanks to the DuMuX development team for providing the framework that made this implementation possible.
