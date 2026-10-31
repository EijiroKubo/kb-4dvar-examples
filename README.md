# KB-4DVar Examples

Test simulations of Koopman operator-based four-dimensional variational data assimilation (KB-4DVar).

## Overview

This repository contains example codes for testing KB-4DVar, a data assimilation approach that uses an extended dynamic mode decomposition (EDMD)-based approximate adjoint operator.

The method aims to perform four-dimensional variational data assimilation without explicitly developing tangent linear and adjoint models for the original physical simulation code.

The repository includes two numerical examples:

| Example | Description |
|---|---|
| [Kármán vortex](./kb4dvar_karman_vortex) | Data assimilation for two-dimensional incompressible flow, including initial-state correction using KB-4DVar. |
| [Shallow-water equation](./kb4dvar_shallow_water_equation) | Data assimilation for a two-dimensional shallow-water model, including water-depth observations and correction of state variables. |

## Code and Usage

Each example contains Python and Fortran code for the corresponding numerical simulation and data assimilation procedures.

Please refer to the README files in the individual example directories for implementation details, dependencies, compilation instructions, and execution procedures.

## Notes

These codes are provided as research examples to demonstrate the basic implementation and behavior of KB-4DVar.

The examples are not intended to serve as a complete reproduction package for all numerical results presented in the associated research paper. Some simulation settings, output procedures, and supporting files may differ from those used in the original experiments.

## Third-Party Software

The Kármán vortex simulation includes code based on [NS2D_4DVar](https://github.com/DAE-Code/NS2D_4DVar), developed by Takashi Misaka.

Please refer to the license files and third-party notices included in the corresponding example directories for further information.

## License

See the [LICENSE](./LICENSE) file and any additional license notices in the example directories. Third-party components remain subject to their respective license terms.
