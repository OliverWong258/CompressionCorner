# Project Description

This project implements a numerical simulation with flux splitting schemes for solving hyperbolic partial differential equations in MATLAB. The program is designed to run in the MATLAB console with the function `main(Nx, Ny, splitting_type)`.

## Usage

To run the program, open MATLAB and enter the following command in the console:

```matlab
main(Nx, Ny, splitting_type)
```

Where:
- `Nx` is the number of grid points in the x-direction.
- `Ny` is the number of grid points in the y-direction.
- `splitting_type` specifies the flux splitting scheme to use. It can be one of the following:
  - `'LF'` for Lax-Friedrichs scheme.
  - `'SW'` for Steger-Warming scheme.
  - `'VL'` for Van Leer scheme.

### Example
```matlab
main(100, 100, 'LF')
```

This will run the simulation with 100 grid points in both the x and y directions, using the Lax-Friedrichs flux splitting scheme.

## Output

During the simulation, the current calculation time (not the program runtime) will be displayed in the MATLAB console. After the simulation finishes, the program will display the total time it took for the execution.