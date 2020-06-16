# circuit-simulations
Modeling the behavior of the voltage and current for various transmission line configurations.

This repository contains MATLAB simulations of voltage(V) and current(I) for different curcuit configurations of a transmission line. 
The telegraph's equation can be used to model a transmission line, and it can be written in terms of first and second derivatives
of V and I. The resulting system of PDEs can be readily solved for using numerical methods. 

For this project, finite differences were used to model the derivatives, and finite-difference time-domain (FDTD) methods were
applied to the numerical solution of the telegraph's equation. 4 circuit configurations were selected and modeled, and their MATLAB
programs can be found in the ./src folder:

1. oneresistor.m represents the simplest case possible: a circuit with one resistor. The resistor is added at the load end
of the transmission line. 

2. RCLload.m represents a circuit with parallel RCL elements added to the load. 

3. MatNetwork.m is a matched network. This means that the RC element at the load will perfectly match the network at the
simulated frequency, eliminating any reflections of V and I along the transmission line.

4. standing.m shows a standing wave pattern for both V and I. Lumped resistors are added at the load end so that the reflected 
wave interferes with incident wave.

The simulations contains frames showing the reflection patterns of voltage and current waves. The reflection index S11 
is computed and plotted.
