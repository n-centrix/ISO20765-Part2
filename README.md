# ISO20765-Part2
ISO20765-Part 2 Equation of State (GERG2008) coded in Fortran g95

The GERG-2008 equation of state was developed by the University of Bochum in Germany as a new wide-range equation of state for the volumetric and caloric properties of natural gases and other mixtures. It is the subject of the ISO 20765 Part 2.

The ranges of temperature, pressure, and composition to which the GERG-2008 equation of state applies are much wider than the ISO 20765 Part 1 (AGA-8) equation and cover an extended range of application.

ISO-20765 methods use a standardized 21-component gas system in which all of the major and minor components of natural gas are included. Trace component present but not identified as one of the 21 specified components may be reassigned (see ISO-20765 for details).

The SRK equation of state is implemented as part of the ISO 20765 routines (density solver) and is applied in order to obtain an initial approximation of the density roots (vapor phase). More specifically, SRK equation of state is used to narrow down the root search interval of the density solver for ISO20765 part II.
For reference, SRK equation of state formulation, i.e. alpha function, mixing rules and binary interaction coefficients are those from the API Technical Data Book 6th edition.

The present implementation applies to gas phase (single phase). 
It is the user responsibility to verify the gas condition vs. envelope boundaries using an appropriate process simulator.
