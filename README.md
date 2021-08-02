Codes to compute structural properties from computer or experimental data

The codes are mainly developed in FORTRAN, contributions in other languages are welcome.


Code: Sk.f90

Description: Performs the 3D Fourier transform of the radial distribution
function contained in the file Gr.dat. The number of data in aforementioned 
file is especified by means of the volume fraction within the code


Code: PressGrv2.f90

Description: Computes the pressure for a square-well fluid using the radial
distribution function and the approaches of Henderson and Tago. This procedure 
requeries the contact value of g(r), and the values at the potential
discontinuities, this values are computed with an extrapolation 
of the g(r) present data employing a third degree Lagrange polynomial function


Code: PressTW.f90

Description: Computes the pressure for a triangular well fluid using the radial
distribution function separating the potential contributions on the hard-sphere 
and triangular contribs, the first one is computed with the Carnahan-Starling 
expresion and the second one explicitly computing the integral of potential 
derivative times the radial distribution function. This procedure requeries the 
contact value of g(r), which is computed using a third degree Lagrange polynom 
to extrapolate the value. Hence, provide the file with the g(r) data as well as
interaction range (lambda), density and temperatures is mandatory.


Code: nCoordination

Description: Determines the first coordination number in three dimensions using 
the radial distribution function. Program requieres temperature, number density 
and particle number at which the g(r) was determined, however, temperature it is 
not used at any computation stage, is used only to give context to results
