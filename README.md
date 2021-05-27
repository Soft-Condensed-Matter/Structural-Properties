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
