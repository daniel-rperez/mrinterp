# mrinterp
This is a set of functions written in Octave that implement a multigrid/multiresolution algorithm for 2-D interpolation of scalar data.

- mrinterp.m
This is the simplest multigrid/multiresolution interpolation algorithm. The function receives the X,Y,Z coordinates and the spatial resolution dxy, and returns a matrix with the interpolated values and a spatial domain.
- mrinterpfract.m
This is a multigrid/multiresolution interpolation algorithm with fractal extrapolation.
