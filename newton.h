#ifndef __newton
#define __newton

#include <complex.h>
#include <math.h>
#include <cerf.h>


/* Searches for roots to dielectric function
    2D iso dist with 2 vrms^2 = vperp^2, and warm in 
    dougnut axis with temp vpar = sqrt(2 k T / m)
    Fourier electric field at reciprocal k
    searches with res * res points
    in box [0, wmaxr] x [0, i wmaxi]
    starts at w and does one iteration
*/
int newton( double complex* w,
            double vperp, double vpar, 
            double k);


complex double Dxx(complex double * w, double vperp, double vpar, double k);

double complex pdisp(double complex zeta);
#endif
