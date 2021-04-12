#include "newton.h"
#include <stdio.h>


//One iteration of newton
int newton( double complex *w,
            double vperp, double vpar, 
            double k){
    double kcwp = k * k;
    kcwp = 1. + kcwp;
    double complex w2 = *w * *w;
    double tmp = vpar * k;
    double complex zeta = *w / tmp;
    double complex zeta2 = zeta * zeta;


    double A = vperp / vpar;
    A = A * A;

    double complex Z = pdisp(zeta);

    if(isfinite(cabs(Z))){
        //double complex num = w2 + kcwp - A * (2. * zeta + 2. * (zeta2 - zeta) * Z - 1.);
        //double complex den = 2. * w2 - A * zeta * (2. + (2 * zeta - 1) * Z);
        //*w *= num / den;
        double complex num = w2 - kcwp + A * (1. + zeta * Z);
        double complex den = 2. * *w - A / tmp * (2. * zeta + (2. * zeta2 - 1.) * Z);
        *w -= num / den;
        //double complex den = 2. * w2 - A * zeta * (2. * zeta + (2. * zeta2 - 1.) * Z);
        //*w *= 1. - num / den;
        if(isfinite(cabs(*w))){
            return 0;
        }else{
            return -1;
        }
    }
    return -1;
}


complex double Dxx(complex double * w, double vperp, double vpar, double k){
    double A = vperp / vpar;
    A = A * A;
    complex double w2 = *w * *w;
    double kcwp = 1. + k*k;


    double tmp = vpar * k;
    complex double zeta = *w / tmp;
    complex double Z = pdisp(zeta); 

    return w2 - kcwp + A * (1. + zeta * Z);
}


/**
   * calculates the plasma dispersion function a double complex zeta
   * It is possible it will be implemented to complex.h and thus the name cerf is reserved 
**/
double complex pdisp(double complex zeta){
    double complex tmp = _Complex_I * (M_PI_2 * M_2_SQRTPI); //i sqrt(pi)

    //Expand at infinity to remove inf * 0
    // as expression is with 1 + cerfi(zeta), 
    // and cerfi is generally better accuracy than 1E-13, then with 1+cerfi close to 1
    // it is correct generally correct to to 1e-3, i.e. we should expand when exp(-zeta^2) > 1e3
    // which is about zeta > 2.627 
    //if (creal(zeta * zeta) < 7.){
    //    double complex invzeta = 1. / zeta;
    //    double complex invzeta2 = invzeta * invzeta;
    //    double complex invzeta3 = invzeta2 * invzeta;
    //    double complex invzeta5 = invzeta3 * invzeta2;
    //    double complex invzeta7 = invzeta5 * invzeta2;
    //    double complex exp2 = cexp(- zeta * zeta);
    //    exp2 = 2. * tmp * exp2;
    //    return - invzeta - 0.5 * invzeta3 - 0.75 * invzeta5 - 1.875 * invzeta7 + exp2;
    //}

    //no expansion default case
    double complex tmp2 = cexp(- zeta * zeta);
    double complex tmp3 = cerfi(zeta);
    tmp3 = 1. + _Complex_I * tmp3;
    return tmp * tmp2 * tmp3;
}
