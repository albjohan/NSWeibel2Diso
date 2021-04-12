#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>
#include <time.h>


#include "newton.h"


int main(int argc, char * argv[]){
    double kmax;
    double vperp;
    double vpar;
    double wmaxr;
    double wmaxi;
    int res;
    int nk = 1;
    int fnum;

    int ierr, num_procs, id;
    MPI_Status status;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

    time_t seconds;
    if(id == 0){
        time(&seconds);
    }
    
    //PARSE INPUT
    {
        int opt;
        while ( (opt = getopt(argc,argv,"k:y:z:R:I:r:s:")) != -1){
            switch(opt){
                case 'k':
                    kmax = atof(optarg);
                    break;
                case 'y':
                    vperp = atof(optarg);
                    break;
                case 'z':
                    vpar = atof(optarg);
                    break;
                case 'R':
                    wmaxr = atof(optarg);
                    break;
                case 'I':
                    wmaxi = atof(optarg);
                    break;
                case 'r':
                    res = atoi(optarg);
                    break;
                case 's':
                    nk = atoi(optarg);
                    break;
            }
        }
    }

    /********* Newton iterations ***********/
    double tol = 1e-12;
    double dwr = wmaxr / (double)(res - 1);
    double dwi = wmaxi / (double)(res - 1);
    double dk = kmax / nk;
    double complex * w = malloc(sizeof(double complex) * (res/num_procs + 1) * res);
    int maxiter = 10000;


    for(int ik = 0; ik < nk; ++ik){
        memset(w, 0x00, sizeof(double complex) * (res/num_procs + 1) * res/sizeof(char));
        double k = dk * (ik + 1);
        if(id == 0){
            printf("Starting %ith iter of k out of %i\n", ik, nk);
        }

        for(int ir = 0; ir * num_procs + id < res; ++ir){
            double wr = dwr * (ir * num_procs + id);
            for(int ii = 0; ii < res; ++ii){
                double wi = dwi * ii;
                w[ir * res + ii] = wr + _Complex_I * wi;
                int iter = 0;
                while(1){
                    if( -1 == newton(w + ir * res + ii, vperp, vpar, k) ){
                        w[ir * res + ii] = 0. + 0. * _Complex_I;
                        break;
                    }
                    double tmp = cabs(Dxx(w + ir * res + ii, vperp, vpar, k));
                    if (isfinite(tmp)){
                        if (tmp < tol){
                            break;
                        }
                    }else{
                        w[ir * res + ii] = 0. + 0. * _Complex_I;
                        break;
                    }
                    if(cabs(w[ir * res + ii]) > 1e2){
                        w[ir * res + ii] = 0. + 0. * _Complex_I;
                        break;
                    }
                    if(maxiter < ++iter){
                        w[ir * res + ii] = 0. + 0. * _Complex_I;
                        break;
                    }
                }
            }
        }
        

        /************* write data to file *************/
        if(id == 0){//<----------- MASTER RECEIVE AND WRITE TO FILE ---------->
            char fname[] = "newton_convergence_000.dat";
            fname[19] = (ik % 1000) / 100 + '0';
            fname[20] = (ik % 100) / 10 + '0';
            fname[21] = ik % 10 + '0';
            FILE* fconv = fopen(fname, "w");
            for(int ir = 0; ir < res; ++ir){
                if(ir % num_procs == 0){
                    fwrite(w + (ir / num_procs) * res, sizeof(double complex), res, fconv);
                }else{
                    MPI_Recv(w, 2 * res, MPI_DOUBLE, ir % num_procs, 0, MPI_COMM_WORLD, &status);
                    fwrite(w, sizeof(double complex), res, fconv);
                }
            }
            fclose(fconv);

        }else{//<------------------ SLAVE SEND DATA ----------------->
            for(int ir = 0; ir * num_procs + id < res; ++ir){
                MPI_Send(w + ir * res, 2 * res, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }

    }//<------------ END k LOOP ------------------>
            
    if(id == 0){
        double * karr = malloc(sizeof(double) * nk);
        for(int ik = 0; ik < nk; ++ik){
            karr[ik] = dk * (ik + 1);
        }
        FILE* filek = fopen("k.dat", "w");
        fwrite(karr, sizeof(double), nk, filek);
        fclose(filek);
        memset(karr, 0x00, (sizeof(double) * nk)/sizeof(char));
        free(karr);
        printf("Runtime: %li\n", time(NULL) - seconds);
    }

    memset(w, 0x00, (sizeof(double complex) * (res/num_procs + 1) * res)/sizeof(char));
    free(w);
    MPI_Finalize();
    return 0;
}
