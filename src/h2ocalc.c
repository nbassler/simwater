#include <stdio.h>
#include <math.h>
#include <unistd.h> /* for getopt */
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "h2ocalc.h"
//#include"h2ocalc_test.h"

/* gcc h2ocalc.c -o h2ocalc -lgsl -lgslcblas -lm -Wall  */

/* struct ODE_param{ */
/*   int a[NEQ];  /\* dont worry let just run with global vars for now *\/ */
/*   float b[NEQ]; */
/* }; */

int nbprint(double x, double y[], int n) {
    int i;
    double sum = 0;
    printf("%6.4f   ",x);
    for (i = 0; i < n; i++) {
        printf("%6.2f ",y[i]);
        sum += y[i];
    }
    printf("  %6.2f\n",sum);
    return 0;
}


int nbeprint(double x, double y[], int n) {
    int i;
    double sum = 0;

    printf("%6.4e   ",x);
    for (i = 0; i < n; i++) {
        printf("%6.2e ",y[i]);
        sum += y[i];
    }
    printf("  %6.2e\n",sum);
    return 0;
}


int func (double t, const double y[],double f[], void *params) {
    // struct ODE_param *p = (struct ODE_params *) params;
    int i,j,k;
    double rate[NEQ];

    //printf("func\n");

    /* for each species ... */
    for (k = 0; k < NSPECIES; k++) {
        //printf("SPECIES %i\n",k);

        f[k] = 0;

        /* ...find all contributing and removing things */

        for (j = 0; j < NEQ; j++) {

            if (nmatrix[j][k] != 0) {

                rate[j] = nmatrix [j][k]; /* build the final equation */
                //printf("%3i ",nmatrix[j][k]);

                /* add generate I equations */
                rate[j] *= rconst[j]; /* get rate for reaction I_j */

                //printf(" * k%i ",j);
                for (i = 0; i < NSPECIES; i++) {
                    if (nmatrix[j][i] == -1 ) {
                        rate[j] *= y[i]; /* first order kinetics */
                        //printf("* y[%i] ",i);
                    }
                    if (nmatrix[j][i] == -2 ) {
                        rate[j] *= y[i]*y[i]; // second order kinetics
                        //printf("* y[%i]^2 ",i);
                    }
                }

                f[k] += rate[j];
                //printf(" + \n");
            }
        }


        //printf(" ----  k f[k]: %i     %e ---- \n",k,f[k]);

        //    nbeprint(t,y,NSPECIES);
        //    if (k ==8)
        //      exit(0);
    }

    //nbeprint(t,y,NSPECIES);
    //nbeprint(t,f,NSPECIES);


    /* add function for radiation, assumed in last bin */
    //f[NEQ-1] = 0;

    return GSL_SUCCESS;
}


int main(int argc, char **argv) {
    //  struct ODE_param *p = (struct ODE_params *) params;

    double dummy_param = 0;
    gsl_odeiv2_system sys = {func, NULL, NSPECIES, &dummy_param};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
                                                         gsl_odeiv2_step_rk8pd,
                                                         1e-6,
                                                         1e-6,
                                                         0.0);
    int i,j,status;
    double t = NSTART;
    //  double t1 = NSTOP;
    //double foobar; /* trash parameter */
    double ti;
    double y[NSPECIES];

    double freq = FREQ;
    double doser = DOSER;
    double simtime = SIMTIME;
    double resol = RESOL;
    double period = 0;
    double dpulse = 0;
    double tick = 0;
    double ta = 1e-3; /* plotting clocks*/
    int pulse_left = 0;
    int pulse_counter = 0;
    signed long int tick_counter = 0;
    char *fname;
    int c;
    int flagp = 0, flagl = 0, flagc = 1; /* default: print each tick */

    /* parse options */
    //  strcpy(fname, "foobar.dat");
    opterr = 0;
    while ((c = getopt(argc, argv, "f:d:o:r:t:")) != -1) {
        switch(c) {
        case 'l':
            flagl = 1;
            flagp = 0;
            flagc = 0;
            break;
        case 'p':
            flagl = 0;
            flagp = 1;
            flagc = 0;
            break;
        case 'c':
            flagl = 0;
            flagp = 0;
            flagc = 1;
            break;
        case 'f':
            sscanf(optarg,"%lf",&freq);
            break;
        case 'd':
            sscanf(optarg,"%lf",&doser);
            break;
        case 'o':
            fname = optarg; /* not implemented */
            break;
        case 'r':
            sscanf(optarg,"%lf",&resol);
            break;
        case 't':
            sscanf(optarg,"%lf",&simtime);
            break;
        case '?':
            printf("Options:\n");
            printf("  -p output per pulse\n");
            printf("  -c output per tick\n");
            printf("  -l only end of output\n");
            printf("  -f frequency [Hz]\n");
            printf("  -d dose rate [Gy/min]\n");
            printf("  -o output file\n");
            printf("  -r tick resolution [sec]\n");
            printf("  -t simulation time [sec]\n");
            printf("\n");
            exit(0);
            break;
        default:
            printf ("?? no handle ?? %i\n", c);
            exit(-1);
        }
    }


    period = 1 / freq;
    dpulse = (doser / 60.0) / freq * EVJ; // in eV per pulse per liter


    if (period < resol) {
        printf(" *** Error: ");
        printf("Period must be larger than tick resolution.\n");
        exit(-1);
    }

    /* find ticksize closest to the requested one */
    tick = period/(round(period/ resol));

    /* number of pulses to be simulated */
    pulse_left = (int) ((simtime-RSTART) * freq);

    /* print a header */
    printf("# Frequency    : %.3e Hz \n",freq);
    printf("# Period       : %.3e sec\n",period);
    printf("#\n");
    printf("# Tick size    : %.3e sec\n",tick);
    printf("# Time for sim : %.3e sec\n",simtime);
    printf("#\n");
    printf("# Dose rate    : %.3e Gy/min\n",doser );
    printf("# Pulse size   : %.3e eV/l/pulse\n",dpulse );
    printf("# Pulse count  : %i pulses\n", pulse_left);
    printf("#\n");
    printf("# Total delivered dose for this simulation\n");
    printf("#              : %.6e Gy \n", pulse_left * dpulse / EVJ);
    printf("#\n");

    /* copy start conditions into working array */
    for (i = 0; i < NSPECIES; i++)
        y[i] = ystart[i];

    /* Print first line with starting conditions. */
    nbeprint(t,y,NSPECIES);

    /* loop over all pulses */
    while (pulse_left != 0) {

        /* check if we are at a pulse time step */
        if (t >= ((pulse_counter * period) + RSTART)) {
            if (flagp) /* print per pulse (pre pulse) */
                nbeprint(t,y,NSPECIES);

            printf("# Pulse! %i \n", pulse_left);
            for (j=0; j < NSPECIES; j++)
                y[j] += dpulse * gval[j] * 0.01 / NA;
            pulse_left--;
            pulse_counter++;
        }

        tick_counter++;
        ti = tick_counter * tick;

        status = gsl_odeiv2_driver_apply(d, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf("error, return value = %d\n", status);
            break;
        }

        if (flagc) /* print per tick */
            nbeprint(t,y,NSPECIES);

    } /* end pulses left iterator */

    if (flagl) /* print last */
        nbeprint(t,y,NSPECIES);

    gsl_odeiv2_driver_free(d);
    return 0;
}

