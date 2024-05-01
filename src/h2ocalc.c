#include <stdio.h>
#include <math.h>
#include <unistd.h> /* for getopt */
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "defaults.h"
#include "h2ocalc_2002.h"
//#include"h2ocalc_test.h"

/* gcc h2ocalc.c -o h2ocalc -lgsl -lgslcblas -lm -Wall  */

int nbprint(double x, double y[], int n);
int nbeprint(double x, double y[], int n);
int func (double t, const double y[],double f[], void *params);
int load_default_conc(const char *fname, double *y, int n);

int main(int argc, char **argv) {
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
    double doserate = DOSER;
    double simtime = SIMTIME;
    double rstart = RSTART;
    double resol = RESOL;
    double period = 0;
    double dose_per_pulse = 0;
    double tick = 0;
    // double ta = 1e-3; /* plotting clocks*/
    int pulses_left = 0;
    int pulse_counter = 0;
    signed long int tick_counter = 0;
    char *fname;
    int c;
    int flagp = 0, flagl = 0, flagc = 1; /* default: print each tick */

    /* parse options */
    opterr = 0;
    while ((c = getopt(argc, argv, "f:d:o:r:t:s:")) != -1) {
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
            sscanf(optarg, "%lf", &freq);
            break;
        case 'd':
            sscanf(optarg, "%lf", &doserate);
            break;
        case 'o':
            fname = optarg; /* not implemented */
            break;
        case 'r':
            sscanf(optarg, "%lf", &resol);
            break;
        case 't':
            sscanf(optarg, "%lf", &simtime);
            break;
        case 's':
            sscanf(optarg, "%lf", &rstart);
            break;
        case '?':
            printf("Options:\n");
            printf("  -p output per pulse\n");
            printf("  -c output per tick (default)\n");
            printf("  -l only end of output\n");
            printf("  -f frequency [Hz]\n");
            printf("  -d dose rate [Gy/min]\n");
            printf("  -o output file\n");
            printf("  -r tick resolution [sec]\n");
            printf("  -t simulation time [sec]\n");
            printf("  -s start irradiation at [sec] (default %.3e sec)\n", RSTART);
            printf("\n");
            exit(0);
            break;
        default:
            printf ("?? no handle ?? %i\n", c);
            exit(-1);
        }
    }

    period = 1 / freq;
    dose_per_pulse = (doserate / 60.0) / freq * EVJ; /* in eV per pulse per liter */


    if (period < resol) {
        printf(" *** Error: ");
        printf("Period must be larger than tick resolution.\n");
        exit(-1);
    }

    /* find tick size closest to the requested one */
    tick = period / (round(period / resol));

    /* number of pulses to be simulated */
    pulses_left = (int) ((simtime - rstart) * freq) + 1;

    /* print a header */
    printf("# Frequency    : %.3e Hz \n", freq);
    printf("# Period       : %.3e sec\n", period);
    printf("#\n");
    printf("# Tick size    : %.3e sec\n", tick);
    printf("# Time for sim : %.3e sec\n", simtime);
    printf("#\n");
    printf("# Average dose rate    : %.3e Gy/min\n", doserate );
    printf("# Pulse size   : %.3e eV/l/pulse\n", dose_per_pulse );
    printf("# Pulse count  : %i pulses\n", pulses_left);
    printf("#\n");
    printf("# Total delivered dose for this simulation\n");
    printf("#              : %.6e Gy \n", pulses_left * dose_per_pulse / EVJ);
    printf("#\n");


    /* if input file is given with starting conditions, then load these,
       otherwise use the model default values. */
    if (optind < argc) {
        load_default_conc(argv[optind], y, NSPECIES);
    }
    else {
        /* copy start conditions into working array */
        for (i = 0; i < NSPECIES; i++)
            y[i] = ystart[i];
    }

    /* Print first line with starting conditions. */
    nbeprint(t, y, NSPECIES);

    /* loop over all pulses */
    while (pulses_left != -1 && ti < simtime) {

        /* check if we are at a pulse time step */
        if (t >= ((pulse_counter * period) + rstart)) {
            if (flagp) /* print per pulse (pre pulse) */
                nbeprint(t, y, NSPECIES);

            printf("# Pulse! %i \n", pulses_left);
            for (j=0; j < NSPECIES; j++)
                y[j] += dose_per_pulse * gval[j] * 0.01 / NA;
            pulses_left--;
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
            nbeprint(t, y, NSPECIES);

    } /* end pulses left iterator */

    if (flagl) /* print last */
        nbeprint(t, y, NSPECIES);

    gsl_odeiv2_driver_free(d);
    return 0;
}


int nbprint(double x, double y[], int n) {
    int i;
    double sum = 0;
    printf("%6.4f   ", x);
    for (i = 0; i < n; i++) {
        printf("%6.2f ", y[i]);
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
        printf("%6.2e ", y[i]);
        sum += y[i];
    }
    printf("  %6.2e\n", sum);
    return 0;
}


int func (double t, const double y[],double f[], void *params) {
    // struct ODE_param *p = (struct ODE_params *) params;
    int i,j,k;
    double rate[NEQ];

    /* for each species ... */
    for (k = 0; k < NSPECIES; k++) {

        f[k] = 0;

        /* ...find all contributing and removing things */

        for (j = 0; j < NEQ; j++) {

            if (nmatrix[j][k] != 0) {

                rate[j] = nmatrix [j][k]; /* build the final equation */

                /* add generate I equations */
                rate[j] *= rconst[j]; /* get rate for reaction I_j */

                for (i = 0; i < NSPECIES; i++) {
                    if (nmatrix[j][i] == -1 ) {
                        rate[j] *= y[i]; /* first order kinetics */
                    }
                    if (nmatrix[j][i] == -2 ) {
                        rate[j] *= y[i]*y[i]; // second order kinetics
                    }
                }

                f[k] += rate[j];
            }
        }
    }

    return GSL_SUCCESS;
}


int load_default_conc(const char *fname, double *y, int n) {
    FILE *fp;
    char line[256]; // Buffer to store each line
    int i = 0;
    double val;

    fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", fname);
        exit(-1);
    }

    while (i < n && fgets(line, sizeof(line), fp) != NULL) {
        // Check if the line starts with a digit (simple filter for comments)
        if (sscanf(line, "%lf", &val) == 1) {
            y[i++] = val;
        }
    }

    fclose(fp);

    // Check if we have loaded enough values
    if (i != n) {
        fprintf(stderr,"Warning: expected %d values, but only %d were loaded from file %s.\n", n, i, fname);
        return 1; // Return error code if not all data could be loaded
    }

    return 0;
}