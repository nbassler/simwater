#include <stdio.h>
#include <math.h>

#include "defaults.h"
#include "h2ocalc.h"


double nbpulse(double tau, double T, double t) {

    int n;
    double temp;

    temp = tau / T;
    for (n = 1; n < 13; n++)
        temp +=  2.0 / (n * M_PI) * sin(M_PI * n * tau / T) * cos(2.0 * M_PI * n * t / T);
    return temp;
}


int main(int argc, char **argv) {

    double time = 0;
    double dt = 0.001;
    //int i;

    //double phase = 0;
    double a = 0;

    double ff = 3; // 10 Hz
    double duty = 0.5;

    double tau, period;

    period = 1 / ff;
    tau = duty * period;

    //  for (i = 0; i < 1e6; i++) {
    while (time < 1.0) {

        //    phase = modf(time,&a);
        a = nbpulse(tau, period, time);
        printf ("%.3f %.3f\n", time, a);

        time += dt;

    }
    return 0;
}

