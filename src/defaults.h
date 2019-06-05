#ifndef DEFAULTS
#define DEFAULTS

/*

   Typical cellular values:
   pH  = 7.4, intra cellular
   pO2 = ? TODO
 */

//#define NSTEP    100  /* Number of steps between pulses */
//#define NPULSE   100   /* Number of pulses */
#define NSTART   0
#define SIMTIME  1.0e-2 /* time to be simulated in sec */
//#define NSTOP    0.00002
#define FREQ     1e6   /* 1 MHz blips */
#define RSTART   5e-6  /* start after 5 musec */
#define DOSER    1.00  /* dose rate in Gy/min */
#define RESOL    2e-8  /* 20 nsec resolution */

#define NA       6.02214129e23      /* Avogadro constant */
#define EVJ      6.24150934e18      /* eV per joule */

#endif /* DEFAULTS */
