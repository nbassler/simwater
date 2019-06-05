#ifndef H2OCALC_TEST
#define H2OCALC_TEST

/*
   Simplified tests for ODE solver.
 */

#define NEQ 4 /* Number of equations */
#define NSPECIES 5 /* Number of species */
#define NSTEP 100
#define NSTART 0
#define NSTOP  10

float const ystart[NSPECIES] = {
    1,
    0,
    0,
    0,
    0
};


/* sochiometric matrix */

/*     A0 A1 A2 A3 .... */
/* v0                   */
/* v1                   */
/* v2                   */
/* v3                   */
/* ...                  */

int nmatrix[NEQ][NSPECIES] = {
    {-1, 1, 0, 0, 0},
    { 0,-2, 1, 0, 0},
    { 0, 0,-1, 2, 0},
    { 0, 0, 0,-1, 1}
};


/* define order of laws */
/* first or second order */
/* int rlaws[NEQ] = { */
/*   1, */
/*   2, */
/*   1, */
/*   1, */
/* }; */

double rconst[NEQ] = {
    1.0,
    3.0,
    1e-1,
    1.03,
};


#endif /* H2OCALC_TEST */
