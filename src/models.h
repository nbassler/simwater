struct model
{
    double *cstart; /* start condition array, nspecies long */
    double *gvals; /* gvals array, nspecies long */
    double *rconst; /* rconst array, neq long */
    int **nmatrix; /* nmatrix array, neq x nspecies long */
    int neq;
    int nspecies;
};
