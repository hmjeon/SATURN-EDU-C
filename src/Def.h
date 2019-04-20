//
// Problem description for plane stress condition
//

//#define     _FILEINPUT
#define     DOMAINSIZE      20

#define     YOUNG           1.7472e7            // Young's modulus
#define     POSSION         0.3                 // Poisson's ratio
#define     BXFORCE         0.0                 // x-direction body force
#define     BYFORCE        -1.0                 // y-direction body force
#define     THICK           1.0                 // thickness
#define     X_WIDTH         1.0                 // x length
#define     Y_WIDTH         1.0                 // y length
#define     N_I_NODE        (DOMAINSIZE + 1)    // the number of nodes in the x-direction
#define     N_J_NODE        (DOMAINSIZE + 1)    // the number of nodes in the y-direction

#define     nDIM            2                   // problem dimension
#define     nNPE            4                   // # of nodes per element
#define     nDPN            2                   // # of DOFs per node(u, v)