#include <stdio.h>
#include <stdlib.h>
#include <math.h>       /* pow */
#include <assert.h>

/*#define zero 0.0
#define one 1.0
#define two 2.0
#define three 3.0
#define four 4.0
#define six 6.0
#define half 0.5
#define fourth 0.25
#define fifth 0.2
#define sixth 1.0/6.0
#define big 1.7976931348623157E+308
#define small 2.2250738585072014E-308*/


int x_nodes, y_nodes;
double xmin, xmax, ymin, ymax, dx, dy, length;

double re, rho, u_lid, p_guage, mu, nu, dtd;
// character (len=8)
char solver[8];

int max_iter, resid_out, write_solution;
//real(dp) ::
double cfl, k, c2, Cx, Cy, rkappa, conv_toler, visc_eps;

double zero;// = 0.0;
double one;// = 1.0;
double two;// = 2.0;
double three;// = 3.0;
double four;// = 4.0;
double six;// = 6.0;
double half;// = 0.5;
double fourth;// = 0.25;
double fifth;// = 0.2;
double sixth;// = 1.0/6.0;
double big;// = 1.7976931348623157E+308;
double small;// = 2.2250738585072014E-308;

double* dt;
double* beta;
double* soln; //real(c_double), allocatable, target :: soln(:,:,:)
double* soln_new; // real(c_double), allocatable, target :: soln_new(:,:,:)
