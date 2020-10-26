//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>       /* pow */
//#include <assert.h>
#include "global_vars.h"

/*real(dp), parameter :: zero   = 0.0_dp
real(dp), parameter :: sixth  = 1.0_dp/6.0_dp
real(dp), parameter :: fifth  = 0.2_dp
real(dp), parameter :: fourth = 0.25_dp
real(dp), parameter :: third  = 1.0_dp/3.0_dp
real(dp), parameter :: half   = 0.5_dp
real(dp), parameter :: one    = 1.0_dp
real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: three  = 3.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: six    = 6.0_dp
  real(dp), parameter :: big    = huge(zero)
  real(dp), parameter :: small  = tiny(zero)
  

double zero = 0.0;
double one = 1.0;
double two = 2.0;
double three = 3.0;
double four = 4.0;
double six = 6.0;
double half = 0.5;
double fourth = 0.25;
double fifth = 0.2;
double sixth = 1.0/6.0;
double big = 1.7976931348623157E+308;
double small = 2.2250738585072014E-308;
*/
//double Cx = 0.01;



/*int x_nodes, y_nodes;
double xmin, xmax, ymin, ymax, dx, dy, length;

double re, rho, u_lid, p_guage, mu, nu, dtd;
// character (len=8)
char solver[8];

int max_iter, resid_out, write_solution;
//real(dp) ::
//double cfl, k, c2, Cx, Cy, rkappa, conv_toler, visc_eps;

//real(dp), allocatable, dimension(:,:)   ::
//double** dt;
//double** beta;
double* dt;
double* beta;
//real(dp), allocatable, dimension(:,:,:) ::
//double*** soln;
//double*** soln_new;
//double* soln;
//double* soln_new;
*/

void ldc_allocate(double* soln, double* soln_new, double* dt, double* beta){

//    printf("[C_Setup_File] ## Allocating Memory for soln and soln_new...\n");

    //use set_precision,  only : dp
    //use set_constants,  only : zero, one, four, big

    //real(dp)               :: x, y, xmaxr, ymaxr, xminr, yminr
    double x, y, xmaxr, ymaxr, xminr, yminr;
    //real(dp), dimension(2) :: r
    double r[2];
    //integer                :: i, j
    int i, j, k;

    //test_dbg
//    double* test_dbg;
//    int test_cnt = 1;

    //continue
    //--------------------------------------------------------------------------//
    //------------------ allocating Three dimensional Arrays -------------------//
    //soln = (double*) malloc(sizeof(double)*3*y_nodes*x_nodes);
    //soln_new = (double*) malloc(sizeof(double)*3*y_nodes*x_nodes);
//    printf("[C_Setup_File] ## Trace 1\n");
    //assert(soln != NULL);
    //assert(soln_new != NULL);
    /*****************************************************************************/
    //--------------------------------------------------------------------------//
    //------------------- allocating Two dimensional Arrays --------------------//
    //dt = (double*) malloc(sizeof(double)*x_nodes*y_nodes);
    //beta = (double*) malloc(sizeof(double)*x_nodes*y_nodes);
    //assert(dt != NULL);
    //assert(beta != NULL);
    /*****************************************************************************/
	
    //allocate(soln(x_nodes,y_nodes,3), soln_new(x_nodes,y_nodes,3))
    //allocate(dt(x_nodes,y_nodes), beta(x_nodes,y_nodes))

// Set BC's and initialize matrices

zero = 0.0;
one = 1.0;
two = 2.0;
three = 3.0;
four = 4.0;
six = 6.0;
half = 0.5;
fourth = 0.25;
fifth = 0.2;
sixth = 1.0/6.0;
big = 1.7976931348623157E+308;
small = 2.2250738585072014E-308;



// Test_dbg
/*test_dbg = (double*) malloc(2*4*10*sizeof(double));
for (j=0; j<4; j=j+1){
   for (i=0; i<10; i=i+1){
      test_dbg[(j*10)+i] = test_cnt;
      test_dbg[(4*10) + (j*10) + i] = test_cnt + 1;
      test_cnt = test_cnt + 1;
   }
}
for (j=0; j<4; j=j+1){
   for (i=0; i<10; i=i+1){
      printf("test_dbg[%d] = %lf\n", (j*10)+i, test_dbg[(j*10)+i]);
   }
}
for (j=0; j<4; j=j+1){
   for (i=0; i<10; i=i+1){
      printf("test_dbg[%d] = %lf\n", ((4*10)+(j*10)+i), test_dbg[(4*10)+(j*10)+i]);
   }
}*/

/*for (i=0; i<80; i=i+1){
   printf("test_dbg[%d] = %lf\n", i, test_dbg[i]);
}*/



	for (j=0; j<y_nodes; j=j+1){
		for (i=0; i<x_nodes; i=i+1){
		   xmaxr = one;
        	   ymaxr = one;
            	   xminr = zero;
            	   yminr = zero;
			
		   x = (xmaxr - xminr)*(i+1)/(x_nodes - 1);
            	   y = (ymaxr - yminr)*(j+1)/(y_nodes - 1);
			
		   r[0] = x - ( (xmaxr - xminr)*0.5 + xminr );
            	   r[1] = y - ( (ymaxr - yminr)*0.5 + yminr );
				
            	   //soln[0][j][i] = p_guage;
                   soln[(j*x_nodes) + i] = p_guage;
            	   //soln[1][j][i] = r[1];
                   soln[(y_nodes*x_nodes) + (j*x_nodes) + i] = r[1];
              	   //soln[2][j][i] = -r[0];
                   soln[(2*x_nodes*y_nodes) + (j*x_nodes) + i] = -r[0];
                   //printf("soln(%d,%d,%d) = %1.20f\n", i+1, j+1, 1, soln[(j*x_nodes) + i]);
                   //printf("soln(%d,%d,%d) = %1.20f\n", i+1, j+1, 2, soln[(y_nodes*x_nodes) + (j*x_nodes) + i]);
                   //printf("soln(%d,%d,%d) = %1.20f\n", i+1, j+1, 3, soln[(2*x_nodes*y_nodes) + (j*x_nodes) + i]); 
		}
	}
//        printf("[C_Setup_File] ## Trace 8\n");
	
	for (i=0; i<x_nodes; i=i+1){
		//soln[1][0][i] = zero;  =>  soln[(1*x_nodes*y_nodes) + (0*x_nodes) + i]
                soln[(x_nodes*y_nodes) + i] = zero; 

	}
//       printf("[C_Setup_File] ## Trace 9\n");
	for (i=0; i<y_nodes; i=i+1){
		//soln[1][i][0] = zero;
                soln[(x_nodes*y_nodes) + (i*x_nodes)] = zero;
                soln[(x_nodes*y_nodes) + (i*x_nodes) + (x_nodes-1)] = zero;

	}

        for (i=0; i<x_nodes; i=i+1){
           //soln[2][0][i] = zero;  ==> soln[(2*x_nodes*y_nodes) + (0*x_nodes) + i]
           soln[(2*x_nodes*y_nodes) + i] = zero;
        }
        for (i=0; i<y_nodes; i=i+1){
           //soln[2][i][0] = zero;
           soln[(2*x_nodes*y_nodes) + (i*x_nodes)] = zero;
           //soln[2][i][x_nodes-1] = zero;
           soln[(2*x_nodes*y_nodes) + (i*x_nodes) + (x_nodes-1)] = zero;

        }
        for (i=0; i<x_nodes; i=i+1){
           //soln[2][y_nodes-1][i] = zero;
           soln[(2*x_nodes*y_nodes) + ((y_nodes-1)*x_nodes) + i] = zero;
           //soln[1][y_nodes-1][i] = u_lid;
           soln[(x_nodes*y_nodes) + ((y_nodes-1)*x_nodes) + i] = u_lid;
        }

        //printf("[C_Setup_File] ## Trace 10\n");
	//soln_new = soln
	for (i = 0; i < (3*y_nodes*x_nodes); i=i+1)
    	{
            soln_new[i] = soln[i];
        }
        //printf("[C_Setup_File] ## Trace 11\n");
	
    
    // Look up head at the triple ***beta/dt Initialize***
    //dt    = big;
    //beta  = big;

// Set constant variables

    mu = rho*u_lid*length/re;  //Dynamic viscosity, Pa*s
    nu = mu/rho;               //Kinematic viscosity, m^2/s
    dx = (xmax-xmin)/(x_nodes-1);
    dy = (ymax-ymin)/(y_nodes-1);
    //dtd = dx**2/(four*nu)
    dtd = pow(dx,2.0)/(four*nu);
    //printf("soln(1,1,3) = %1.10f, soln(2,1,3) = %1.10f\n", soln[(2*x_nodes*y_nodes) + (0*x_nodes) + 0], soln[(2*x_nodes*y_nodes) + (0*x_nodes) + 1]);
}
//  end subroutine ldc_allocate

void ldc_deallocate(){
	//deallocate(dt, beta, soln, soln_new)

//    free(soln);	// Old Body
//    free(soln_new);	// Old Body
//    free(dt);		// Old Body
//    free(beta);	// Old Body

/*	int i,j;

    for(i=0;i<3;i++)
    {
    	for(j=0;j<y_nodes;j++)
    	{
    		//free(soln[i][j]);
		free(soln[i][j]);
		free(soln_new[i][j]);
    	}
    	//free(soln[i]);
	free(soln[i]);
	free(soln_new[j]);
    }
    //free(soln);
    free(soln);
    free(soln_new);
	
    // Freeing dt and beta
    for(i=0;i<y_nodes;i++)
    {
    	free(dt[i]);
	free(beta[i]);
    }
    free(dt);
    free(beta);*/
}

