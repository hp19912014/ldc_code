#include "global_vars.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })




void ldc_explicit_iter(int x_nodes, int y_nodes, double dx, double dy, const double* __restrict__ soln, double* __restrict__ soln_new, double L2[3]){
   double R1, R2, R3;
   double dt, uvel2, beta2;
   double dpdx, dudx, dvdx, dpdy, dudy, dvdy, d2udx2, d2vdx2, d2udy2, d2vdy2, d4pdx4, d4pdy4; 
   double lambda_x, lambda_y, lambda_max, dtconv, dtvisc, viscx, viscy;
   double vel2ref,rhoinv, rmu;

   double r_twodx, r_twody, r_dxdx, r_dydy, r_dxdxdxdx, r_dydydydy;

   double L2_1, L2_2, L2_3;

   int i, j, array_size, array_init, dbg, y, x;
   double Pweightfactor, internal_nodes;

   L2_1   = zero;
   L2_2   = zero;
   L2_3   = zero;

   //reciprocal factors for finite differences
   r_twodx    = one / (two*dx);
   r_twody    = one / (two*dy);
   r_dxdx     = one / (dx*dx);

   r_dydy     = one / (dy*dy);
   r_dxdxdxdx = one / pow(dx, 4);
   r_dydydydy = one / pow(dy, 4);

   //Nodes in interior of domain (not including boundary nodes)
   internal_nodes = (x_nodes-2)*(y_nodes-2); 

   vel2ref = pow(u_lid, 2);
   rhoinv = one / rho;
   rmu = rho*u_lid*length/re;          //Viscosity (N*s/m^2)

   //Diffusive/viscous limit - constant over domain
   dtvisc = fourth*(dx*dy)*(rho/rmu);  

   array_init = 0;
   array_size = 3*x_nodes*y_nodes - 1;
//printf("Iteration:\n");
//printf("[Trace MSG] ## 1\n");

//   printf("##Initial Data##:: x_nodes= %d, y_nodes= %d, dx= %lf, dy= %lf, zero= %lf, two= %lf, four= %lf, six= %lf, fourth= %lf, half= %lf, one= %lf, dtd= %lf, cfl= %lf, k= %lf, u_lid= %lf, p_guage= %lf, re= %lf, rkappa= %lf, conv_toler= %lf, resid_out= %d, rho= %lf, nu= %lf, visc_eps= %lf, c2= %lf, length= %lf, Cx= %lf, Cy= %lf\n", x_nodes, y_nodes, dx, dy, zero, two, four, six, fourth, half, one, dtd, cfl, k, u_lid, p_guage, re, rkappa, conv_toler, resid_out, rho, nu, visc_eps, c2, length, Cx, Cy);
   
   #ifdef _OPENACC
      #if defined(ITHREADS) && defined(JTHREADS)
         //!$acc kernels present(soln(:,:,:), soln_new(:,:,:)) 
         #pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size])  
      #else
         //!$acc kernels loop independent present(soln(:,:,:), soln_new(:,:,:))
         #pragma acc kernels loop independent present(soln[array_init:array_size], soln_new[array_init:array_size])
      #endif 
   #elif defined _OPENMP
      //!$omp parallel do private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
  //    #pragma omp parallel for private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
   #endif
   
   #if defined _OPENACC && defined(ITHREADS) && defined(JTHREADS)
      //!$acc loop independent   vector(JTHREADS)
      #pragma acc loop independent  vector(JTHREADS) 
   #endif    
    
   //do j = 3, y_nodes-2
   for (i = 2; i < x_nodes - 2; i = i+1){

   #if defined _OPENACC && defined(ITHREADS) && defined(JTHREADS)
      //!$acc loop independent vector(ITHREADS)
      #pragma acc loop independent vector(ITHREADS)
   #endif
      for(j = 2; j < y_nodes-2; j = j+1){

         //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		 
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)

         //d4pdx4 = ( soln(i+2,j,1) - four*soln(i+1,j,1) + six*soln(i,j,1) - four*soln(i-1,j,1) + soln(i-2,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+2)] - four*soln[(j*x_nodes)+(i+1)] + six*soln[(j*x_nodes)+i] - four*soln[(j*x_nodes)+(i-1)] + soln[(j*x_nodes)+(i-2)]) * r_dxdxdxdx;
         //d4pdy4 = ( soln(i,j+2,1) - four*soln(i,j+1,1) + six*soln(i,j,1) - four*soln(i,j-1,1) + soln(i,j-2,1) )*r_dydydydy
         d4pdy4 = (soln[((j+2)*x_nodes)+ i] - four*soln[((j+1)*x_nodes)+i] + six*soln[(j*x_nodes)+i] - four*soln[((j-1)*x_nodes)+i] + soln[((j-2)*x_nodes)+i]) * r_dydydydy;

         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;

         // Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 

         R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

      }		// For on i
   }		// For on j
   
#ifdef _OPENACC
//!$acc end kernels
#elif defined _OPENMP
//!$omp end parallel do
#endif

//cudaDeviceSynchronize();


   L2[0] = L2_1;
   L2[1] = L2_2;
   L2[2] = L2_3;
   L2_1  = zero;
   L2_2  = zero;
   L2_3  = zero;

#ifdef _OPENACC
//!$acc kernels present(soln(:,:,:), soln_new(:,:,:)) 
#pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size]) 
//!$acc loop independent reduction(+:L2_1,L2_2,L2_3)
#pragma acc loop independent reduction(+:L2_1,L2_2,L2_3)
#elif defined _OPENMP
//!$omp parallel do private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
//#pragma omp parallel for private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
#endif

   for(j=2; j < y_nodes-2; j=j+1){
      
	  i = 1;
	     //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
  
         //uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		 
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)

         //d4pdx4 = ( soln(i+3,j,1) - four*soln(i+2,j,1) + six*soln(i+1,j,1) - four*soln(i,j,1) + soln(i-1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+3)] - four*soln[(j*x_nodes)+(i+2)] + six*soln[(j*x_nodes)+(i+1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i-1)]) * r_dxdxdxdx;
         //d4pdy4 = ( soln(i,j+2,1) - four*soln(i,j+1,1) + six*soln(i,j,1) - four*soln(i,j-1,1) + soln(i,j-2,1) )*r_dydydydy
         d4pdy4 = (soln[((j+2)*x_nodes)+ i] - four*soln[((j+1)*x_nodes)+i] + six*soln[(j*x_nodes)+i] - four*soln[((j-1)*x_nodes)+i] + soln[((j-2)*x_nodes)+i]) * r_dydydydy;
         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;

         // Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 

         R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);  
  
                        
         //Extrapolate pressure
         i = 0;
         //soln_new(i,j,1) = two*soln_new(i+1,j,1) - soln_new(i+2,j,1)
	 soln_new[(j*x_nodes)+i] = two*soln_new[(j*x_nodes)+(i+1)] - soln_new[(j*x_nodes)+(i+2)];

         i = x_nodes-2;
         //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
  
         //uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		 
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
		 
	 //d4pdx4 = ( soln(i-3,j,1) - four*soln(i-2,j,1) + six*soln(i-1,j,1) - four*soln(i,j,1) + soln(i+1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i-3)] - four*soln[(j*x_nodes)+(i-2)] + six*soln[(j*x_nodes)+(i-1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i+1)]) * r_dxdxdxdx;
         //d4pdy4 = ( soln(i,j+2,1) - four*soln(i,j+1,1) + six*soln(i,j,1) - four*soln(i,j-1,1) + soln(i,j-2,1) )*r_dydydydy
         d4pdy4 = (soln[((j+2)*x_nodes)+ i] - four*soln[((j+1)*x_nodes)+i] + six*soln[(j*x_nodes)+i] - four*soln[((j-1)*x_nodes)+i] + soln[((j-2)*x_nodes)+i]) * r_dydydydy;
         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;

         // Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 

         R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);
         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);
//printf("[Trace_Hit2] soln_new[38] = %2.15lf\n", soln_new[38]);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

         //Extrapolate pressure
         i = x_nodes - 1;
         //soln_new(i,j,1) = two*soln_new(i-1,j,1) - soln_new(i-2,j,1)
	 soln_new[(j*x_nodes)+i] = two*soln_new[(j*x_nodes)+(i-1)] - soln_new[(j*x_nodes)+(i-2)];

   }	// End For

#ifdef _OPENACC
//!$acc end kernels
#elif defined _OPENMP
//!$omp end parallel do
#endif

    L2[0] = L2[0] + L2_1;
    L2[1] = L2[1] + L2_2;
    L2[2] = L2[2] + L2_3;
    L2_1  = zero;
    L2_2  = zero;
    L2_3  = zero;

#ifdef _OPENACC
//!$acc kernels present(soln(:,:,:), soln_new(:,:,:)) 
#pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size]) 
//!$acc loop independent reduction(+:L2_1,L2_2,L2_3)
#pragma acc loop independent reduction(+:L2_1,L2_2,L2_3)
#elif defined _OPENMP
//!$omp parallel do private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
//#pragma omp parallel for private(i,j,dpdx,dudx,dvdx,dpdy,dudy,dvdy,d2udx2,d2vdx2,d2udy2,d2vdy2,uvel2,beta2,lambda_x,lambda_y,lambda_max,dtconv,dt,d4pdx4,d4pdy4,viscx,viscy,R1,R2,R3) reduction(+:L2_1, L2_2, L2_3)
#endif
   for(i = 2; i < x_nodes-2; i=i+1){

        j = 1;
        //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
  
         //uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		 
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)

         //d4pdx4 = ( soln(i+2,j,1) - four*soln(i+1,j,1) + six*soln(i,j,1) - four*soln(i-1,j,1) + soln(i-2,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+2)] - four*soln[(j*x_nodes)+(i+1)] + six*soln[(j*x_nodes)+i] - four*soln[(j*x_nodes)+(i-1)] + soln[(j*x_nodes)+(i-2)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j+3,1) - four*soln(i,j+2,1) + six*soln(i,j+1,1) - four*soln(i,j,1) + soln(i,j-1,1) )*r_dydydydy
         d4pdy4 = (soln[((j+3)*x_nodes)+ i] - four*soln[((j+2)*x_nodes)+i] + six*soln[((j+1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j-1)*x_nodes)+i]) * r_dydydydy;
         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;

         // Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 

         R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

         //Extrapolate pressure
         j = 0;
         //soln_new(i,j,1) = two*soln_new(i,j+1,1) - soln_new(i,j+2,1)         
	 soln_new[(j*x_nodes)+i] = two*soln_new[((j+1)*x_nodes)+i] - soln_new[((j+2)*x_nodes)+i];

         j = y_nodes-2;		 
         //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
  
           
         //uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		   
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
           
         //d4pdx4 = ( soln(i+2,j,1) - four*soln(i+1,j,1) + six*soln(i,j,1) - four*soln(i-1,j,1) + soln(i-2,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+2)] - four*soln[(j*x_nodes)+(i+1)] + six*soln[(j*x_nodes)+i] - four*soln[(j*x_nodes)+(i-1)] + soln[(j*x_nodes)+(i-2)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j-3,1) - four*soln(i,j-2,1) + six*soln(i,j-1,1) - four*soln(i,j,1) + soln(i,j+1,1) )*r_dydydydy         
         d4pdy4 = (soln[((j-3)*x_nodes)+ i] - four*soln[((j-2)*x_nodes)+i] + six*soln[((j-1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j+1)*x_nodes)+i]) * r_dydydydy;
         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;

         // Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 

         R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);
		 
         //Extrapolate pressure
         j = y_nodes - 1;
	 //soln_new(i,j,1) = two*soln_new(i,j-1,1) - soln_new(i,j-2,1)
	 soln_new[(j*x_nodes)+i] = two*soln_new[((j-1)*x_nodes)+i] - soln_new[((j-2)*x_nodes)+i];
            
   }	// End For

   
#ifdef _OPENACC
//!$acc end kernels
#elif defined _OPENMP
//!$omp end parallel do
#endif

    L2[0] = L2[0] + L2_1;
    L2[1] = L2[1] + L2_2;
    L2[2] = L2[2] + L2_3;
    //printf("S2#L2:%lf, %lf, %lf\n", L2[0], L2[1], L2[2]);

    L2_1   = zero;
    L2_2   = zero;
    L2_3   = zero;


//!$acc kernels present(soln(:,:,:), soln_new(:,:,:)) 
//!Corners-- (i,j) = (2,2), (2,y_nodes-1), (x_nodes-1,2), (x_nodes-1,y_nodes-1)
#pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size])  
{
    i = 1;
    j = 1;
         //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
        
	//uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		   
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         //lambda_max = max(lambda_x,lambda_y)                                    !Maximum eigenvalue
	 lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
                  
	 //d4pdx4 = ( soln(i+3,j,1) - four*soln(i+2,j,1) + six*soln(i+1,j,1) - four*soln(i,j,1) + soln(i-1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+3)] - four*soln[(j*x_nodes)+(i+2)] + six*soln[(j*x_nodes)+(i+1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i-1)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j+3,1) - four*soln(i,j+2,1) + six*soln(i,j+1,1) - four*soln(i,j,1) + soln(i,j-1,1) )*r_dydydydy
         d4pdy4 = (soln[((j+3)*x_nodes)+ i] - four*soln[((j+2)*x_nodes)+i] + six*soln[((j+1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j-1)*x_nodes)+i]) * r_dydydydy;
         //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;
                        
        //! Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 
	 R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

        //!Extrapolate pressure
        //soln_new(1,2,1) = two*soln_new(2,2,1) - soln_new(3,2,1)
	 soln_new[x_nodes] = two*soln_new[x_nodes+1] - soln_new[x_nodes+2];
        //soln_new(2,1,1) = two*soln_new(2,2,1) - soln_new(2,3,1)
	 soln_new[1] = two*soln_new[x_nodes+1] - soln_new[(2*x_nodes)+1];

        i = x_nodes-2;
        j = 1;
		
         //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
        
	//uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		   
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
		
	//d4pdx4 = ( soln(i-3,j,1) - four*soln(i-2,j,1) + six*soln(i-1,j,1) - four*soln(i,j,1) + soln(i+1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i-3)] - four*soln[(j*x_nodes)+(i-2)] + six*soln[(j*x_nodes)+(i-1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i+1)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j+3,1) - four*soln(i,j+2,1) + six*soln(i,j+1,1) - four*soln(i,j,1) + soln(i,j-1,1) )*r_dydydydy
         d4pdy4 = (soln[((j+3)*x_nodes)+ i] - four*soln[((j+2)*x_nodes)+i] + six*soln[((j+1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j-1)*x_nodes)+i]) * r_dydydydy;
        
	

        //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
        //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;
                        
        //! Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 
	 R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

        //!Extrapolate pressure
	//soln_new(x_nodes,2,1) = two*soln_new(x_nodes-1,2,1) - soln_new(x_nodes-2,2,1)
	 soln_new[x_nodes + (x_nodes-1)] = two*soln_new[x_nodes+(x_nodes-2)] - soln_new[x_nodes+(x_nodes-3)];
	//soln_new(x_nodes-1,1,1) = two*soln_new(x_nodes-1,2,1) - soln_new(x_nodes-1,3,1)
	 soln_new[x_nodes-2] = two*soln_new[x_nodes+(x_nodes-2)] - soln_new[(2*x_nodes)+(x_nodes-2)];


	 i = 1;
         j = y_nodes-2;
        //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
        
	//uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
        //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		   
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
		
	//d4pdx4 = ( soln(i+3,j,1) - four*soln(i+2,j,1) + six*soln(i+1,j,1) - four*soln(i,j,1) + soln(i-1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i+3)] - four*soln[(j*x_nodes)+(i+2)] + six*soln[(j*x_nodes)+(i+1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i-1)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j-3,1) - four*soln(i,j-2,1) + six*soln(i,j-1,1) - four*soln(i,j,1) + soln(i,j+1,1) )*r_dydydydy
         d4pdy4 = (soln[((j-3)*x_nodes)+ i] - four*soln[((j-2)*x_nodes)+i] + six*soln[((j-1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j+1)*x_nodes)+i]) * r_dydydydy;

        //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
        //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
        viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;
                        
        //! Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 
	R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
        //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

        //!Extrapolate pressure
	//soln_new(1,y_nodes-1,1) = two*soln_new(2,y_nodes-1,1) - soln_new(3,y_nodes-1,1)
	 soln_new[(y_nodes-2)*x_nodes] = two*soln_new[((y_nodes-2)*x_nodes) + 1] - soln_new[((y_nodes-2)*x_nodes) + 2];
	//soln_new(2,y_nodes,1) = two*soln_new(2,y_nodes-1,1) - soln_new(2,y_nodes-2,1)
	 soln_new[((y_nodes-1)*x_nodes) + 1] = two*soln_new[((y_nodes-2)*x_nodes) + 1] - soln_new[((y_nodes-3)*x_nodes) + 1];


        i = x_nodes-2;
        j = y_nodes-2;
        //dpdx = ( soln(i+1,j,1) - soln(i-1,j,1) )*r_twodx                  !2nd Order Central of pressure w.r.t. x
         dpdx = (soln[(j*x_nodes) + i+1] - soln[(j*x_nodes) + (i-1)]) * r_twodx;
         //dudx = ( soln(i+1,j,2) - soln(i-1,j,2) )*r_twodx                  !2nd Order Central of x velocity w.r.t. x
         dudx = (soln[(y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dvdx = ( soln(i+1,j,3) - soln(i-1,j,3) )*r_twodx                  !2nd Order Central of y velocity w.r.t. x
         dvdx = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i+1] - soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i-1]) * r_twodx;
         //dpdy = ( soln(i,j+1,1) - soln(i,j-1,1) )*r_twody                  !2nd Order Central of pressure w.r.t. y
         dpdy = (soln[((j+1)*x_nodes) + i] - soln[((j-1)*x_nodes) + i]) * r_twody;
         //dudy = ( soln(i,j+1,2) - soln(i,j-1,2) )*r_twody                  !2nd Order Central of x velocity w.r.t. y
         dudy = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //dvdy = ( soln(i,j+1,3) - soln(i,j-1,3) )*r_twody                  !2nd Order Central of y velocity w.r.t. y
         dvdy = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_twody;
         //d2udx2 = ( soln(i+1,j,2) - two*soln(i,j,2) + soln(i-1,j,2) )*r_dxdx   !2nd Order Central 2nd of x velocity w.r.t. x
         d2udx2 = (soln[(y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2vdx2 = ( soln(i+1,j,3) - two*soln(i,j,3) + soln(i-1,j,3) )*r_dxdx   !2nd Order Central 2nd of y velocity w.r.t. x
         d2vdx2 = (soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i+1)] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + (j*x_nodes) + (i-1)]) * r_dxdx;
         //d2udy2 = ( soln(i,j+1,2) - two*soln(i,j,2) + soln(i,j-1,2) )*r_dydy   !2nd Order Central 2nd of x velocity w.r.t. y
         d2udy2 = (soln[(y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
         //d2vdy2 = ( soln(i,j+1,3) - two*soln(i,j,3) + soln(i,j-1,3) )*r_dydy   !2nd Order Central 2nd of y velocity w.r.t. y
         d2vdy2 = (soln[(2*y_nodes*x_nodes) + ((j+1)*x_nodes) + i] - (two * soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i]) + soln[(2*y_nodes*x_nodes) + ((j-1)*x_nodes) + i]) * r_dydy;
        
	//uvel2 = soln(i,j,2)*soln(i,j,2) + soln(i,j,3)*soln(i,j,3)          !velocity squared at node
         uvel2 = (soln[(y_nodes*x_nodes)+(j*x_nodes)+i] * soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + (soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]);
         //beta2 = max(uvel2,rkappa*vel2ref)                                  !Beta squared at node
	 beta2 = max(uvel2, (rkappa*vel2ref));
		   
         //lambda_x = half*(fabs( soln(i,j,2) ) + sqrt( soln(i,j,2)*soln(i,j,2) + four*beta2 ))     !X eigenvalue
	 lambda_x = half*(fabs(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));
         //lambda_y = half*(fabs( soln(i,j,3) ) + sqrt( soln(i,j,3)*soln(i,j,3) + four*beta2 ))     !Y eigenvalue
	 lambda_y = half*(fabs(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]) + sqrt(soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] + four*beta2));

         lambda_max = max(lambda_x,lambda_y);
         dtconv = min(dx,dy)/lambda_max;                                         //Convective Limit for node
         dt = cfl*min(dtvisc,dtconv);                                            //Timestep for node (minimum)
		
	//d4pdx4 = ( soln(i-3,j,1) - four*soln(i-2,j,1) + six*soln(i-1,j,1) - four*soln(i,j,1) + soln(i+1,j,1) )*r_dxdxdxdx
	 d4pdx4 = (soln[(j*x_nodes)+(i-3)] - four*soln[(j*x_nodes)+(i-2)] + six*soln[(j*x_nodes)+(i-1)] - four*soln[(j*x_nodes)+i] + soln[(j*x_nodes)+(i+1)]) * r_dxdxdxdx;
	 //d4pdy4 = ( soln(i,j-3,1) - four*soln(i,j-2,1) + six*soln(i,j-1,1) - four*soln(i,j,1) + soln(i,j+1,1) )*r_dydydydy
         d4pdy4 = (soln[((j-3)*x_nodes)+ i] - four*soln[((j-2)*x_nodes)+i] + six*soln[((j-1)*x_nodes)+i] - four*soln[(j*x_nodes)+i] + soln[((j+1)*x_nodes)+i]) * r_dydydydy;

        //viscx = -((lambda_x*Cx*dx**3)/beta2)*d4pdx4                    !Set x artificial viscosity with d4pdx4
	 viscx = -((lambda_x*Cx*pow(dx,3))/beta2)*d4pdx4;
         //viscy = -((lambda_y*Cy*dy**3)/beta2)*d4pdy4                    !Set y artificial viscosity with d4pdy4
	 viscy = -((lambda_y*Cy*pow(dy,3))/beta2)*d4pdy4;
                        
        //! Calculate the new values of pressure (1), u velocity (2) and v velocity (3) at all interior nodes 
	 R1 = beta2*(rho*dudx+rho*dvdy-(viscx + viscy));
         //soln_new(i,j,1) = soln(i,j,1)-dt*R1
	 soln_new[(j*x_nodes) + i] = soln[(j*x_nodes) + i] - (dt*R1);

         //R2 = rhoinv*(rho*soln(i,j,2)*dudx+rho*soln(i,j,3)*dudy+dpdx-rmu*d2udx2-rmu*d2udy2)
	 R2 = rhoinv*(rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dudx+rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dudy+dpdx-rmu*d2udx2-rmu*d2udy2);
         //soln_new(i,j,2) = soln(i,j,2)-dt*R2
	 soln_new[(y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R2);

         //R3 = rhoinv*(rho*soln(i,j,3)*dvdy+rho*soln(i,j,2)*dvdx+dpdy-rmu*d2vdy2-rmu*d2vdx2)
	 R3 = rhoinv*(rho*soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i]*dvdy + rho*soln[(y_nodes*x_nodes)+(j*x_nodes)+i]*dvdx + dpdy - rmu*d2vdy2 - rmu*d2vdx2);
         //soln_new(i,j,3) = soln(i,j,3)-dt*R3
	 soln_new[(2*y_nodes*x_nodes)+(j*x_nodes)+i] = soln[(2*y_nodes*x_nodes)+(j*x_nodes)+i] - (dt*R3);

         L2_1   = L2_1 + pow(R1,2);
         L2_2   = L2_2 + pow(R2,2);
         L2_3   = L2_3 + pow(R3,2);

        //!Extrapolate pressure
	//soln_new(1,y_nodes-1,1) = two*soln_new(2,y_nodes-1,1) - soln_new(3,y_nodes-1,1)
	//soln_new(x_nodes,y_nodes-1,1) = two*soln_new(x_nodes-1,y_nodes-1,1) - soln_new(x_nodes-2,y_nodes-1,1)
	 soln_new[((y_nodes-2)*x_nodes) + (x_nodes-1)] = two*soln_new[((y_nodes-2)*x_nodes) + (x_nodes-2)] - soln_new[((y_nodes-2)*x_nodes) + (x_nodes-3)];
	//soln_new(x_nodes-1,y_nodes,1) = two*soln_new(x_nodes-1,y_nodes-1,1) - soln_new(x_nodes-1,y_nodes-2,1)
	 soln_new[((y_nodes-1)*x_nodes) + (x_nodes-2)] = two*soln_new[((y_nodes-2)*x_nodes) + (x_nodes-2)] - soln_new[((y_nodes-3)*x_nodes) + (x_nodes-2)];

         Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
        // printf("Pweightfactor= %lf\n", Pweightfactor);
//[D]         soln_new[(j*x_nodes) + i] = soln_new[(j*x_nodes) + i] - Pweightfactor;
//         printf("DEBUG## Pweightfactor = %lf\n", Pweightfactor);

        //!Pressure rescaling at the center point of the bottom floor
        //Pweightfactor = soln_new((x_nodes-x_nodes/2),1,1) - p_guage
//}       // End Scope of the acc kernel

//#pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size])
//        Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
	//soln_new(:,:,1)   = soln_new(:,:,1) - Pweightfactor
//#pragma acc kernels loop independent present(soln[array_init:array_size], soln_new[array_init:array_size])
//        #pragma acc loop independent vector(JTHREADS)
//	for(j=0; j<y_nodes; j=j+1){
//           #pragma acc loop independent vector(ITHREADS)
//	   for(i=0; i<x_nodes; i=i+1){
//	      soln_new[(j*x_nodes) + i] = soln_new[(j*x_nodes) + i] - Pweightfactor;
//	   }
//	}

}       // End Scope of the acc kernel

//#pragma acc wait(dbg == 1)
#pragma acc kernels present(soln[array_init:array_size], soln_new[array_init:array_size])
{
//   Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
#pragma acc loop independent vector(JTHREADS)
   for(y=0; y<y_nodes; y=y+1){
#pragma acc loop independent vector(ITHREADS)
      for(x=0; x<x_nodes; x=x+1){
//         Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
         soln_new[(y*x_nodes) + x] = soln_new[(y*x_nodes) + x] - Pweightfactor;
//         printf("Pweightfactor= %lf\n", Pweightfactor);
//         soln_new[(j*x_nodes) + i] = soln_new[(j*x_nodes) + i] - (soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage);
//           printf("DEBUG## Pweightfactor = %lf\n", Pweightfactor);
      }
   }
}

/*#pragma acc kernels present(soln_new[array_init:array_size])
{
   //Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
#pragma acc loop independent
   for(j=0; j<y_nodes; j=j+1){
#pragma acc loop independent
      for(i=0; i<x_nodes; i=i+1){
         Pweightfactor = soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage;
         soln_new[(j*x_nodes) + i] = soln_new[(j*x_nodes) + i] - Pweightfactor;
//         printf("Pweightfactor= %lf\n", Pweightfactor);
//         soln_new[(j*x_nodes) + i] = soln_new[(j*x_nodes) + i] - (soln_new[(x_nodes-(x_nodes/2)) - 1] - p_guage);
//           printf("DEBUG## Pweightfactor = %lf\n", Pweightfactor);
      }
   }
}*/


//!$acc end kernels

      L2[0] = sqrt(L2[0]) / internal_nodes;
//printf("[solver_c##Trace##] L[0] = %lf\n", L2[0]);
      L2[1] = sqrt(L2[1]) / internal_nodes;
      L2[2] = sqrt(L2[2]) / internal_nodes;
      //printf("S3#L2:%lf, %lf, %lf\n", L2[0], L2[1], L2[2]);

	
}

 void ldc_explicit( int x_nodes, int y_nodes, double dx, double dy, double* dt, double* beta, const double* __restrict__ soln, double* __restrict__ soln_new, int* iter_out ){

   
   double L1[3];
   double L2[3];
   double Linf[3];

    int iter, i, j, eq, array_init, array_size;

   array_init = 0;
   array_size = 3*x_nodes*y_nodes - 1;
    
//300     format(1X,i8,3e15.6)
//!$acc data copy(soln(:,:,:), soln_new(:,:,:)) 
//#pragma acc data copy(soln[array_init:array_size], soln_new[array_init:array_size])

    iter = 1;
    //do while (iter < max_iter)
#pragma acc data copy(soln[array_init:array_size], soln_new[array_init:array_size])
   while(iter < max_iter){
    //  !Iteration stage 1
    //  !Writes from soln to soln_new
       for( i=0; i<3; i=i+1){  
          L1[i]   = zero;
          L2[i]   = zero;
          Linf[i] = zero;
       }
      //printf("the x_nodes = %d, the y_nodes = %d, the dx = %lf, the dy = %lf\n", *x_nodes, *y_nodes, *dx, *dy);
      //printf("Iter: %d\n", iter);
      //printf("[Trace*] soln(3000,2001,2) = %2.15lf\n", soln[(x_nodes*y_nodes + 2000*x_nodes + 2999)]);
      //printf("[Trace*] soln[24982408] = %2.15lf\n", soln[24982408]);
      ldc_explicit_iter(x_nodes, y_nodes, dx, dy, soln, soln_new, L2 );
      //printf("\t%d, %e, %e, %e\n", iter, L2[0], L2[1], L2[2]);
      //printf("[Trace**] soln[38] = %2.15lf\n", soln[38]);

//      !Residual Calculations 1   
      if ( (iter%resid_out) == 0 )
        printf("\t%d, %lf, %lf, %lf\n", iter, L2[0], L2[1], L2[2]);
        // printf("\t%d, %e, %e, %e\n", iter, L2[0], L2[1], L2[2]);
      
      iter = iter + 1;

      //!Iteration stage 2
      //!Writes from soln_new to soln
      for( i=0; i<3; i=i+1){  
          L1[i]   = zero;
          L2[i]   = zero;
          Linf[i] = zero;
      }
      //printf("[Trace*1] soln[38] = %2.15lf\n", soln[38]);	   
      ldc_explicit_iter( x_nodes, y_nodes, dx, dy, soln_new, soln, L2 );
      //printf("\t%d, %e, %e, %e\n", iter, L2[0], L2[1], L2[2]);
      //printf("[Trace**1] soln[38] = %2.15lf\n", soln[38]);

      //!Residual Calculations 2
      if ((iter%resid_out) == 0 )
         printf("\t%d, %lf, %lf, %lf\n", iter, L2[0], L2[1], L2[2]);
        // printf("\t%d, %e, %e, %e\n", iter, L2[0], L2[1], L2[2]);

      iter = iter + 1;


   }	// End While

//!$acc end data
//soln = soln_new 
//cudaDeviceSynchronize();

   for (i = 0; i < (3*y_nodes*x_nodes); i=i+1){
      soln_new[i] = soln[i];
   }

   *iter_out = iter;
}
