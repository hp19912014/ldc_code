
#include "omp.h"
#include "global_vars.h"

int main(){

   double dtime, starttime, endtime;
   int niter;
   #pragma acc host_data
   {
   cudaDeviceSetCacheConfig(2);
   }   
//   printf("Reading input file LDC.nml...\n");
   read_input();
//   printf("Allocating soln, soln_new, beta and bt arrays...\n");
   dt = malloc(y_nodes*x_nodes*sizeof(double));
   beta = malloc(y_nodes*x_nodes*sizeof(double));
   soln = malloc(3*y_nodes*x_nodes*sizeof(double));
   soln_new = malloc(3*y_nodes*x_nodes*sizeof(double));
//   double *dt,*beta,*soln,*soln_new;
//   cudaMallocManaged(&dt,y_nodes*x_nodes*sizeof(double));
//   cudaMallocManaged(&beta,y_nodes*x_nodes*sizeof(double));
//   cudaMallocManaged(&soln,3*y_nodes*x_nodes*sizeof(double));
//   cudaMallocManaged(&soln_new,3*y_nodes*x_nodes*sizeof(double));
//   cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
//   printf("Initializing the arrays...\n");
     ldc_allocate(soln, soln_new,dt, beta);
//   printf("Arrays Are Ready! ...\n");
//   printf("------------------------------------------------\n");
  
   starttime = omp_get_wtime();  
//   printf("[Trace MSG] ## solver = %s\n", solver);
   
   if (strncmp(solver,"explicit", 8) == 0){
//      printf("Beginning explict solve...\n");
      ldc_explicit(x_nodes, y_nodes, dx, dy, dt, beta, soln, soln_new, &niter);
   }
   else if (strncmp(solver,"implicit",8) == 0){
//      printf("Beginning implict solve...\n");
      //ldc_implicit(x_nodes, y_nodes, dx, dy, dt, beta, soln, soln_new);
   }
   else if (strncmp(solver,"sgs",3) == 0){
//      printf("Beginning SGS solve...\n");
      //ldc_sgs(x_nodes, y_nodes, dx, dy, dt, beta, soln);
   }
   else{
//      printf("Unknown solver type... Freeing Arrays & stopping\n");
	  free(dt);
	  free(beta);
	  free(soln);
	  free(soln_new);
	  exit(1);
   }
   endtime = omp_get_wtime();
   
   // Calculate solver runtime
   dtime = endtime - starttime;

   if (write_solution==1){
//      write_output(x_nodes, y_nodes, soln);
   }

   // TODO: Freeing the Arrays instead of the main
   ldc_deallocate();
   free(dt);
   free(beta);
   free(soln);
   free(soln_new);

   //Display solver runtime
//   printf("Solver runtime(s) = %lf\n", dtime);
    printf("%lf\n",dtime);
     
   // Display gflops

//   printf("x:%d,y:%d,niter:%d,dtime:%f\n",x_nodes,y_nodes,niter,dtime);
   double tmp;
   tmp=130.0*x_nodes*y_nodes/dtime;
//   printf("tmp:%f\n",tmp);
   tmp=tmp/1000000;
     printf("%f\n",tmp);

// printf("GFLOPS=%f\n",  (( x_nodes*y_nodes*niter*130.0) / (dtime*1000000000.0)) ); 
  
   return 0;
}

