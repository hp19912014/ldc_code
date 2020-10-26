#include <stdio.h>
#include <stdlib.h>
/*Comments*/
#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_
/* Place to put all of my definitions etc. */
#include "global_vars.h"

#endif
//#include "global_vars.h"


void read_input(){
   double tmp_1;
   //bool exit_loop = 0;
   char filler_buf[100];		// Holds unneeded data -- Used to advance the fscanf cursor over file fillers (e.g. skip "domain" line)
   char file_name[20] = "ldc_c.nml";
   //char line[100];
   FILE* fp;
   
   fp = fopen(file_name, "r");
   
   if (fp == NULL){
      printf("[Trace#fileio#] Error While Opening the (ldc.nml) file ... \n");
	  exit(1);
   }
   
   /////////////////////////////// Parsing LDC.NML File ///////////////////////////
//   printf("[Trace#fileio#] Start of parsing the Input (ldc.nml) file ... \n");
   fscanf(fp, "%s", filler_buf);	// skip "domain"
   
   fscanf(fp, "%s", filler_buf);	// skip "x_nodes" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%d", &x_nodes);		// reading the x_node value"
//   printf("[Trace#fileio#] x_nodes = %d\n", x_nodes);
   
   fscanf(fp, "%s", filler_buf);	// skip "y_nodes" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%d", &y_nodes);
//   printf("[Trace#fileio#] y_nodes = %d\n", y_nodes);
   
   fscanf(fp, "%s", filler_buf);	// skip "length" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &length);
//   printf("[Trace#fileio#] length = %lf\n", length);
   
   fscanf(fp, "%s", filler_buf);	// skip "xmin" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &xmin);
//   printf("[Trace#fileio#] xmin = %lf\n", xmin);
   
   fscanf(fp, "%s", filler_buf);	// skip "xmax" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &xmax);
//   printf("[Trace#fileio#] xmax = %lf\n", xmax);
   
   fscanf(fp, "%s", filler_buf);	// skip "ymin" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &ymin);
//   printf("[Trace#fileio#] ymin = %lf\n", ymin);
   
   fscanf(fp, "%s", filler_buf);	// skip "ymax" 
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &ymax);
//   printf("[Trace#fileio#] ymax = %lf\n", ymax);
   
   fscanf(fp, "%s", filler_buf);	// skip "physical_properties"
   
   fscanf(fp, "%s", filler_buf);	// skip "re"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &re);
//   printf("[Trace#fileio#] re = %lf\n", re);
   
   fscanf(fp, "%s", filler_buf);	// skip "rho"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &rho);
//   printf("[Trace#fileio#] rho = %lf\n", rho);
   
   fscanf(fp, "%s", filler_buf);	// skip "u_lid"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &u_lid);
//   printf("[Trace#fileio#] u_lid = %lf\n", u_lid);
   
   fscanf(fp, "%s", filler_buf);	// skip "p_guage"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &p_guage);
//   printf("[Trace#fileio#] p_guage = %lf\n", p_guage);
   
   fscanf(fp, "%s", filler_buf);	// skip "solver_properties"
   
   fscanf(fp, "%s", filler_buf);	// skip "solver"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%s", solver);
//   printf("[Trace#fileio#] solver = %s\n", solver);
   
   fscanf(fp, "%s", filler_buf);	// skip "max_iter"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%d", &max_iter);
//   printf("[Trace#fileio#] max_iter = %d\n", max_iter);
   
   fscanf(fp, "%s", filler_buf);	// skip "cfl"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &cfl);
//   printf("[Trace#fileio#] cfl = %lf\n", cfl);
   
   fscanf(fp, "%s", filler_buf);	// skip "k"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &k);
//   printf("[Trace#fileio#] k = %lf\n", k);
   
   fscanf(fp, "%s", filler_buf);	// skip "c2"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &c2);
//   printf("[Trace#fileio#] c2 = %lf\n", c2);
   
   fscanf(fp, "%s", filler_buf);	// skip "Cx"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &Cx);
//   printf("[Trace#fileio#] Cx = %lf\n", Cx);
   
   fscanf(fp, "%s", filler_buf);	// skip "Cy"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &Cy);
//   printf("[Trace#fileio#] Cy = %lf\n", Cy);
   
   fscanf(fp, "%s", filler_buf);	// skip "rkappa"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &rkappa);
//   printf("[Trace#fileio#] rkappa = %lf\n", rkappa);
   
   fscanf(fp, "%s", filler_buf);	// skip "conv_toler"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &conv_toler);
//   printf("[Trace#fileio#] conv_toler = %1.10f\n", conv_toler);
   
   fscanf(fp, "%s", filler_buf);	// skip "visc_eps"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%lf", &visc_eps);
//   printf("[Trace#fileio#] visc_eps = %1.10lf\n", visc_eps);
   
   fscanf(fp, "%s", filler_buf);	// skip "output"
   
   resid_out = 1;					// It is just default, but it is overwritten by the file reading below
   fscanf(fp, "%s", filler_buf);	// skip "resid_out"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%d", &resid_out);
//   printf("[Trace#fileio#] resid_out = %d\n", resid_out);
   
   fscanf(fp, "%s", filler_buf);	// skip "write_solution"
   fscanf(fp, "%s", filler_buf);	// skip "="
   fscanf(fp, "%d", &write_solution);
//   printf("[Trace#fileio#] write_solution = %d\n", write_solution);
   /************************** End Reading Variables & File **********************/
   
   fclose(fp);
}

void write_output(int x, int y, double* soln){

   int i, j;
   char file_name[20] = "ldc.dat";
   FILE* fp;
   
   fp = fopen(file_name, "w");
   if (fp == NULL){
      printf("[Trace#fileio#] Error While Opening the (ldc.dat) file ... \n");
	  exit(1);
   }
   
//   printf("TITLE=\"Lid Driven Cavity Solution\"\n");
//   printf("VARIABLES = \"X\", \"Y\", \"Pressure\", \"U\", \"V\"\n");
//   printf("ZONE DATAPACKING=POINT, I=%d, J=%d\n", x_nodes, y_nodes);
   
   for(j = 0; j < y_nodes; j = j+1){
      for(i = 0; i < x_nodes; i = i+1){
	     printf("%lf %lf %lf %lf %lf\n", ((xmax-xmin)*(i-1)/(x_nodes-1)), ((ymax-ymin)*(i-1)/(y_nodes-1)),
		        soln[(y_nodes*x_nodes) + (j*x_nodes) + i], 
		        soln[(2*y_nodes*x_nodes) + (j*x_nodes) + i],
		        soln[(3*y_nodes*x_nodes) + (j*x_nodes) + i]);
	  }
   }
   
   fclose(fp);   

}

