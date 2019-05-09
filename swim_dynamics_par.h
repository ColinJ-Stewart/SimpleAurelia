#ifndef SWIM_DYNAMICS_H
#define SWIM_DYNAMICS_H

#include <udf.h>
/* #include <unsteady.h> */

#include "global_var_par.h"
#include "timing_and_counters_par.h"	/* needed for NewTimeStepForThisZone() */
#include "mesh_motion_par.h"	 /* needed for Rezero() in Get_BodyVelocity() */ 

/* variables */
/* extern double del_x;*/
/* extern double del_xS;*/
/* extern double nextVel;*/
/* extern double nextVelS;*/

/* prototypes */
double Get_Force_x(Thread *tf);
double Get_Force2_x(Thread *tf);
double Get_BodyVelocity(Thread *tf, char meshZone, double fx_tot, double *del_x_ptr);
double Get_Mass(void);
double Get_Vel_At_End(char meshZone);
double Get_Del_At_End(char meshZone);
double Get_BodyVelocity_Implicit(Thread *tf, char meshZone, double fx_tot, double *del_x_ptr);
/* double Rezero (Thread *tf); */

#endif
