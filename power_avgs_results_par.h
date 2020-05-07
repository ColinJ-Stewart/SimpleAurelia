#ifndef POWER_AVGS_RESULTS_H
#define POWER_AVGS_RESULTS_H

#include <udf.h>
/* #include <unsteady.h> */

#include "global_var_par.h"		
#include "swim_dynamics_par.h"			/* needed for Get_Mass() */
#include "timing_and_counters_par.h"	/* needed for Get_CycleNumber(), New_Cycle(), isPaused */

/* prototypes */
double Compute_Power_In(Thread *tf, double velx, double vely);
double Compute_Power_Loss(Thread *tf, double velx, double vely);
void Compute_Power_Out(double velx, double vely, 
						double *f_thrust_SMC_ptr, double *f_thrust_B_ptr, 
						double *f_drag_SMC_ptr, double *f_drag_B_ptr,
						double *P_thrust_SMC_ptr, double *P_thrust_B_ptr,
						double *P_drag_SMC_ptr, double *P_drag_B_ptr);
void Compute_Averages(double next_vel, double f_net, 
						double f_thrust_SMC, double f_thrust_B, double f_drag_SMC, 
						double f_drag_B, double P_in, double P_loss, double P_thrust_SMC, 
						double P_thrust_B, double P_drag_SMC, double P_drag_B);
void Create_Disp_OutputFile(char side);
void Write_to_Disp_OutputFile(char side, double x, double vel, double del_x, 
	double f_net, double f_thrust_SMC, double f_thrust_B, double f_drag_SMC, double f_drag_B, 
	double P_in, double P_thrust_SMC, double P_thrust_B, double P_drag_SMC, double P_drag_B,
	double P_loss);

#endif
