/*****************************************************************
Title: SimpleAurelia_4b_parallel.c

Purpose: Aurelia aurita mesh motion UDF with kinematics determined by functions similar
	to those found in Herschlag and Miller, 2011. This is for use on a follow up paper to Gemmell
	et al 2013 (PNAS), and specifically addresses points raised in Gemmell et al 2018 (JEB).
	
Date: 2/20/2019
Author: Colin Stewart
Program versions: 
	OS: W7x64
	Compiler: MS Visual Studio 2010 Express
	Fluent: 14.0

Input: 
	Text files with cubic spline coefficients for kinematic functions sa, sb, sd, and 
	margin flexibility function sg:
		SimpleAureliaSplineCoeff_sa.txt
		SimpleAureliaSplineCoeff_sb.txt
		SimpleAureliaSplineCoeff_sd.txt
		SimpleAureliaSplineCoeff_sg.txt
	
	Cubic spline functions are defined by:
		s(t) = C_i(t) = c3_i*t^3 + c2_i*t^2 + c1_i*t^1 + c0_i  		on t_(i-1) < t <= t_i, where i = 1,2,...n

	 Data is formatted as follows:
	 4 header rows
	 c3_1 c2_1 c1_1 c0_1 
	 c3_2 c2_2 c1_2 c0_2 
	 ...
	 c3_numSplines c2_numSplines c1_numSplines c0_numSplines 
	 
	 i.e., Row order: spline1 to numSplines, Column order: c3 c2 c1 c0 


	
Change log:
	v1 -- (11/14/2018) Started with the code from UDF file "SarsiaMiller_swim_velumPivot_2PIVOT.c". 
			Velar kinematics and related code will be removed. The improvements to the old Aurelia (pause) simulation
			UDF will be imported here when relevent. Notably, need to implement (1) pause, (2) power calculation,
			(3) implicit solving of EOM.

	v2 -- picking up where I left off before APS, so I felt like it needed a new version number even though v1 doesn't
			work yet! This version can succesfully move w/ flexibility, but cannot swim.
	
	v3 -- fixing code so the tip node is included in the calculation of arc length for both top and bottom. Also,
			trying refined meshes w/ boundary layer refinement.
	
	v4 -- Fixed tip flex. Calculations (e.g. arc length) now based on unflexed profile.
	
	v5 -- Efficiency equations fixed. Many patches to handle remeshed nodes (assigning tparm, etc.).
	
	v6 -- Fixed pause part of code (was not passing tauSim to moveNodes before).
	
	v7 -- Split code in several libraries and header files.
	
	v8 -- Parallelized UDF
	
			
***************************************************************** */
#include <udf.h>
#include <unsteady.h>

#include "global_var_par.h"
#include "load_kinematics_par.h"
#include "mesh_motion_par.h"
#include "power_avgs_results_par.h"
#include "quicksort_par.h"
#include "swim_dynamics_par.h"
#include "timing_and_counters_par.h"
#include "user_mem_mgmt_par.h"

#define NUM_INPUT_FILES 6
#define BUFFER_SIZE 1024

/* Simulation settings  */
const int LET_IT_SWIM = 1;				/* flag to allow the jellyfish to swim */
const double PAUSE_FRAC = 0.0;			/* creates a pause after relaxation = PAUSE_FRAC * PERIOD */

/* Model animal measurements  */
const double DIAMETER = 0.155;									/* relaxed bell diameter [m] */
const double thickness = 0.155 * (53.0/450.0);				/* relaxed bell center axis thickness [m]. The second value is the ratio of thickness [px] to D [px] in Dabiri 2011 Aurelia video */
const double t_contract = (66.0-42.0)/30.0; 				/* contraction duration [s] */
const double t_relax = (128.0-66.0)/30.0; 					/* relaxation duration [s] */
const double PERIOD = (66.0-42.0)/30.0 + (128.0-66.0)/30.0;	/* period = 2.86667 [s] */

/* Fluid properties  */
const double NU = 1e-6;								/* [m^2/s] */
const double RHO = 1000;							/* [kg/m^3] */

int counter_top = 0;
int counter_bot = 0;

#if !RP_HOST							/* Compiler directive for parallel computing */
Domain *domain;
#endif

/*************************************************************** */
/*********************** Macros ******************************** */
/*************************************************************** */

/********************** Get domain ***************************** */
DEFINE_EXECUTE_ON_LOADING(get_domain, libname)
{
	#if !RP_HOST
		domain = Get_Domain(1);
		if(domain == NULL)
			Message("Something wrong with your domain id!\n");
		else
			Message("\n Domain set.\n");
	#endif
}
	
DEFINE_INIT(init_node_mem, domain)
{
	int total_bytes_copied;
	int i;
	char buffer[NUM_INPUT_FILES][BUFFER_SIZE];
	
	
	/* Initialize memory  */
	#if !RP_HOST
		Thread *tf_ex = Lookup_Thread(domain, EX_ZONE);			/* might not be needed with the getForces function */
		Thread *tf_sub = Lookup_Thread(domain, SUB_ZONE);
		
		Init_node_mem(tf_ex);
		Init_node_mem(tf_sub); 
		
		Message("Memory slots initialized.\n");
	#endif
	
	
	/* Copy input files to all nodes  */
	#if PARALLEL
		sprintf(buffer[0], "SimpleAureliaSplineCoeff_s%c.txt", 'a');
		sprintf(buffer[1], "SimpleAureliaSplineCoeff_s%c.txt", 'b');
		sprintf(buffer[2], "SimpleAureliaSplineCoeff_s%c.txt", 'd');
		sprintf(buffer[3], "SimpleAureliaSplineCoeff_s%c.txt", 'g');
		sprintf(buffer[4], "SimpleAureliaSplineCoeff_asub.txt");
		sprintf(buffer[5], "SimpleAureliaSplineCoeff_bsub.txt");

		for (i = 0; i < NUM_INPUT_FILES; i++)
		{
			#if RP_HOST
				total_bytes_copied = host_to_node_sync_file(buffer[i]);
			#endif
			
			#if RP_NODE
				mkdir("/tmp", 0700);
				total_bytes_copied = host_to_node_sync_file("/tmp");
			#endif
			Message("Total bytes copied by node %i is %d\n", myid, total_bytes_copied);
		}
	#endif
	
}

/********************** Jelly motion ***************************** */
/* Should be called by both exumbrella and subumbrella zones for  */
/* dynamic mesh motion (both kinematic and swimming). */
/*************************************************************** */
DEFINE_GRID_MOTION(Jelly_motion,domain,dt,time,dtime)
{
	#if RP_HOST
		/* --------- INITIALIZE VARS ---------  */
		char meshZone = '\0';
		int meshZoneTemp;
		
		/* ------ PRINT TIME STEP HEADER ---------  */
		Print_Timestep_Header();
	#endif
	
	
	#if !RP_HOST
		/* --------- INITIALIZE VARS ---------  */
		Thread *tf;
		Thread *tf_ex;
		Thread *tf_sub;
		double	tauSim, fx_ex, fx_sub, fx_tot;
		double vel = 0.0, del_x = 0.0;
		double *del_x_ptr = &del_x;
		char meshZone;
		int meshZoneTemp;
		double Sat;
		int zone_ID;
		
		/* Message("Node%i initializing vars... \n", myid); */
		tf = DT_THREAD(dt);
		zone_ID = THREAD_ID(tf);
		tf_ex = Lookup_Thread(domain, EX_ZONE);
		tf_sub = Lookup_Thread(domain, SUB_ZONE);
		SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));				/* MUST CALL to set deforming flag on adjacent cell zone  */
		
		/* --------- MESH ZONE INFO ---------- */ 
		/* Message("Node%i getting mesh zone... \n", myid); */
		if (zone_ID == EX_ZONE) meshZone = 'e';
		else if (zone_ID == SUB_ZONE) meshZone = 's';
		else Error("Mesh zone not defined properly!\n");
		
		/* ------ PRINT HEADERS ---------  */
		/* Message("Node%i printing timestep header... \n", myid); */
		#if !PARALLEL
			Print_Timestep_Header();
		#endif
		Print_Zone_Header(meshZone);
		
		/* ------------- MEMORY ------------- */ 
		if ( NewTimeStep() )
		{
			Message("New time step. Reinitializing node memory\n");
			Reinit_node_mem_int(tf_ex, 0, 0);	/* Re-Initialize node memory slot 0 to 0 */
			Reinit_node_mem_int(tf_sub, 0, 0);	/* Re-Initialize node memory slot 0 to 0 */
		}
			
		/* ------------- TIME  --------------  */	
		tauSim = Get_SimTimeInCycle(meshZone);	/* SIMULATED time [s] within the cycle, i.e. in the range [0, period). Pause time will cause the simulated time to stop */

		/* ------------ FORCES -------------  */
		/* fx_ex = Get_Force_x(tf_ex); */			
		/* fx_sub = Get_Force_x(tf_sub); */
		fx_ex = Get_Force2_x(tf_ex);			
		fx_sub = Get_Force2_x(tf_sub);
		fx_tot = fx_ex + fx_sub;
		if ( NewTimeStepForThisZone(meshZone) )
			Message("Force_x: %lf, fx_ex: %lf, fx_sub: %lf\n", fx_tot, fx_ex, fx_sub);	
		
		/* ------ SWIMMING DYNAMICS -------  */
		if((LET_IT_SWIM) && (N_TIME > 2)) 
		{
			if ( NewTimeStepForThisZone(meshZone) )
				Message("Calculating swimming speed... \n");	
			vel = Get_BodyVelocity(tf, meshZone, fx_tot, del_x_ptr);
			if ( NewTimeStepForThisZone(meshZone) )
				Message("\tvel_x = %lf, del_x = %lf\n", vel, del_x);	/** make Message0 after DEBUG **/
		}

		/* ------ MESH MOTION/KINEMATICS -------  */
		if ( NewTimeStepForThisZone(meshZone) )
			Message("Moving jellyfish (node %i)... \n", myid);
		Calc_Mesh_Movement(tf, meshZone, tauSim, del_x, vel);
	#endif 
	
	
	/* send the meshzone info from node0 to host  */
	meshZoneTemp = meshZone;
	node_to_host_int_1(meshZoneTemp);
	meshZone = meshZoneTemp;
	
	
	#if RP_HOST
		/* ------ PRINT ZONE HEADER ------  */
		Print_Zone_Header(meshZone);
	#endif
	
	#if !RP_NODE
		/* ------ INITIALIZE RESULT FILES -------  */
		/* If first iteration, create result files, write out headers, etc.  */
		
		/* Message("N_TIME = %i, NewTimeStepForThisZone = %i\n", N_TIME, NewTimeStepForThisZone(meshZone)); */
		if((N_TIME == 0) && (NewTimeStepForThisZone(meshZone))) 
		{
			Message("Creating result files\n");
			Create_Disp_OutputFile(meshZone);
			
			/* write out the results at t = 0s because the DEFINE_EXECUTE_AT_END macro isn't called until step 1  */
			Message("Writing to result files\n");
			Write_to_Disp_OutputFile(meshZone, 
									0, 0, 0, 0,
									0, 0, 0, 0, 
									0, 0, 0, 0, 
									0);
		}
	#endif
	
	
	/* ------ STORE THIS TIME STEP # -------  */
	/* Message("Updating counters... "); */
	UpdateCounters(meshZone);
	/* Message("Done.\n"); */
	if ( NewTimeStepForThisZone(meshZone) )
		Message("---FINISHED GRID MOTION UDF---\n");
}


/********************** axis_top ***************************** */
/* Slide the two axis zones, axis_top = x < jellyfish and  */
/* axis_bot = x > jellyfish, along the x-axis.  */
/* NOTE: This can simply be done by specifying the deforming  */
/* motion as "plane" with normal in (0,1) direction and (0,0)  */
/* defining a point on the plane. */
/*************************************************************** */
DEFINE_GEOM(axis_top, domain, dt, position)
{
	/* #if !RP_HOST */
	/* int counter_top = 0; */
	/* Thread *tf = DT_THREAD(dt); */
	/* int zone_ID = THREAD_ID(tf); */
	/* Thread *tf_ex = Lookup_Thread(domain, EX_ZONE);			// might not be needed with the getForces function */
	/* Thread *tf_sub = Lookup_Thread(domain, SUB_ZONE); */
	
	/* Message("(node%i) Call %i\n", myid, ++counter_top); */
    position[1] = 0;
	/* #endif */
}

/************************ axis_bot ***************************** */
DEFINE_GEOM(axis_bot, domain, dt, position)
{
	/* #if !RP_HOST */
	/* int counter_bot = 0; */
	/* Message("(node%i) Call %i\n", myid, ++counter_bot); */
    position[1] = 0;
	/* #endif */
}

/************************ CG_motion_intfc **********************/
/* Move the interfaces between nearfield (tet) and farfield (hex) meshes */
/****************************************************************/
DEFINE_CG_MOTION(CG_motion_intfc, dt, vel, omega, time, dtime)
{
	double bodyVel = 0.0;
	double del_x = 0.0;
	double *del_x_ptr = &del_x;
	char meshZone = 'i';
	
	#if !RP_HOST
		Thread *tf_ex	= Lookup_Thread(domain, EX_ZONE);
		Thread *tf_sub	= Lookup_Thread(domain, SUB_ZONE);
		
		/* Calculate forces */
		double fx_ex 	= Get_Force2_x(tf_ex);					/* net force across jelly exum [N] */
		double fx_sub 	= Get_Force2_x(tf_sub);		/* net force across jelly sub [N] */
		double fx_tot = fx_ex + fx_sub;
		
		/* Calculate current position, body velocity */
		if((LET_IT_SWIM) && (N_TIME > 2)) 
		{
			bodyVel = Get_BodyVelocity(tf_ex, 'i', fx_tot, del_x_ptr);
		}
		
		vel[0] = bodyVel;
		UpdateCounters('i');
	
	#endif
}

/********************** Forces_at_end *************************** */
/* Compute and write out results at end of time step. This includes */
/* force, power, averages, and efficiencies.  */
/*************************************************************** */
DEFINE_EXECUTE_AT_END(forces_at_end)
{
	
	/* Initialize/calc variables on SERIAL or NODE processes  */
	#if !RP_HOST  
		Thread *tf_ex_end = Lookup_Thread(domain, EX_ZONE);
		Thread *tf_sub_end = Lookup_Thread(domain, SUB_ZONE);
		
		double x = Rezero(tf_ex_end);
		double xS = Rezero(tf_sub_end);
		
		double vel = Get_Vel_At_End('e');
		double velS = Get_Vel_At_End('s');
		double del_x = Get_Del_At_End('e');
		double del_xS = Get_Del_At_End('s');
		
		double fx_ex = Get_Force2_x(tf_ex_end);				/* net force on exumbrellar (top) surface [N] */
		double fx_sub = Get_Force2_x(tf_sub_end);			/* net force on subumbrellar (bottom) surface [N] */
		double fx_tot = fx_ex + fx_sub;
		
		double P_in_ex = Compute_Power_In(tf_ex_end, vel, 0.0);
		double P_in_sub = Compute_Power_In(tf_sub_end, velS, 0.0);
		double P_in = P_in_ex + P_in_sub;
		
		double f_thrust_SMC, f_thrust_B;
		double f_drag_SMC, f_drag_B;
		double P_thrust_SMC, P_thrust_B;
		double P_drag_SMC, P_drag_B;
		
		double *f_thrust_SMC_ptr = &f_thrust_SMC;
		double *f_thrust_B_ptr = &f_thrust_B;
		double *f_drag_SMC_ptr = &f_drag_SMC;
		double *f_drag_B_ptr = &f_drag_B;
		
		double *P_thrust_SMC_ptr = &P_thrust_SMC;
		double *P_thrust_B_ptr = &P_thrust_B;
		double *P_drag_SMC_ptr = &P_drag_SMC;
		double *P_drag_B_ptr = &P_drag_B;
		/* Message0("\n---NODE%i EXECUTE AT END---\n", myid); */
		Message0("Force_x at end: %f\n",fx_tot);	
		
		Compute_Power_Out(vel, 0.0, f_thrust_SMC_ptr, 	f_thrust_B_ptr, f_drag_SMC_ptr, 	f_drag_B_ptr, 
									P_thrust_SMC_ptr,	P_thrust_B_ptr, P_drag_SMC_ptr, 	P_drag_B_ptr);
		#endif 
	
	
	/* Initialize variables on HOST to be sent from NODE0  */
	#if RP_HOST
		double x, vel, del_x;
		double xS, velS, del_xS;
		double fx_tot, f_thrust_SMC, f_thrust_B, f_drag_SMC, f_drag_B;
		double P_in, P_thrust_SMC, P_thrust_B, P_drag_SMC, P_drag_B;
		Message("---HOST EXECUTE AT END---\n");
	#endif
	
	
	/* Commands for ALL (i.e. SERIAL, NODE, and HOST)  */
	/* - Send data from node0 to host for printing 	   */
	node_to_host_real_3(x, vel, del_x);
	node_to_host_real_3(xS, velS, del_xS);
	node_to_host_real_5(fx_tot, f_thrust_SMC, f_thrust_B, f_drag_SMC, f_drag_B);
	node_to_host_real_5(P_in, P_thrust_SMC, P_thrust_B, P_drag_SMC, P_drag_B);
	
	
	/* SERIAL/HOST writes out results for instant displacement, vel, force, power, etc.  */
	#if !RP_NODE
		Compute_Averages(vel, fx_tot, 	f_thrust_SMC, 	f_thrust_B, f_drag_SMC, 	f_drag_B,
								P_in, 	P_thrust_SMC, 	P_thrust_B, P_drag_SMC, 	P_drag_B);
		
		Write_to_Disp_OutputFile('e', x, vel, del_x, 	
								fx_tot, 	f_thrust_SMC, 	f_thrust_B, f_drag_SMC, 	f_drag_B,
								P_in, 	P_thrust_SMC, 	P_thrust_B, P_drag_SMC, 	P_drag_B);
													
		Write_to_Disp_OutputFile('s', xS, velS, del_xS,
								fx_tot, 	f_thrust_SMC, 	f_thrust_B, f_drag_SMC, 	f_drag_B,
								P_in, 	P_thrust_SMC, 	P_thrust_B, P_drag_SMC, 	P_drag_B);							
	#endif
	
	
	/* Commands for ALL (i.e. SERIAL, NODE, and HOST)  */
	if (NewCycle()) Store_OldCycleN(); 	/* i.e. cycleN_prev = cycleN; */
	
	
	/* SERIAL/NODES store the old values of tparm, tip coords, etc.  */
	#if !RP_HOST
		Store_OldKinematicVars();
	#endif
	
	Message0("---FINISHED EXECUTE AT END UDF---\n");
}



