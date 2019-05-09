#include "power_avgs_results_par.h"

/** global variables * */
#if !RP_NODE	/* SERIAL or HOST */
static double vel_tot = 0.0;				/* running sum of velocity in this cycle [m/s] */
static double f_tot = 0.0;					/* running sum of net force [N] */
static double f_thrust_SMC_tot = 0.0;		/* running sum of thrust force [N] as defined by Sahin, Mohseni, and Colin */
static double f_thrust_B_tot = 0.0;			/* running sum of thrust force [N] as defined by Borazjani et al */
static double f_drag_SMC_tot = 0.0;			/* running sum of drag (friction) force [N] as defined by Sahin, Mohseni, and Colin */
static double f_drag_B_tot = 0.0;			/* running sum of drag force [N] as defined by Borazjani et al */
static double P_in_tot = 0.0;				/* running sum of propulsive power input [W] */
static double P_thrust_SMC_tot = 0.0;		/* running sum of thrust power [W] as defined by Sahin, Mohseni, and Colin */
static double P_thrust_B_tot = 0.0;			/* running sum of thrust power [W] as defined by Borazjani et al */
static double P_drag_SMC_tot = 0.0;			/* running sum of drag (friction) power [W] as defined by Sahin, Mohseni, and Colin */
static double P_drag_B_tot = 0.0;			/* running sum of drag power [W] as defined by Borazjani et al */
static int kavg = 0;						/* number of measurements/timesteps */
#endif

/** Function definitions * */
/*---------- Compute_Power_In ------------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
double Compute_Power_In(Thread *tf, double velx, double vely) 
{
	double power = 0.0;
	
	#if !RP_HOST	/*SERIAL or NODE */
	Thread *tf0;
	face_t f;
	cell_t c;
	double A[ND_ND];
	double powerInviscid = 0.0;
	double powerDebug = 0.0;
	double thrustDebugInv = 0.0;
	double thrustDebugVisc1 = 0.0;
	double thrustDebugVisc2 = 0.0;
	

	/* loop through cells on  surface, calculating power from each  */
	begin_f_loop(f,tf) {  
		c = F_C0(f,tf);
		tf0 = THREAD_T0(tf);
		F_AREA(A,f,tf);
		
		/* power = power +  */
			/* A[0]*2*M_PI*((-F_P(f,tf) + 2*C_MU_L(c,tf)*C_DUDX(c,tf))*C_U(c,tf) +  */
					/* C_MU_L(c,tf)*(C_DVDX(c,tf) + C_DUDY(c,tf))*C_V(c,tf)) + */
			/* A[1]*2*M_PI*((-F_P(f,tf) + 2*C_MU_L(c,tf)*C_DVDY(c,tf))*C_V(c,tf) +  */
					/* C_MU_L(c,tf) * (C_DVDX(c,tf) + C_DUDY(c,tf))*C_U(c,tf)); */
		
		/* powerInviscid = powerInviscid + A[0]*2*M_PI*(-F_P(f,tf))*F_U(f,tf)  + A[1]*2*M_PI*(-F_P(f,tf))*F_V(f,tf);		// Inviscid power */
										
										
										
		/* _______________________________________________________  */
		
		/* thrustDebugVisc2 = thrustDebugVisc2 + 2*M_PI*(A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]); // exactly equivalent to Compute_Force_And_Moment in x-direction */
		power = power + -2*M_PI*( (A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]) * F_U(f,tf) 
								+ (A[1]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[1]) * F_V(f,tf) );
		/* _______________________________________________________  */

		
	} end_f_loop(f,tf)
	
	#endif
	
	return power;
	/* return powerInviscid; */
	/* return powerDebug; */
	/* return thrustDebugInv; */
	/* return thrustDebugVisc1; */
	/* return thrustDebugVisc2; */

}


/*---------- Compute_Power_Out -----------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
void Compute_Power_Out(double velx, double vely, 
						double *f_thrust_SMC_ptr, double *f_thrust_B_ptr, 
						double *f_drag_SMC_ptr, double *f_drag_B_ptr,
						double *P_thrust_SMC_ptr, double *P_thrust_B_ptr,
						double *P_drag_SMC_ptr, double *P_drag_B_ptr) 
{
	#if !RP_HOST	/*SERIAL or NODE */
	Thread *tf_ex = Lookup_Thread(domain, EX_ZONE);
	Thread *tf_sub = Lookup_Thread(domain, SUB_ZONE);
	Thread *tf_array[NUM_ZONES];
	Thread *tf, *tf0;
	face_t f;
	cell_t c;
	double A[ND_ND];
	double f_i, f_ip, f_if;
	double f_p = 0.0, f_f = 0.0;
	double f_thrust_SMC = 0.0;			/* thrust force [N] as defined by Sahin, Mohseni, and Colin */
	double f_thrust_B = 0.0;			/* thrust force [N] as defined by Borazjani et al */
	double f_drag_SMC = 0.0;			/* drag (friction) force [N] as defined by Sahin, Mohseni, and Colin */
	double f_drag_B = 0.0;				/* drag force [N] as defined by Borazjani et al */
	double P_thrust_SMC = 0.0;
	double P_thrust_B = 0.0;
	double P_drag_SMC = 0.0;
	double P_drag_B = 0.0;
	int i;

	tf_array[0] = tf_ex;
	tf_array[1] = tf_sub;
	
	/* loop through cells on  surface, calculating power from each  */
	for (i = 0; i < NUM_ZONES; i++)
	{
		tf = tf_array[i];
		begin_f_loop(f,tf) {  
			c = F_C0(f,tf);
			tf0 = THREAD_T0(tf);
			F_AREA(A,f,tf);
			
			/* power = power +  */
				/* A[0]*2*M_PI*((-F_P(f,tf) + 2*C_MU_L(c,tf)*C_DUDX(c,tf))*C_U(c,tf) +  */
						/* C_MU_L(c,tf)*(C_DVDX(c,tf) + C_DUDY(c,tf))*C_V(c,tf)) + */
				/* A[1]*2*M_PI*((-F_P(f,tf) + 2*C_MU_L(c,tf)*C_DVDY(c,tf))*C_V(c,tf) +  */
						/* C_MU_L(c,tf) * (C_DVDX(c,tf) + C_DUDY(c,tf))*C_U(c,tf)); */
			
			/* powerInviscid = powerInviscid + A[0]*2*M_PI*(-F_P(f,tf))*F_U(f,tf)  + A[1]*2*M_PI*(-F_P(f,tf))*F_V(f,tf);		// Inviscid power */
											
											
											
			/* _______________________________________________________ */
			/*		The following code has been verified			   */
			/* _______________________________________________________ */
			/* thrustDebugVisc2 = thrustDebugVisc2 + 2*M_PI*(A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]); // exactly equivalent to Compute_Force_And_Moment in x-direction */
			
			f_i = 2*M_PI*(A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]);
			f_ip = 2*M_PI*A[0]*F_P(f,tf);									/* pressure force */
			f_if = 2*M_PI* (-F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]);		/* friction force */
			
			f_p = f_p + f_ip;		/* sum up individual pressure contributions */
			f_f = f_f + f_if;		/* sum up individual friction contributions */

			/* f_i = -f_i; 	// Simulation makes the jellyfish swim in negative x-direction */
			/* f_thrust = f_thrust + (f_i + fabs(f_i) )/2; */
			/* f_drag = f_drag + (f_i - fabs(f_i) )/2; */
			/* power = power + ((f_i + fabs(f_i) )/2) * -F_U(f,tf);			// thrust power. Note that F_U has negative sign in front b/c jellyfish swims in the negative x-direction */
			/* _______________________________________________________  */
			
		} end_f_loop(f,tf)
	}
	
	
	/* if total pressure force helps the jellyfish swim then add it to thrust  */
	if (f_p < 0)
	{
		f_thrust_SMC = f_p;
		f_thrust_B = f_p;
	}
	/* only Borazjani's formulation includes pressure term in drag  */
	else
	{
		f_drag_B = f_p;
	}
	
	/* skin friction is always part of SMC's friction term (obviously)  */
	f_drag_SMC = f_f;
	/* if total skin friction retards the jellyfish then add it to drag  */
	if (f_f > 0)
	{
		f_drag_B = f_drag_B + f_f;
	}
	/* total skin friction helps the jellyfish swim forward (basically impossible, 
	but here for completeness)  */
	else
	{
		f_thrust_B = f_thrust_B + f_f;
	}

	P_thrust_SMC = f_thrust_SMC * velx;
	P_drag_SMC = f_drag_SMC * velx;
	P_thrust_B = f_thrust_B * velx;
	P_drag_B = f_drag_B * velx;
	
	*f_thrust_SMC_ptr = f_thrust_SMC;
	*f_thrust_B_ptr = f_thrust_B;
	*f_drag_SMC_ptr = f_drag_SMC;
	*f_drag_B_ptr = f_drag_B;
	
	*P_thrust_SMC_ptr = P_thrust_SMC;
	*P_thrust_B_ptr = P_thrust_B;
	*P_drag_SMC_ptr = P_drag_SMC;
	*P_drag_B_ptr = P_drag_B;
	
	/* return power_thrust; */
	/* return powerInviscid; */
	/* return powerDebug; */
	/* return thrustDebugInv; */
	/* return thrustDebugVisc1; */
	/* return thrustDebugVisc2; */
	
	#endif
}


/*---------- Compute_Averages -----------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
void Compute_Averages(double next_vel, double f_net, 
	double f_thrust_SMC, double f_thrust_B, double f_drag_SMC, 
	double f_drag_B, double P_in, double P_thrust_SMC, 
	double P_thrust_B, double P_drag_SMC, double P_drag_B) 
{
	#if !RP_NODE	/*SERIAL or HOST */
	double vel_avg, Re_avg;
	double f_avg, f_thrust_SMC_avg, f_thrust_B_avg, f_drag_SMC_avg, f_drag_B_avg;
	double P_in_avg, P_thrust_SMC_avg, P_thrust_B_avg, P_drag_SMC_avg, P_drag_B_avg;
	double eff_SMC, eff_B, COT;
	double mass_jelly = Get_Mass();
	int cycleN = Get_CycleNumber();
	FILE *favgs;
	
	
	/* add new values to running totals  */
	vel_tot = vel_tot + next_vel;
	
	f_tot = f_tot + f_net;
	f_thrust_SMC_tot = f_thrust_SMC_tot + f_thrust_SMC;
	f_thrust_B_tot = f_thrust_B_tot + f_thrust_B;
	f_drag_SMC_tot = f_drag_SMC_tot + f_drag_SMC;
	f_drag_B_tot = f_drag_B_tot + f_drag_B;
	
	P_in_tot = P_in_tot + P_in;
	
	P_thrust_SMC_tot = P_thrust_SMC_tot + P_thrust_SMC;
	P_thrust_B_tot = P_thrust_B_tot + P_thrust_B;
	P_drag_SMC_tot = P_drag_SMC_tot + P_drag_SMC;
	P_drag_B_tot = P_drag_B_tot + P_drag_B;
	
	kavg = kavg + 1;
	
	/* at the end of each swimming cycle, compute cycle average values and print out to file  */
	Message("New_Cycle = %i\n", NewCycle());
	if (NewCycle()) {
		Message("Computing averages... ");
		
		vel_avg = vel_tot/kavg;	/* go back to using the above code after confirming the denominator  */
		Re_avg = (vel_avg * DIAMETER)/NU;
		f_avg = f_tot/kavg;
		f_thrust_SMC_avg = f_thrust_SMC_tot/kavg;
		f_thrust_B_avg = f_thrust_B_tot/kavg;
		f_drag_SMC_avg = f_drag_SMC_tot/kavg;
		f_drag_B_avg = f_drag_B_tot/kavg;
		P_in_avg = P_in_tot/kavg;
		P_thrust_SMC_avg = P_thrust_SMC_tot/kavg;
		P_thrust_B_avg = P_thrust_B_tot/kavg;
		P_drag_SMC_avg = P_drag_SMC_tot/kavg;
		P_drag_B_avg = P_drag_B_tot/kavg;
		
		eff_SMC = (vel_avg*f_thrust_SMC_avg)/P_in_avg;	/* propulsive efficiency w/ Sahin, Mohseni, Colin thrust definition */
		eff_B = (vel_avg*f_thrust_B_avg)/P_in_avg;		/* propulsive efficiency w/ Borazjani thrust definition */
		
		COT = P_in_avg/(mass_jelly*-vel_avg);			/* cost of transport [W/(kg*m/s)] */
		
		Message("Writing to averages_and_eff.txt... ");
		/* write out averages and efficiencies  */
		if (cycleN == 2) {		/* create file and write header */
			favgs = fopen("averages_and_eff.txt","w+");
			fprintf(favgs,"cycleN\t");
			fprintf(favgs,"Time (s)\t");
			fprintf(favgs,"Vel_avg (m/s)\t");
			fprintf(favgs,"Re_body\t");
			fprintf(favgs,"Force_avg (N)\t");
			fprintf(favgs,"Thrust_SMC_avg (N)\t");
			fprintf(favgs,"Thrust_B_avg (N)\t");
			fprintf(favgs,"Drag_SMC_avg (N)\t");
			fprintf(favgs,"Drag_B_avg (N)\t");
			fprintf(favgs,"Power_in_avg (W)\t");
			fprintf(favgs,"Power_thrust_SMC_avg (W)\t");
			fprintf(favgs,"Power_thrust_B_avg (W)\t");
			fprintf(favgs,"Power_drag_SMC_avg (W)\t");
			fprintf(favgs,"Power_drag_B_avg (W)\t");
			fprintf(favgs,"Efficiency_SMC\t");
			fprintf(favgs,"Efficiency_B\t");
			fprintf(favgs,"COT (J/kg/m)\t");
			fprintf(favgs,"kavg\n");
			fclose(favgs);
		}
		
		favgs = fopen("averages_and_eff.txt","a+");
		fprintf(favgs,"%i\t", cycleN-1);
		fprintf(favgs,"%f\t", CURRENT_TIME);
		fprintf(favgs,"%16.12f\t", vel_avg);
		fprintf(favgs,"%16.12f\t", Re_avg);
		fprintf(favgs,"%16.12f\t", f_avg);
		fprintf(favgs,"%16.12f\t", f_thrust_SMC_avg);
		fprintf(favgs,"%16.12f\t", f_thrust_B_avg);
		fprintf(favgs,"%16.12f\t", f_drag_SMC_avg);
		fprintf(favgs,"%16.12f\t", f_drag_B_avg);
		fprintf(favgs,"%16.12f\t", P_in_avg);
		fprintf(favgs,"%16.12f\t", P_thrust_SMC_avg);
		fprintf(favgs,"%16.12f\t", P_thrust_B_avg);
		fprintf(favgs,"%16.12f\t", P_drag_SMC_avg);
		fprintf(favgs,"%16.12f\t", P_drag_B_avg);
		fprintf(favgs,"%16.12f\t", eff_SMC);
		fprintf(favgs,"%16.12f\t", eff_B);
		fprintf(favgs,"%f\t", COT);
		fprintf(favgs,"%i\n", kavg);
		fclose(favgs);
		
		/* reset running totals and count  */
		vel_tot = 0.0;
		f_tot = 0.0;
		f_thrust_SMC_tot = 0.0;
		f_thrust_B_tot = 0.0;
		f_drag_SMC_tot = 0.0;
		f_drag_B_tot = 0.0;
		P_in_tot = 0.0;
		P_thrust_SMC_tot = 0.0;
		P_thrust_B_tot = 0.0;
		P_drag_SMC_tot = 0.0;
		P_drag_B_tot = 0.0;
		kavg = 0;
		
		Message("Done.\n");
	}
	#endif
}


	
/*---------- Create_Disp_OutputFile -------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
void Create_Disp_OutputFile(char side)  
{
	#if !RP_NODE	/*SERIAL or HOST */
	FILE *fdisp;
	
	if (side == 'e') {
		fdisp = fopen("displacement.txt","w+");
		fprintf(fdisp,"Time (s)\t");
		fprintf(fdisp,"Cycle #\t");
		fprintf(fdisp,"x (m)\t");
		fprintf(fdisp,"del x (m)\t");
		fprintf(fdisp,"V_t (m/s)\t");
		fprintf(fdisp,"Net force (N)\t");
		fprintf(fdisp,"Thrust_SMC (N)\t");
		fprintf(fdisp,"Thrust_B (N)\t");
		fprintf(fdisp,"Drag_SMC (N)\t");
		fprintf(fdisp,"Drag_B (N)\t");			
		fprintf(fdisp,"Power in (W)\t");
		fprintf(fdisp,"Power_thrust_SMC (W)\t");
		fprintf(fdisp,"Power_thrust_B (W)\t");
		fprintf(fdisp,"Power_drag_SMC (W)\t");
		fprintf(fdisp,"Power_drag_B (W)\t");
		fprintf(fdisp,"Paused\n");
		
		fclose(fdisp);
	}
	else if (side == 's') {
		fdisp = fopen("displacement_sub.txt","w+");
		fprintf(fdisp,"Time (s)\t");
		fprintf(fdisp,"Cycle #\t");
		fprintf(fdisp,"x (m)\t");
		fprintf(fdisp,"del x (m)\t");
		fprintf(fdisp,"V_t (m/s)\t");
		fprintf(fdisp,"Net force (N)\t");
		fprintf(fdisp,"Thrust_SMC (N)\t");
		fprintf(fdisp,"Thrust_B (N)\t");
		fprintf(fdisp,"Drag_SMC (N)\t");
		fprintf(fdisp,"Drag_B (N)\t");			
		fprintf(fdisp,"Power in (W)\t");
		fprintf(fdisp,"Power_thrust_SMC (W)\t");
		fprintf(fdisp,"Power_thrust_B (W)\t");
		fprintf(fdisp,"Power_drag_SMC (W)\t");
		fprintf(fdisp,"Power_drag_B (W)\t");
		fprintf(fdisp,"Paused\n");
		fclose(fdisp);
	}
	else Error("\n Incorrect side identifier used in create_Disp_OutputFile() \n");
	
	#endif
}


/*---------- write_to_Disp_OutputFile -------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
void Write_to_Disp_OutputFile(char side, double x, double vel, double del_x, 
	double f_net, double f_thrust_SMC, double f_thrust_B, double f_drag_SMC, double f_drag_B, 
	double P_in, double P_thrust_SMC, double P_thrust_B, double P_drag_SMC, double P_drag_B)  
{
	#if !RP_NODE	/*SERIAL or HOST */
	FILE *fdisp;
	int cycleN = Get_CycleNumber();
	
	if (side == 'e') {
		fdisp = fopen("displacement.txt","a+");
		fprintf(fdisp,"%f\t",CURRENT_TIME);
		fprintf(fdisp,"%i\t",cycleN);
		fprintf(fdisp,"%16.12f\t",x);
		fprintf(fdisp,"%16.12f\t",del_x);
		fprintf(fdisp,"%16.12f\t",vel);
		fprintf(fdisp,"%16.12f\t",f_net);
		fprintf(fdisp,"%16.12f\t",f_thrust_SMC);
		fprintf(fdisp,"%16.12f\t",f_thrust_B);
		fprintf(fdisp,"%16.12f\t",f_drag_SMC);
		fprintf(fdisp,"%16.12f\t",f_drag_B);
		fprintf(fdisp,"%16.12f\t",P_in);
		fprintf(fdisp,"%16.12f\t",P_thrust_SMC);
		fprintf(fdisp,"%16.12f\t",P_thrust_B);
		fprintf(fdisp,"%16.12f\t",P_drag_SMC);
		fprintf(fdisp,"%16.12f\t",P_drag_B);
		fprintf(fdisp,"%i\n",isPaused);
		fclose(fdisp);
	}
	else if (side == 's') {
		fdisp = fopen("displacement_sub.txt","a+");
		fprintf(fdisp,"%f\t",CURRENT_TIME);
		fprintf(fdisp,"%i\t",cycleN);
		fprintf(fdisp,"%16.12f\t",x);
		fprintf(fdisp,"%16.12f\t",del_x);
		fprintf(fdisp,"%16.12f\t",vel);
		fprintf(fdisp,"%16.12f\t",f_net);
		fprintf(fdisp,"%16.12f\t",f_thrust_SMC);
		fprintf(fdisp,"%16.12f\t",f_thrust_B);
		fprintf(fdisp,"%16.12f\t",f_drag_SMC);
		fprintf(fdisp,"%16.12f\t",f_drag_B);
		fprintf(fdisp,"%16.12f\t",P_in);
		fprintf(fdisp,"%16.12f\t",P_thrust_SMC);
		fprintf(fdisp,"%16.12f\t",P_thrust_B);
		fprintf(fdisp,"%16.12f\t",P_drag_SMC);
		fprintf(fdisp,"%16.12f\t",P_drag_B);
		fprintf(fdisp,"%i\n",isPaused);
		fclose(fdisp);
	}
	else Error("\n Incorrect side identifier used in write_to_Disp_OutputFile() \n");
	
	#endif
}

