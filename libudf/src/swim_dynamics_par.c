#include "swim_dynamics_par.h"

/** global variables * */
static double x_prev = 0.0;
static double x_target = 0.0;
static double del_x_prev = 0.0;
static double preVel = 0.0;		        /* x-velocity one iteration ago (needed for finite diff) */
static double fx_prev = 0.0;

static double xS_prev = 0.0;
static double xS_target = 0.0;
static double del_xS_prev = 0.0;
static double preVelS = 0.0;
static double fxS_prev = 0.0;

int ke = 0;
int ks = 0;

const double ALPHA_F = 0.05;		/* Settings for explicit scheme */
const double ALPHA_X = 1.0;			/* Settings for explicit scheme */

/* const double ALPHA_F = 1.0;		// Settings for implicit scheme */
/* const double ALPHA_X = 0.4;		// Settings for implicit scheme */


/** Function definitions * */

/*---------- GetForce_x -------------------------------------------
Purpose: 

Input:	

Output: Returns force in the x-direction
----------------------------------------------------------------------- */
double Get_Force_x(Thread *tf) 
{	
	double force_x = 0.0;
	
	#if !RP_HOST
		double x_cg[ND_ND], force[ND_ND], m_glob[ND_ND];
		int i;
		
		for (i = 0; i<ND_ND; i++)
			force[i] = 0.0;
		
		#if !PARALLEL /* if serial */
			Compute_Force_And_Moment(domain, tf, x_cg, force, m_glob, TRUE);	
		#endif
		
		#if RP_NODE
			Compute_Force_And_Moment(domain, tf, x_cg, force, m_glob, FALSE);
			Message("FORCE = %f (I am node %i)\n", force[0], myid);
		#endif
		
		force_x = force[0];
	#endif
		
	return force_x;
}




/*---------- Get_Force2_x -------------------------------------------
Purpose: 

Input:	

Output: Returns force in the x-direction
----------------------------------------------------------------------- */
double Get_Force2_x(Thread *tf) 
{	
	face_t f;
	cell_t c;
	double A[ND_ND];
	double f_i, f_ip, f_if;
	double f_p = 0.0, f_f = 0.0, fx = 0;
	int i;
	
	/* loop through cells on  surface, calculating power from each  */
	begin_f_loop(f,tf) 
	{  
		F_AREA(A,f,tf);
					
		/* _______________________________________________________ */
		/*		The following code has been verified			   */
		/* _______________________________________________________ */
		/* thrustDebugVisc2 = thrustDebugVisc2 + 2*M_PI*(A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]); // exactly equivalent to Compute_Force_And_Moment in x-direction */
		
		f_i = 2*M_PI*(A[0]*F_P(f,tf) - F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]);
		f_ip = 2*M_PI*A[0]*F_P(f,tf);									/* pressure force */
		f_if = 2*M_PI* (-F_STORAGE_R_N3V(f,tf,SV_WALL_SHEAR)[0]);		/* friction force */
		
		f_p = f_p + f_ip;		/* sum up individual pressure contributions */
		f_f = f_f + f_if;		/* sum up individual friction contributions */
		fx = fx + f_i;
		/* _______________________________________________________  */
		
	} end_f_loop(f,tf)

	return fx;
}
	


/*---------- GetBodyVelocity -------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
double Get_BodyVelocity(Thread *tf, char meshZone, double fx_tot, double *del_x_ptr)
{
	double vel = 0.0;
	
	#if !RP_HOST
		double mass_jelly = Get_Mass();
		double fx, fxS, fx_raw, fxS_raw; 
		double nextVel = 0.0, nextVelS = 0.0;
		double x, xS, x_raw, xS_raw;
		double x_position;
		double del_x = 0.0, del_xS = 0.0;
		
		
		x_position = Rezero(tf);
		
		if ((myid == 0) || (myid == 1000000))	
		{
			
			/* Message("\tx_position = %16.12f, ", x_position); */
			/* Message("x_prev = %16.12f, x_prevS = %16.12f\n", x_prev, xS_prev); */
			
			
			if ((meshZone == 'e') || (meshZone == 'i')) 
			{
				if (NewTimeStepForThisZone('e') && NewTimeStepForThisZone('i')) 
				{
					Message("\t-------- calc -------- \n");
					ke = 1;
					
					if (N_TIME == 0) 	/* first time step */
					{
						x = x_position;
						x_prev = x;
						fx = fx_tot;
						fx_prev = fx_tot;
						del_x = 0.0;
						del_x_prev = 0.0;
						nextVel = 0.0;
						preVel = 0.0;
					}
					else				/* all other time steps */
					{
						/* Message("mass = %lf\t", mass_jelly); */
				
						/* Calc new velocity and x-displacement  */		
						fx_raw = fx_tot;
						fx = ALPHA_F*fx_raw + (1-ALPHA_F)*fx_prev;
						/* Message("fx = %16.12f\n",fx); */
						nextVel = preVel + (fx/mass_jelly) * CURRENT_TIMESTEP; 		/* forward diff, 1st order accurate */
						/* Message("nextVel = %lf, preVel = %lf, (fx/mass_jelly) = %lf, CURRENT_TIMESTEP = %lf\n", nextVel, preVel, fx/mass, CURRENT_TIMESTEP); */
						x_raw = x_prev + nextVel * CURRENT_TIMESTEP;
						/* Message("x_raw = %lf, x_prev = %lf, nextVel*CURRENT_TIMESTEP = %lf\n", x_raw, x_prev, nextVel*CURRENT_TIMESTEP); */
						x = ALPHA_X*x_raw + (1-ALPHA_X)*x_prev;
						/* Message("x = %lf, ALPHA_X*x_raw = %lf, (1-ALPHA_X)*x_prev = %lf\n", x, ALPHA_X*x_raw, (1-ALPHA_X)*x_prev); */
					}

					del_x = x - x_prev;		/* this is what is used by moveJellyNodes() */
					/* Message("del_x = %lf\n", del_x); */
					
					/* store values for next time step/iteration  */
					x_prev = x;
					del_x_prev = del_x;
					preVel = nextVel;
					fx_prev = fx;
				}
				else 
				{
					ke = ke + 1;
					x = x_prev;
					del_x = del_x_prev;
					nextVel = preVel;
				}
				
				/* Message("Iteration %i through EXUM: del_x = %16.12f, vel_x = %16.12f\n", ke, del_x, preVel); */
				vel = nextVel; 
				*del_x_ptr = del_x;
			}
			
			else if (meshZone == 's') {
				if (NewTimeStepForThisZone(meshZone)) 
				{
					Message("\t-------- calc -------- \n");
					ks = 1;
					
					if (N_TIME == 0) 	/* first time step */
					{
						xS = x_position;
						xS_prev = xS;
						fxS = fx_tot;
						fxS_prev = fx_tot;
						del_xS = 0.0;
						del_xS_prev = 0.0;
						nextVelS = 0.0;
						preVelS = 0.0;
					}
					else				/* all other time steps */
					{
						/* Message("mass = %lf\t", mass_jelly); */
						/* Calc new velocity and x-displacement  */
						fxS_raw = fx_tot;
						fxS = ALPHA_F*fxS_raw + (1-ALPHA_F)*fxS_prev;
						/* Message("fx = %16.12f\n",fxS); */
						nextVelS = preVelS + (fxS/mass_jelly) * CURRENT_TIMESTEP; 		/* forward diff, 1st order accurate */
						xS_raw = xS_prev + nextVelS * CURRENT_TIMESTEP;
						xS = ALPHA_X*xS_raw + (1-ALPHA_X)*xS_prev;	
					}
					
					del_xS = xS - xS_prev;	
					/* Message("del_xS = %lf\n", del_xS); */
					
					/* store values for next time step/iteration  */
					xS_prev = xS;
					preVelS = nextVelS;	
					del_xS_prev = del_xS;
					fxS_prev = fxS;
				}
				else 
				{
					ks = ks + 1;
					xS = xS_prev;
					nextVelS = preVelS;
					del_xS = del_xS_prev;
				}
				/* Message("Iteration %i through SUB: del_x = %16.12f, vel_x = %16.12f\n", ks, del_xS, preVelS); */
				
				vel = nextVelS;
				*del_x_ptr = del_xS;
			}
		}
	#endif
	
	/* (Parallel) Sync velocities across all nodes  */
	vel = PRF_GRSUM1(vel);
	*del_x_ptr = PRF_GRSUM1(*del_x_ptr);
	
	
	return vel;
}



/*---------- Get_Mass -------------------------------------------
Purpose: Calculate the mass of the jellyfish according to its volume,
		assuming the animal is neutrally buoyant (i.e. density of 
		animal = density of the surrounding fluid)

Input:	None
Output: Mass in [kg]
---------------------------------------------------------------- */
double Get_Mass() 
{	
	double density = 998.2;				/* [kg/m^3] fluent default density for liquid water */
	double volume = 0.0001376821611;		/* 0.0001376821611 [m^3] Measured from initial geometry (t=0) w/ SpaceClaim */
	double mass = density*volume;
	/* mass = 0.1675; 			   // [kg] from m[g]=0.29*(D[cm])^2.32 (Omori, Ishii, and Fujinaga, 1995) for 15.5cm  */

	return mass;
}

/*---------- Get_Vel_At_End -------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
double Get_Vel_At_End(char meshZone)
{
	#if !RP_HOST
		if (meshZone == 'e')  return preVel;
		else if (meshZone == 's') return preVelS;
		else 
		{
			Error("Wrong meshZone char specified");
			return 0.0;
		}
	#endif
	
	#if RP_HOST
		return 0.0;
	#endif
}

/*---------- Get_Del_At_End -------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
double Get_Del_At_End(char meshZone)
{
	#if !RP_HOST
		if (meshZone == 'e') return del_x_prev;
		else if (meshZone == 's') return del_xS_prev;
		else
		{
			Error("Wrong meshZone char specified");
			return 0.0;
		}
	#endif
	
	#if RP_HOST
		return 0.0;
	#endif
}


/*---------- Get_BodyVelocity_Implicit -------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
double Get_BodyVelocity_Implicit(Thread *tf, char meshZone, double fx_tot, double *del_x_ptr)
{
	double vel = 0.0;
	
	#if !RP_HOST
		double mass_jelly = Get_Mass();
		double nextVel = 0.0, nextVelS = 0.0;
		double x, x_raw, xS_raw;
		double x_position;
		double del_x = 0.0, del_xS = 0.0;
		int endFlag = 0;
		
		
		x_position = Rezero(tf);
		
		if ((myid == 0) || (myid == 1000000))	
		{
			
			/* Message("\tx_position = %16.12f, ", x_position); */
			/* Message("x_prev = %16.12f, x_prevS = %16.12f\n", x_prev, xS_prev); */
			
			
			if (meshZone == 'e') 
			{
				if (N_TIME == 0) 
				{
					/* Message("N_TIME = %i, preTime_ex = %i\n", N_TIME, preTime_ex); */
					if (NewTimeStepForThisZone(meshZone)) 
					{
						ke = 1;
						x = x_position;
						x_prev = x;
						fx_prev = fx_tot;
						del_x = 0.0;
						del_x_prev = 0.0;
						nextVel = 0.0;
						preVel = 0.0;
					}
					else {
						/* do nothing  */
					}
				}
				else /* N_TIME > 0 */
				{
					if (NewTimeStepForThisZone(meshZone)) 
					{
						ke = 1;
						
						/* Calc INITIAL new velocity and x-displacement  */		
						nextVel = preVel + (fx_prev/mass_jelly) * CURRENT_TIMESTEP; 		/* forward diff w/o smoothing, 1st order accurate */
						x_target = x_prev + nextVel*CURRENT_TIMESTEP;
					}
					else if ((ke%2 == 1) && (endFlag == 0))
					{
						/* Calc new velocity and x-displacement  */		
						Message("Recalculating velocity... ");
						nextVel = preVel + (fx_tot/mass_jelly) * CURRENT_TIMESTEP; 		/* backwards diff, 1st order accurate */
						x_raw = x_prev + nextVel*CURRENT_TIMESTEP;
						x_target = ALPHA_X*x_raw + (1-ALPHA_X)*x_target;		/* <--- THIS IS THE EQ THAT IS DIFFERENT FOR IMPLICIT */
						Message("Done. \n");
					}
					else {
						/* Do nothing when UDF is being called by Laplace mesh smoothing  */
					}
				}
				
				/* calculate delta_x at beginning of time step and every other
					subsequent iteration  */
				if ((ke%2 == 1) && (endFlag == 0)) 
				{
					del_x = x_target - x_prev;
					vel = nextVel;
					*del_x_ptr = del_x;
				}
				
				/* Message("Iteration %i through EXUM: del_x = %16.12f, vel_x = %16.12f\n", ke, del_x, preVel); */
				Message("ke = %i\n", ke);
				ke = ke + 1;
			}
			
			if (meshZone == 's')
			{
				if (N_TIME == 0) 
				{
					if (NewTimeStepForThisZone(meshZone)) 
					{
						ks = 1;
						xS_target = x_position;
						xS_prev = x_position;
						fxS_prev = 0.0;
						del_xS = 0.0;
						nextVelS = 0.0;
						preVelS = 0.0;
					}
					else {
						/* do nothing  */
					}
				}
				else /* N_TIME > 0 */
				{
					if(NewTimeStepForThisZone(meshZone)) 
					{
						ks = 1;
						
						/* Calc new velocity and x-displacement  */
						nextVelS = preVelS + (fxS_prev/mass_jelly) * CURRENT_TIMESTEP; 		/* forward diff w/o smoothing, 1st order accurate */
						xS_target = xS_prev + nextVelS*CURRENT_TIMESTEP;
						
					}
					else if ((ks%2 == 1) && (endFlag == 0)) {
						/* Calc new velocity and x-displacement  */
						Message("Recalculating velocity... ");
						nextVelS = preVelS + (fx_tot/mass_jelly) * CURRENT_TIMESTEP; 		/* backwards diff, 1st order accurate */
						xS_raw = xS_prev + nextVelS*CURRENT_TIMESTEP;
						xS_target = ALPHA_X*xS_raw + (1-ALPHA_X)*xS_target;			/* <--- THIS IS THE EQ THAT IS DIFFERENT FOR IMPLICIT */
						Message("Done. \n");
					}
					else {
						/* Do nothing when UDF is being called by Laplace mesh smoothing  */
					}
				}
					
				if ((ks%2 == 1) && (endFlag == 0))
				{
					del_xS = xS_target - xS_prev;
					vel = nextVelS;
					*del_x_ptr = del_xS;
				}
					
				Message("ks = %i\n", ks);
				ks = ks + 1;
			}
		}
	#endif
	
	/* (Parallel) Sync velocities across all nodes  */
	vel = PRF_GRSUM1(vel);
	*del_x_ptr = PRF_GRSUM1(*del_x_ptr);
	
	
	return vel;
}
