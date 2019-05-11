#include "timing_and_counters_par.h"

/** Global var **/
int isPaused = 0;
int isRamping = 0;

/* const int N_RAMPSTEPS = 10;							// number of iterations to ramp down/up into/from pause*/

static double tSim_ex = 0.0;
static double tSim_sub = 0.0;
static double tSim_intfc = 0.0;

static int preIt_ex = -1;
static int preIt_sub = -1;
static int preIt_intfc = -1;

static int preTime_ex = -1;
static int preTime_sub = -1;
static int preTime_intfc = -1;

static int cycleN_prev = 1;


/** Function definitions **/
/*---------- Print_Timestep_Header ------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
void Print_Timestep_Header()
{
	if (NewTimeStep()) Message("\n\n**********\nTime step: %i\n**********",N_TIME);
}


/*---------- Print_Zone_Header ------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
void Print_Zone_Header(char meshZone)
{
	if (meshZone == 'e') 
	{
		#if !PARALLEL
			Message("---EXUMBRELLA (time step %i)---\n", N_TIME);
		#endif
			
		#if RP_NODE
			/* Message0("---EXUMBRELLA (node%i, time step %i)---\n", myid, N_TIME);*/
		#endif
		
		#if RP_HOST
			Message("---EXUMBRELLA (HOST, time step %i)---\n", N_TIME);
		#endif
	}
	else if (meshZone == 's') 
	{
		#if !PARALLEL
			Message("---SUBUMBRELLA (time step %i)---\n", N_TIME);
		#endif
			
		#if RP_NODE
			/* Message0("---SUBUMBRELLA (node%i, time step %i)---\n", myid, N_TIME);*/
		#endif
		
		#if RP_HOST
			Message("---SUBUMBRELLA (HOST, time step %i)---\n", N_TIME);
		#endif
	}
}


/*---------- Get_EffectivePeriod() --------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
double Get_EffectivePeriod()
{
	double period_eff = PERIOD*(1 + PAUSE_FRAC);
	return period_eff;
}


/*---------- Get_N_RAMPSTEPS() --------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int Get_N_RAMPSTEPS(double ramp_frac)
{
	int N_RAMPSTEPS = 0;
	
	if (ramp_frac > PAUSE_FRAC)
		Error("ramp duration exceeds the time the jellyfish is paused!");
	else
		N_RAMPSTEPS = floor( (ramp_frac*PERIOD)/CURRENT_TIMESTEP );
	
	return N_RAMPSTEPS;
}


/*---------- Smootherstep() --------------------------------------
Purpose: transition from edge0 to edge1 values with 1st and 2nd
derivatives = zero. 

From AMD/Wikipedia
----------------------------------------------------------------*/
float Smootherstep(float edge0, float edge1, float x) 
{
  /* Scale, and clamp x to 0..1 range*/
  x = Clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  /* Evaluate polynomial*/
  return x * x * x * (x * (x * 6 - 15) + 10);
}


/*---------- Clamp() --------------------------------------
Purpose: transition from edge0 to edge1 values with 1st and 2nd
derivatives = zero. 

From AMD/Wikipedia
----------------------------------------------------------------*/
float Clamp(float x, float lowerlimit, float upperlimit) 
{
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
}


/*---------- Get_SimTimeInCycle() --------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
double Get_SimTimeInCycle(char meshZone)
{
	double period_eff = Get_EffectivePeriod();
	double tau = fmod(CURRENT_TIME,period_eff);
	double gamma;
	double tauSim = 0.0;
	int cycleN = Get_CycleNumber();
	double timestepSim;
	int N_RAMPSTEPS = 0;
	
	/* ---------------------------------------------------------------------------------------  */
	/* Normal simulation 																		*/
	/* ---------------------------------------------------------------------------------------  */
	if (PAUSE_FRAC <= 0.0) 
	{
		isPaused = 0;
		isRamping = 0;
		gamma = 1.0;
		tauSim = tau;
	}
	
	/* ---------------------------------------------------------------------------------------  */
	/* Pause simulation 																		*/
	/* ---------------------------------------------------------------------------------------  */
	/* Create a simulated time vector that will "slow down" and "speed up" in order to slow the */
	/* kinematics to a stop when ramping down to pause or bring the jellyfish back to motion by */
	/* ramping up from pause.  																	*/
	/* ---------------------------------------------------------------------------------------- */
	else
	{
		N_RAMPSTEPS = Get_N_RAMPSTEPS(0.24);	/* ramp lasts 5% of period*/
		
		/* first time step */
		if (N_TIME == 0) 
		{ 		
			isPaused = 0;
			isRamping = 0;
			gamma = 1.0;
			tauSim = tau;
		}
		
		/* all other time steps */
		else
		{
			/* time is PAUSED or SLOWING DOWN (ramping down) into pause */
			/* if ((tau > PERIOD) || (tau == 0)) 		// Pause lasts from (period, period_eff] */
			if ((tau >= PERIOD - 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS) &&
				(tau < period_eff - 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS)) 	
			{
				isPaused = 1; 		
				
				/* SLOWING DOWN (ramping down) */
				/* if ((tau <= PERIOD + CURRENT_TIMESTEP*N_RAMPSTEPS) && (tau != 0.0)) */
				if (tau <= PERIOD + 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS)
				{		
					isRamping = 1;
					/* gamma = 1 - (tau - PERIOD)/(CURRENT_TIMESTEP*N_RAMPSTEPS);  // linear ramp */
					/* gamma = 0.5*(cos((tau-PERIOD)*M_PI/(CURRENT_TIMESTEP*N_RAMPSTEPS))+1); // cosine ramp*/
					gamma = 1 - Smootherstep(PERIOD - 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS, 
									PERIOD + 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS, tau);
				}
				
				/* FULLY PAUSED */
				else 
				{
					isRamping = 0;
					gamma = 0;
				}
			}
			
			/* time is NORMAL or SPEEDING UP (ramping up) from pause */
			else
			{
				isPaused = 0;		
				
				/* SPEEDING UP (ramping up) */
				/* if ((cycleN > 1) && (tau < CURRENT_TIMESTEP*N_RAMPSTEPS)) */
				if (tau >= period_eff - 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS)
				{
					isRamping = 1;
					/* gamma = tau/(CURRENT_TIMESTEP*N_RAMPSTEPS);	// linear ramp*/
					/* gamma = 0.5*(-cos(tau*M_PI/(CURRENT_TIMESTEP*N_RAMPSTEPS))+1);	// cosine ramp*/
					gamma = Smootherstep(period_eff - 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS,
							period_eff + 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS, tau); 
				}
				else if ((cycleN > 1) && (tau < 0.5*CURRENT_TIMESTEP*N_RAMPSTEPS)) 
				{
					isRamping = 1;
					/* gamma = tau/(CURRENT_TIMESTEP*N_RAMPSTEPS);	// linear ramp*/
					/* gamma = 0.5*(-cos(tau*M_PI/(CURRENT_TIMESTEP*N_RAMPSTEPS))+1);	// cosine ramp*/
					gamma = Smootherstep(-0.5*CURRENT_TIMESTEP*N_RAMPSTEPS,
								0.5*CURRENT_TIMESTEP*N_RAMPSTEPS, tau); 
				}
				
				/* NORMAL TIME */
				else 
				{
					isRamping = 0;
					gamma = 1;
				}
			}
		}
	}

	/* print the results to console (both normal and paused) */
	timestepSim = gamma*CURRENT_TIMESTEP;
	
	if (meshZone == 'e') {
		if (NewTimeStepForThisZone('e') && (N_TIME > 0)) {
			tSim_ex = tSim_ex + timestepSim;		
		}
		tauSim = fmod(tSim_ex,PERIOD);
		
		if ( NewTimeStepForThisZone(meshZone) )
		{
			Message0("Current time: %lf\ttSim_ex: %lf\ttauSim: %lf\ttimestepSim: %lf\n",
				CURRENT_TIME, tSim_ex, tauSim, timestepSim);
		}
	}
	else if (meshZone == 's') {
		if (NewTimeStepForThisZone('s') && (N_TIME > 0)) {
			tSim_sub = tSim_sub + timestepSim;		
		}
		tauSim = fmod(tSim_sub,PERIOD);
		
		if ( NewTimeStepForThisZone(meshZone) )
		{
			Message0("Current time: %lf\ttSim_sub: %lf\ttauSim: %lf\ttimestepSim: %lf\n",
				CURRENT_TIME, tSim_sub, tauSim, timestepSim);
		}
	}
	else if (meshZone == 'i') {
		if (NewTimeStepForThisZone('i') && (N_TIME > 0)) {
			tSim_intfc = tSim_intfc + timestepSim;		
		}
		tauSim = fmod(tSim_intfc,PERIOD);
	}
	
	
	return tauSim;
}

/*---------- NewTimeStepForThisZone ---------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int NewTimeStepForThisZone(char zone)
{
	int new_zone_ts = 0;

	switch (zone)
	{
		case 'e' :
			if (N_TIME != preTime_ex) 
				new_zone_ts = 1;
			break;
		
		case 's' :
			if (N_TIME != preTime_sub) 
				new_zone_ts = 1;
			break;
			
		case 'i' :
			if (N_TIME != preTime_intfc) 
				new_zone_ts = 1;
			break;
		
		default :
			Message0("(timing_and_counters_par.c::NewTimeStepForThisZone) Invalid meshzone specified \n");
	}
	
	return new_zone_ts;
}

/*---------- NewIterationForThisZone ---------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int NewIterationForThisZone(char zone)
{
	int new_zone_it = 0;
	
	switch (zone)
	{
		case 'e' :
			if (N_ITER != preIt_ex) 
				new_zone_it = 1;
			break;
		
		case 's' :
			if (N_ITER != preIt_sub) 
				new_zone_it = 1;
			break;
			
		case 'i' :
			if (N_ITER != preIt_intfc) 
				new_zone_it = 1;
			break;
		
		default :
			Message0("(timing_and_counters_par.c::NewIterationForThisZone) Invalid meshzone specified \n");
	}
	
	return new_zone_it;
}

/*---------- NewIteration ---------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int NewIteration()
{
	int new_it = 0;
	
	if ((N_ITER != preIt_ex) && (N_ITER != preIt_sub) && 
		(N_ITER != preIt_intfc) )
	{
		new_it = 1;
	}

	return new_it;
}


/*---------- NewTimeStep ---------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int NewTimeStep()
{
	int new_ts;
	
	if ((N_TIME != preTime_ex) && (N_TIME != preTime_sub) &&
	    (N_TIME != preTime_intfc) )
	{
		new_ts = 1;
	}

	return new_ts;
}


/*---------- Get_CycleNumber -------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
int Get_CycleNumber() 
{
	double period_eff = Get_EffectivePeriod();			/* function found in timing_jelly.c*/
	int cycleN = (int)(CURRENT_TIME/period_eff) + 1;
	
	return cycleN;
}

/*---------- Store_OldCycleN -------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
void Store_OldCycleN()
{
	int cycleN = Get_CycleNumber();
	
	/* Message("cycleN_prev = %i, cycleN = %i --> new_Cycle = %i\n", cycleN_prev, cycleN, NewCycle());*/
	cycleN_prev = cycleN;
	
}

/*---------- NewCycle -------------------------------------------
Purpose: 

Input:	
----------------------------------------------------------------*/
int NewCycle() 
{
	int cycleN = Get_CycleNumber();
	
	if(cycleN != cycleN_prev)
		return 1;	
	else 
		return 0;
}


/*---------- UpdateCounters -------------------------------------
Purpose: 

Input:	

Output:
----------------------------------------------------------------*/
void UpdateCounters(char zone)
{
	switch (zone)
	{
		case 'e' :
			if (N_TIME != preTime_ex) preTime_ex = N_TIME;
			if (N_ITER != preIt_ex) preIt_ex = N_ITER;
			break;
		
		case 's' :
			if (N_TIME != preTime_sub) preTime_sub = N_TIME;
			if (N_ITER != preIt_sub) preIt_sub = N_ITER;
			break;
			
		case 'i' :
			if (N_TIME != preTime_intfc) preTime_intfc = N_TIME;
			if (N_ITER != preIt_intfc) preIt_intfc = N_ITER;
			break;
			
		default:
			Message0("(timing_and_counters_par.c::UpdateCounters) Invalid meshzone specified \n");
	}
}



	