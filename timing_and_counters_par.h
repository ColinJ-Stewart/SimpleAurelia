#ifndef TIMING_COUNTERS_H
#define TIMING_COUNTERS_H

#include <udf.h>
/* #include <unsteady.h> */

#include "global_var_par.h"

/* global */
extern int isPaused;

/* prototypes */
void Print_Timestep_Header(void);
void Print_Zone_Header(char meshZone);
double Get_EffectivePeriod(void);
int Get_N_RAMPSTEPS(double ramp_frac);
float Smootherstep(float edge0, float edge1, float x);
float Clamp(float x, float lowerlimit, float upperlimit);
double Get_SimTimeInCycle(char meshZone);
int NewTimeStepForThisZone(char zone);
int NewIterationForThisZone(char zone);
int NewIteration(void);
int NewTimeStep(void);
int Get_CycleNumber(void);
void Store_OldCycleN(void);
int NewCycle(void);
void UpdateCounters(char zone);


/* void endUpdateCounters(void);*/

#endif
