#ifndef MESH_MOTION_H
#define MESH_MOTION_H

#include <udf.h>
/* #include <unsteady.h> */
#include <math.h>
#include <float.h>

#include "global_var_par.h"
#include "load_kinematics_par.h"				/* needed for all I/O functions */
#include "timing_and_counters_par.h"			/* needed for NewTimeStepForThisZone() */
#include "user_mem_mgmt_par.h"					/* needed for Reinit_node_mem_int() */
#include "quicksort_par.h"

/* constants */

/* global variables */




/* DEBUGGING prototypes */
void Calc_Mesh_Movement_debug(Thread *tf, char meshZone, double tau, double del_x, double vel);
void Calc_Kinematics_and_Move_debug(Thread *tf, char meshZone, double tau, double del_x);
void Move_Nodes_to_Stored_Positions_debug(Thread *tf, int x_memslot, int y_memslot);

/* prototypes */
void Calc_Mesh_Movement(Thread *tf, char meshZone, double tau, double del_x, double vel);

void Calc_Kinematics_and_Move(Thread *tf, char meshZone, double tau, double del_x);

int Get_ArcLengths(Thread *tf, char meshZone, Node *holdNodes[], 
		double *SmaxUnflexed_ptr, int idArray[][2], double arclengthArray_unflex[]);
					
double Get_Unflexed_ArcLengths(char meshZone, Node *holdNodes[], int nNodes);

double Rezero(Thread *tf);

int Get_NodeDistances(Thread *tf, char meshZone, Node *holdNodes[], double xApex);

void Get_ParArcLengths(char meshZone, Node *holdNodes[], int i, 
		double *SmaxUnflexed_ptr, int idArray[][2], double arclengthArray_unflex[] );
						
 /* void Get_ParArcLengths(char meshZone, Node *holdNodes[], int i,  */
					 /* double distArray[],int idArray[][2],  */
					 /* double coordArray_x_flex[], double coordArray_y_flex[],  */
					 /* double coordArray_x_unflex[], double coordArray_y_unflex[],  */
					 /* double arclengthArray_flex[], double arclengthArray_unflex[],  */
					 /* double *Smax_ptr, double *SmaxUnflexed_ptr) */
					
int NodeIsTip (Node *v);

void Get_tparm(char meshZone, Node *holdNodes[], double a, double b, double b_sub[], int nNodes, 
				double SmaxUnflexed, int idArray[][2]);
				
void Get_a_sub(double *a_sub_ptr, double a, Node *holdNodes[], int nNodes, 
				double SmaxUnflexed, int idArray[][2]);
				
void Get_b_sub(double *b_sub_ptr, double b, Node *holdNodes[], int nNodes, 
				double SmaxUnflexed, int idArray[][2]);

void Get_FlexPoint(char meshZone, Node *holdNodes[], int nNodes, double SmaxUnflexed, 
					double sg, double g[], int idArray[][2]);

void Move_Nodes_to_Stored_Positions(Thread *tf, int x_memslot, int y_memslot);

void Store_OldKinematicVars(void);

#endif

