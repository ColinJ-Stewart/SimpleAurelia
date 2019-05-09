#ifndef LOAD_KIN_H
#define LOAD_KIN_H

#include <udf.h>
/* #include <unsteady.h> */

#include "global_var_par.h"
/* #include "mesh_motion_par.h" */		/* needed for NodeIsTip() call */

/* constants */



/* prototypes */
void Load_sa(double t_norm, double *sa_ptr);
void Load_sb(double t_norm, double *sb_ptr);
void Load_sd(double t_norm, double *sd_ptr);
void Load_sg(double t_norm, double *sg_ptr);
void Load_a_sub(double *a_sub_ptr, double a, Node *holdNodes[], int nNodes, double SmaxUnflexed, int idArray[][2], 
				double coordArray_y_flex[], double arclengthArray_unflex[]);
void Load_b_sub(double *b_sub_ptr, double b, Node *holdNodes[], int nNodes, double SmaxUnflexed, int idArray[][2], 
				double coordArray_x_flex[], double arclengthArray_unflex[]);
#endif
