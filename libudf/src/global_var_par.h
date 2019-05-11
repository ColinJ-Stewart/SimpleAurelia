#ifndef GLOBAL_VAR_H
#define GLOBAL_VAR_H

#include <udf.h>
/* #include <unsteady.h> */

/* constants */
#define MAX_NODES 400			/* maximum number of nodes on the exum or subum surface*/
#define NUM_ZONES 2				/* number of zones that are being moved (i.e. number of zones that comprise the jellyfish)*/
#define EX_ZONE 26				/* mesh zone ID in fluent (found in BC menu)*/
#define SUB_ZONE 27				/* mesh zone ID in fluent (found in BC menu)*/


/* + global variables */
#if !RP_HOST
extern Domain *domain;
#endif

/* - simulation settings */
extern const int LET_IT_SWIM;
extern const double PAUSE_FRAC;

/* - morphometrics */
extern const double DIAMETER;
extern const double thickness;

/* - timing */
/* extern const double t_contract; */
/* extern const double t_relax; */
extern const double PERIOD;


/* - fluid constants */
extern const double NU;
extern const double RHO;


/* prototypes */


#endif

