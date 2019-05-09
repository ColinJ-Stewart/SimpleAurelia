#include "load_kinematics_par.h"

/** global variables * */
static const int numSplines = 100;

static int readInCoeff_sa = 0;
static int readInCoeff_sb = 0;
static int readInCoeff_sd = 0;
static int readInCoeff_sg = 0;
static int readInCoeff_a_sub = 0;
static int readInCoeff_b_sub = 0;

static double ccoeff_sa[100][4];
static double ccoeff_sb[100][4];
static double ccoeff_sd[100][4];
static double ccoeff_sg[100][4];
static double ccoeff_a_sub[100][4];
static double ccoeff_b_sub[100][4];



/** Function definitions * */
/*---------- Load_sa ---------------------------------------
Purpose: 

Input:	
	tau = current time on [0,1]

Output: sa
---------------------------------------------------------------- */
void Load_sa(double t_norm, double *sa_ptr)
{
	int splineN, i, j, c;
	double saT = 0.0;
	double sold, coeffDelta, tDelta;
	char buffer1[1024],temp[400], line[70];
	FILE *finp, *fdebug;
	
	/* Message("\ngot to sa function\n"); */
	
	/*  Read in cubic spline coefficients from text files  */
	if (!readInCoeff_sa)
	{
		/* Open the input file stream  */
		#if !PARALLEL
			sprintf(buffer1, "SimpleAureliaSplineCoeff_s%c.txt", 'a');		
		#endif
		
		#if PARALLEL
			sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_s%c.txt", 'a');		
		#endif
		finp = fopen(buffer1,"r");	
		if(finp == NULL) 
			Error("Error opening files (sa).\n");
		
		/* read until EOF reached in input  */ 
		/* Message("\nSplines: \n");   */
		i = 1;
		while( i < numSplines + 5) 
		{
			if (i < 5) 
				fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
			else 
			{
				fgets(line, sizeof(line), finp);
				c = i-5;
				sscanf(line, "%lf %lf %lf %lf", &ccoeff_sa[c][3], &ccoeff_sa[c][2],
					&ccoeff_sa[c][1], &ccoeff_sa[c][0]);
				/* Message("%lf  %lf  %lf  %lf\n",ccoeff_sa[c][3],ccoeff_sa[c][2],ccoeff_sa[c][1],ccoeff_sa[c][0]);  */
			}
			i++;
		}
		fclose(finp);
		readInCoeff_sa = 1;
	}
	
	/* Message("\ngot to coefficient calc\n");  */
	/* compute s using the spline coefficients  */
	coeffDelta = 1.0/numSplines;
	/* coeffDelta = 0.01; */
	splineN =  floor(t_norm/coeffDelta);
	tDelta = fmod(t_norm, coeffDelta);
	for( j = 0; j < 4; j++) 
	{
		*sa_ptr += ccoeff_sa[splineN][j] * pow(tDelta,j);
	}
		
}



/*---------- Load_sb ---------------------------------------
Purpose: 

Input:	
	t_norm = current time on [0,1]

Output: sb
---------------------------------------------------------------- */
void Load_sb(double t_norm, double *sb_ptr)
{
	int splineN, i, j, c;
	/* double sb = 0; */
	double sold, coeffDelta, tDelta;
	char buffer1[1024],temp[400], line[70];
	FILE *finp;
	
	/* Message("\ngot to sb function\n"); */
	
	/*  Read in cubic spline coefficients from text files  */
	/*! Probably should move this function to be read in on loading (EXECUTE_ON_LOADING)  */
	if (!readInCoeff_sb)
	{
		/* Open the input file stream  */
		#if !PARALLEL
			sprintf(buffer1, "SimpleAureliaSplineCoeff_s%c.txt", 'b');		
		#endif
		
		#if PARALLEL
			sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_s%c.txt", 'b');
		#endif		
		finp = fopen(buffer1,"r");	
		if(finp == NULL) 
			Error("Error opening files. (sb)\n");
		
		/*read until EOF reached in input  */ 
		/* Message("\nSplines: \n"); */
		i = 1;
		while( i < numSplines + 5) 
		{
			if (i < 5) 
				fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
			else 
			{
				fgets(line, sizeof(line), finp);
				c = i-5;
				sscanf(line, "%lf %lf %lf %lf", &ccoeff_sb[c][3], &ccoeff_sb[c][2],
					&ccoeff_sb[c][1], &ccoeff_sb[c][0]);
				/* Message("%lf  %lf  %lf  %lf\n",ccoeff_sb[c][3],ccoeff_sb[c][2],ccoeff_sb[c][1],ccoeff_sb[c][0]); */
			}
			i++;
		}
		fclose(finp);
		readInCoeff_sb = 1;
	}
	
	/* compute s using the spline coefficients  */
	coeffDelta = 1.0 / numSplines;
	splineN =  floor(t_norm/coeffDelta);
	tDelta = fmod(t_norm, coeffDelta);
	for( j = 0; j < 4; j++) 
	{
		*sb_ptr += ccoeff_sb[splineN][j] * pow(tDelta,j);
	}
	/* Message("\ncoeffDelta = %lf, splineN = %i, tDelta = %lf, sb = %lf\n",coeffDelta,splineN,tDelta,*sb_ptr); */

}



/*---------- Load_sd ---------------------------------------
Purpose: 

Input:	
	t_norm = current time on [0,1]

Output: sd
---------------------------------------------------------------- */
void Load_sd(double t_norm, double *sd_ptr)
{
	int splineN, i, j, c;
	/* double sd = 0; */
	double sold, coeffDelta, tDelta;
	char buffer1[1024],temp[400], line[70];
	FILE *finp;
	
	/* Message("\ngot to sd function\n"); */
	
	/*  Read in cubic spline coefficients from text files  */
	/*! Probably should move this function to be read in on loading (EXECUTE_ON_LOADING)  */
	if (!readInCoeff_sd)
	{
		/* Open the input file stream  */
		#if !PARALLEL
			sprintf(buffer1, "SimpleAureliaSplineCoeff_s%c.txt", 'd');		
		#endif
		
		#if PARALLEL
			sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_s%c.txt", 'd');
		#endif	
		finp = fopen(buffer1,"r");	
		if(finp == NULL) 
			Error("Error opening files. (sd)\n");
		
		/*read until EOF reached in input  */ 
		/* Message("\nSplines: \n"); */
		i = 1;
		while( i < numSplines + 5) 
		{
			if (i < 5) 
				fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
			else 
			{
				fgets(line, sizeof(line), finp);
				c = i-5;
				sscanf(line, "%lf %lf %lf %lf", &ccoeff_sd[c][3], &ccoeff_sd[c][2],
					&ccoeff_sd[c][1], &ccoeff_sd[c][0]);
				/* Message("%lf  %lf  %lf  %lf\n",ccoeff_sd[c][3],ccoeff_sd[c][2],ccoeff_sd[c][1],ccoeff_sd[c][0]); */
			}
			
			i++;
		}
		
		fclose(finp);
		readInCoeff_sd = 1;
	}
	
	/* compute s using the spline coefficients  */
	coeffDelta = 1.0 / numSplines;
	splineN =  floor(t_norm/coeffDelta);
	tDelta = fmod(t_norm, coeffDelta);
	for( j = 0; j < 4; j++) 
	{
		*sd_ptr += ccoeff_sd[splineN][j] * pow(tDelta,j);
	}
	/* Message("\ncoeffDelta = %lf, splineN = %i, tDelta = %lf, sd = %lf\n",coeffDelta,splineN,tDelta,*sd_ptr); */

}



/*---------- Load_sg ---------------------------------------
Purpose: 

Input:	
	t_norm = current time on [0,1]

Output: sg
---------------------------------------------------------------- */
void Load_sg(double t_norm, double *sg_ptr)
{
	int splineN, i, j, c;
	/* double sg = 0; */
	double sold, coeffDelta, tDelta;
	char buffer1[1024], temp[400], line[70];
	FILE *finp;
	
	/* Message("\ngot to sg function\n"); */
	
	/*  Read in cubic spline coefficients from text files  */
	/*! Probably should move this function to be read in on loading (EXECUTE_ON_LOADING)  */
	if (!readInCoeff_sg)
	{
		/* Open the input file stream  */
		#if !PARALLEL
			sprintf(buffer1, "SimpleAureliaSplineCoeff_s%c.txt", 'g');		
		#endif
		
		#if PARALLEL
			sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_s%c.txt", 'g');
		#endif	
		finp = fopen(buffer1,"r");	
		if(finp == NULL) 
			Error("Error opening files. (sg)\n");
		
		/*read until EOF reached in input  */ 
		/* Message("\nSplines: \n"); */
		i = 1;
		while( i < numSplines + 5) 
		{
			if (i < 5) 
				fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
			else 
			{
				fgets(line, sizeof(line), finp);
				c = i-5;
				sscanf(line, "%lf %lf %lf %lf", &ccoeff_sg[c][3], &ccoeff_sg[c][2],
					&ccoeff_sg[c][1], &ccoeff_sg[c][0]);
				/* Message("%lf  %lf  %lf  %lf\n",ccoeff_sg[c][3],ccoeff_sg[c][2],ccoeff_sg[c][1],ccoeff_sg[c][0]); */
			}
			
			i++;
		}
		
		fclose(finp);
		readInCoeff_sg = 1;
	}
	
	/* compute s using the spline coefficients  */
	coeffDelta = 1.0 / numSplines;
	splineN =  floor(t_norm/coeffDelta);
	tDelta = fmod(t_norm, coeffDelta);
	for( j = 0; j < 4; j++) 
	{
		*sg_ptr += ccoeff_sg[splineN][j] * pow(tDelta,j);
	}
	/* Message("\ncoeffDelta = %lf, splineN = %i, tDelta = %lf, sg = %lf\n",coeffDelta,splineN,tDelta,*sg_ptr); */

}



/*---------- Load_a_sub ------------------------------------
Purpose: Returns the tparm value for both the ex and sub surfaces

Input:	tparm_ptr	-	pointer to tparm so it can be modified by
						the function.
		meshZone	-	'e' for exumbrella zone, 's' for subumbrella
		holdNodes	-	array of sorted nodes, from apex to margin
		a			-	kinematic variable a for this time step
		a_sub		-	kinematic variable a_sub (only for sub-
						umbrella. Exumbrella will pass zero.
		
Output:		none. All changes are made to tparm via tparm_ptr
---------------------------------------------------------------- */
void Load_a_sub(double *a_sub_ptr, double a, Node *holdNodes[], int nNodes, double SmaxUnflexed, 
				int idArray[][2], double coordArray_y_flex[], double arclengthArray_unflex[]) 
{
    int splineN, i, j, c;
	double coeffDelta, Si, Snorm, del, dif, S_End, S_EndNorm;
	char buffer1[1024], temp[400], line[200];
	Node *v;
	FILE *finp;
	int thisID, localj;
	
	
	if (N_TIME == 0) 
	{
		/*  Read in cubic spline coefficients from text files  */
		/*! Probably should move this function to be read in on loading (EXECUTE_ON_LOADING)  */
		if (!readInCoeff_a_sub)
		{
			/* Open the input file stream  */
			#if !PARALLEL
				sprintf(buffer1, "SimpleAureliaSplineCoeff_%csub.txt",'a');		
			#endif
			
			#if PARALLEL
				sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_%csub.txt",'a');
			#endif	
			finp = fopen(buffer1,"r");	
			if(finp == NULL) 
				Error("Error opening files. (a_b)\n");
			
			/*read until EOF reached in input  */ 
			i = 1;
			while( i < numSplines + 5) 
			{
				if (i < 5)
					fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
				else 
				{
					fgets(line, sizeof(line), finp);
					c = i-5;
					sscanf(line, "%lf %lf %lf %lf", &ccoeff_a_sub[c][3], &ccoeff_a_sub[c][2],
						&ccoeff_a_sub[c][1], &ccoeff_a_sub[c][0]);
				}
				i++;
			}
			fclose(finp);
			readInCoeff_a_sub = 1;
		}
		
		
		/* Loop over nodes, computing a_sub using the spline coefficients  */
		coeffDelta = 1.0 / numSplines;
		for ( j = 0; j < nNodes; j++) 
		{
			/* New code  */
			if (j == nNodes - 1) /* is tip */
			{
				Si = SmaxUnflexed;
				Snorm = 1.0;
				splineN = numSplines-1;
				del = coeffDelta;
			}
			else
			{
				Si = arclengthArray_unflex[j];
				Snorm = Si/SmaxUnflexed;
				splineN =  floor(Snorm/coeffDelta);
				del = fmod(Snorm, coeffDelta);
			}
			
			a_sub_ptr[j] = 0.0;					/* initialize */
				
			for ( i = 0; i < 4; i++) 
				a_sub_ptr[j] += ccoeff_a_sub[splineN][i] * pow(del,i);
			
		}
	}

}


/*---------- Load_b_sub ------------------------------------
Purpose: Returns the tparm value for both the ex and sub surfaces

Input:	tparm_ptr	-	pointer to tparm so it can be modified by
						the function.
		meshZone	-	'e' for exumbrella zone, 's' for subumbrella
		holdNodes	-	array of sorted nodes, from apex to margin
		b			-	kinematic variable b for this time step
		b_sub		-	kinematic variable b_sub (only for sub-
						umbrella. Exumbrella will pass zero.
		
Output:		none. All changes are made to tparm via tparm_ptr
---------------------------------------------------------------- */
void Load_b_sub(double *b_sub_ptr, double b, Node *holdNodes[], int nNodes, double SmaxUnflexed, 
				int idArray[][2], double coordArray_x_flex[], double arclengthArray_unflex[])
{
    int splineN, i, j, c;
	double coeffDelta, Si, Snorm, del, dif, this_input, S_End, S_EndNorm;
	char buffer1[1024], temp[400], line[200];
	Node *v;
	FILE *finp;
	int thisID, localj, j_err, jInterpEnd;
	
	
	/* Message("\ngot to a_sub function\n");  */
	if (N_TIME == 0) 
	{
		/*  Read in cubic spline coefficients from text files  */
		/*! Probably should move this function to be read in on loading (EXECUTE_ON_LOADING)  */
		if (!readInCoeff_b_sub)
		{
			/* Open the input file stream  */
			#if !PARALLEL
				sprintf(buffer1, "SimpleAureliaSplineCoeff_%csub.txt",'b');		
			#endif
			
			#if PARALLEL
				sprintf(buffer1, "/tmp/SimpleAureliaSplineCoeff_%csub.txt",'b');
			#endif
			finp = fopen(buffer1,"r");	
			if(finp == NULL) 
				Error("Error opening files. (b_sub)\n");
			
			/* read until EOF reached in input  */ 
			i = 1;
			while( i < numSplines + 5) 
			{
				if (i < 5) 
					fgets(temp, sizeof(temp), finp); /* read the 4 header lines; do nothing  */
				else 
				{
					fgets(line, sizeof(line), finp);
					c = i-5;
					sscanf(line, "%lf %lf %lf %lf", &ccoeff_b_sub[c][3], &ccoeff_b_sub[c][2],
						&ccoeff_b_sub[c][1], &ccoeff_b_sub[c][0]);

				}
				i++;
			}
			fclose(finp);
			readInCoeff_b_sub = 1;
		}
		
		/** (Loop1) * */
		Message0("\nLoop1\n");
		/* Loop over nodes, computing b_sub using the spline coefficients  */
		coeffDelta = 1.0 / numSplines;
		for (j = 0; j < nNodes; j++) 
		{
			/* New code  */
			if (j == nNodes - 1) /* is tip */
			{
				Si = SmaxUnflexed;
				Snorm = 1.0;
				splineN = numSplines-1;
				del = coeffDelta;
			}
			else
			{
				Si = arclengthArray_unflex[j];
				Snorm = Si/SmaxUnflexed;
				splineN =  floor(Snorm/coeffDelta);
				del = fmod(Snorm, coeffDelta);
			}
			
			b_sub_ptr[j] = 0.0;					/* initialize */
				
			for ( i = 0; i < 4; i++) 
				b_sub_ptr[j] += ccoeff_b_sub[splineN][i] * pow(del,i);
		
		}
	}
	
}





