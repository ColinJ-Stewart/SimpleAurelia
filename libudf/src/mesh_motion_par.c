#include "mesh_motion_par.h"

/** global variables **/
static const double ai = 0.0775;							/* initial/max radius [m] */
static const double bi = 0.039336;							/* initial bell apex-centroid height [m] */
static const double di = 0.00447;							/* initial centroid-bell margin height [m] */

static const double pa = 0.258064516129032;					/* y-contraction % = 1-(D_contract/D_relax) */
static const double pb = 0.198;								/* bell apex-centroid contraciton % = 1.1*0.18 */
static const double pdd = -4.275382550335571;				/* centroid-bell margin contraction %  */

static const double g1 = 1.0;								/* parameters for flex */
static const double g2 = 1.15;
static const double pg = 0.5;
static const double k = 8.0;

static double xTip = 0.0;
static double yTip = 0.0;
static double xTipOld = 0.0;
static double yTipOld = -1.0;
static double xTipUnflexed = 0.0;
static double yTipUnflexed = 0.0;
static double xTipUnflexedOld = 0.0;
static double yTipUnflexedOld = 0.0;

static double xTip_temp = 0.0;
static double yTip_temp = 0.0;
static double xTipUnflexed_temp = 0.0;
static double yTipUnflexed_temp = 0.0;

static double tparmTip = 0.0;
static double tparmTipOld = 0.0;
static double tparm_ex[MAX_NODES];
static double tparm_ex_old[MAX_NODES];
static double tparm_sub[MAX_NODES];
static double tparm_sub_old[MAX_NODES];

static int i1_ex = -1;
static int i1_sub = -1;

/* flags  */
static int newNodeFlag = 0;
static int iHaveTip = 0;


/** Function definitions **/

/*---------- Calc_Mesh_Movement() --------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
void Calc_Mesh_Movement(Thread *tf, char meshZone, double tau, double del_x, double vel) 
{
	/* Message("(before) tip = (%f, %f) tipUnflexed = (%f, %f)\n", xTip, yTip, xTipUnflexed, yTipUnflexed); */
	/* Message("(before) tip_temp = (%f, %f) tipUnflexed_temp = (%f, %f)\n", xTip_temp, yTip_temp, xTipUnflexed_temp, yTipUnflexed_temp); */
	#if !RP_HOST
	
	if (NewTimeStepForThisZone(meshZone))
	{
		Calc_Kinematics_and_Move(tf, meshZone, tau, del_x);
	}
	else
		Move_Nodes_to_Stored_Positions(tf, 4, 5);	/* Node positions (x,y) stored in node memory (4,5) */
	
	#endif 
}


/*---------- Calc_Kinematic_Motion() --------------------------------------
Purpose: 

Input:	

Output:
---------------------------------------------------------------- */
void Calc_Kinematics_and_Move(Thread *tf, char meshZone, double tau, double del_x) 
{
	#if !RP_HOST
	face_t f;
	Node *v;
	Node *holdNodes[MAX_NODES];
	double a, b, d;
	double a_sub[MAX_NODES], b_sub[MAX_NODES], g[MAX_NODES];
	double *a_sub_ptr = a_sub;
	double *b_sub_ptr = b_sub;
	double sa = 0.0;
	double sb = 0.0;
	double sd = 0.0;
	double sg = 0.0;
	double *sa_ptr = &sa;
	double *sb_ptr = &sb;
	double *sd_ptr = &sd;
	double *sg_ptr = &sg;
	
	double xold, yold;
	int i, j, n;
	int i1 = 0;
	FILE *fdebug;
	double SmaxUnflexed;
	double *SmaxUnflexed_ptr = &SmaxUnflexed;
	int total_i = 0;
	double x = Rezero(tf);	
	int idArray[MAX_NODES][2];
	int thisID, localj;
	double arclengthArray_unflex[MAX_NODES];
	
	/* initialize kinematics arrays  */
	for (j=0; j < MAX_NODES; j++) 
	{
		a_sub[j] = 0;
		b_sub[j] = 0;
		g[j] = 0.0;
	}

	/* calculate the current values of kinematic functions sa, sb, sd  */
	/*Message0("\tCalculating kinematics (s-functions):  ");*/
	Load_sa(tau/PERIOD, sa_ptr);
	Load_sb(tau/PERIOD, sb_ptr);
	Load_sd(tau/PERIOD, sd_ptr);
	Load_sg(tau/PERIOD, sg_ptr);
	
	/* calculate the function values for a, b, d  */
	/*Message0("\tCalculating kinematics (a, b, d):  ");*/
	a = ai*(1-sa*pa);
	b = bi*(1-sb*pb);
	d = di*(1-sd*pdd);
	
	/* calculate the new position of the margin tip. This is used to find arc length as well as tparm  */
	tparmTip = acos(d/b);
	
	if (N_TIME == 0) 
	{
		xTipUnflexed = b*cos(tparmTip) + b;
		yTipUnflexed = a*sin(tparmTip);
		xTipUnflexedOld = xTipUnflexed;
		yTipUnflexedOld = yTipUnflexed;
		
		xTip =  xTipUnflexed;
		yTip = yTipUnflexed;
		xTipOld = xTipUnflexed;
		yTipOld = yTipUnflexed;
	}
	Message0("\tReporting tip values: xTipUnflexed = %lf, yTipUnflexed = %lf, tparmTip = %lf \n", xTipUnflexed, yTipUnflexed, tparmTip);
	Message0("\tReporting tip values: xTipUnflexedOld = %lf, yTipUnflexedOld = %lf \n", xTipUnflexedOld, yTipUnflexedOld);
	
	
	/* ----- calculate ARC LENGTHS and SORT NODES from apex to margin ----------  */		
	Message0("\tCalculating surface arc length... \n");
	i = Get_ArcLengths(tf, meshZone, holdNodes, SmaxUnflexed_ptr, idArray, arclengthArray_unflex);	/* i = number of nodes on mesh zone */
	#if PARALLEL
		total_i = PRF_GISUM1(i);
	#endif
	
	/* Message("(i = %i nodes)\n",i); */
	/* Message("Done.\n"); */
	
	/*------------ calculate tparm, asub, and bsub ------------ */
	if (meshZone == 'e') { 
		Message0("\tCalculating initial tparm values... \n");
		Get_tparm(meshZone, holdNodes, b, b_sub, total_i, SmaxUnflexed, idArray);
		Message0("\tDone.\n");
	}
	else if (meshZone == 's') {
		if  (N_TIME == 0) 
		{
			Message0("\tLoading initial a_sub and b_sub values... \n");
			Load_b_sub(b_sub_ptr, b, holdNodes, total_i, SmaxUnflexed, idArray, arclengthArray_unflex);
			Load_a_sub(a_sub_ptr, a, holdNodes, total_i, SmaxUnflexed, idArray, arclengthArray_unflex);
			
			Message0("\tCalculating initial tparm values... \n");
			Get_tparm(meshZone, holdNodes, b, b_sub, total_i, SmaxUnflexed, idArray);
			/* Message("Done.\n"); */
		}
		else 
		{
			/* Message("!! tparmTip = %lf, tparmTipOld = %lf!!\n", tparmTip, tparmTipOld);  */
			/*Message0("\tCalculating tparm... \n");*/
			Get_tparm(meshZone, holdNodes, b, b_sub, total_i, SmaxUnflexed, idArray);
			
			/*Message0("\tCalculating a_sub... \n");*/
			Get_a_sub(a_sub_ptr, a, holdNodes, total_i, SmaxUnflexed, idArray);
			
			/*Message0("\tCalculating b_sub... \n");*/
			Get_b_sub(b_sub_ptr, b, holdNodes, total_i, SmaxUnflexed, idArray);
			/* Message("Done.\n"); */
		}
	}
	else Error("Incorrect meshZone char");
	
	/*------------ find flex point (node) ------------ */	 
	Message0("\tFinding flex points... \n");
	Get_FlexPoint(meshZone, holdNodes, total_i, SmaxUnflexed, sg, g, idArray);
	
	/*--- Now that arclength, tparm, flex, etc. have been assigned to --- */
	/*--- any remeshed nodes, reset the new node flag.				  --- */
	Reinit_node_mem_int(tf, 9, 1);		/* Setting N_UDMI slot 9 to 1 (i.e. N_UDMI(v,9) = 1) */
	
	Message0("\tMove Nodes Loop. ");
	Message0("(x = %lf, del_x = %lf)\n", x, del_x);
	/* Message("xTip = %f, yTip = %f, xTipUnflexed = %f, yTipUnflexed = %f\n", xTip, yTip, xTipUnflexed, yTipUnflexed); */
	for (j = 0; j < total_i; j++) 
	{
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */

		if (myid == thisID)
		{
			v = holdNodes[localj];
			if (NODE_POS_NEED_UPDATE(v) || (NODE_Y(v) <= 1e-7))
			{	
				/* store old position  */
				N_UDMI(v,2) = NODE_X(v); 	/* xold */
				N_UDMI(v,3) = NODE_Y(v);	/* yold */
				
				/* kinematic motion  */
				/* Message("j:%i  v:%i   g = %lf\n", j, v, g[j]); */
				/* g[j] = 1.0; */
				/* calc g for this node  */				
				if (meshZone == 'e') {
					/* flexed position  */
					N_UDMI(v,4) = x + del_x + (g[j]*b*cos(tparm_ex[j]) + b);	/* new x position = current position + forward swimming displacement + kinematic motion */
					N_UDMI(v,5) = g[j]*a*sin(tparm_ex[j]);						/* new y position = current position + kinematic motion	 */
					/* unflexed position  */
					N_UDMI(v,12) = x + del_x + (b*cos(tparm_ex[j]) + b);	/* g[j] term omitted  */
					N_UDMI(v,13) = a*sin(tparm_ex[j]);						/* g[j] term omitted  */
				}
				else if (meshZone == 's') {
					/* flexed position  */
					N_UDMI(v,4) = (x-thickness) + del_x + (g[j]*b_sub[j]*cos(tparm_sub[j]) + b);
					N_UDMI(v,5) = g[j]*a_sub[j]*sin(tparm_sub[j]);
					/* unflexed position  */
					N_UDMI(v,12) = (x-thickness) + del_x + (b_sub[j]*cos(tparm_sub[j]) + b);
					N_UDMI(v,13) = a_sub[j]*sin(tparm_sub[j]);
					
					/* if ((myid == 1) || (myid == 0)) */
					/* { */
						/* Message("j:%i  (x-thickness)=%f, as=%f, bs=%f, tparm=%f\n", */
							/* j, x-thickness, a_sub[j], b_sub[j], tparm_sub[j]); */
						
					/* } */
				}
				
				/* Store new location of flexed tip  */
				/* if (NodeIsTip(v)) */
				if ( NodeIsTip(v) && (meshZone == 'e') )
				{
					xTip_temp = N_UDMI(v,4);
					yTip_temp = N_UDMI(v,5);
					xTipUnflexed_temp = N_UDMI(v,12);
					yTipUnflexed_temp = N_UDMI(v,13);
					iHaveTip = 1;
				}
				
				/* -- Move nodes except apex --  */
				if (NODE_Y(v) > 1e-7)
				{
					NODE_X(v) = N_UDMI(v,4);
					NODE_Y(v) = N_UDMI(v,5);
				}
				NODE_POS_UPDATED(v);
				/** DEBUG * */
				/* if ( (N_TIME == 1)  || (N_TIME >= 63) ) */
				/* { */
				/* Message("\tj:%i  xold = %lf\t\t yold = %lf (node%i)\n",j, N_UDMI(v,2), N_UDMI(v,3), myid); */
				/* Message("\tj:%i  xnew = %lf\t\t ynew = %lf (node%i)\n",j, N_UDMI(v,4), N_UDMI(v,5), myid); */
				/* Message("\tj:%i  delx = %lf\t\t dely = %lf (node%i)\n\n",j, N_UDMI(v,4)-N_UDMI(v,2), N_UDMI(v,5)-N_UDMI(v,3), myid); */ 
				/* }
				
			}
			else
			{
				/* Message("\tAlready updated: j:%i  xold = %lf\t\t yold = %lf\n",j, N_UDMI(v,2), N_UDMI(v,3)); */
				/* Message("\tAlready updated: j:%i  xnew = %lf\t\t ynew = %lf\n",j, N_UDMI(v,4), N_UDMI(v,5)); */
				/* debug */
				/* Message("Already updated: j:%i  delx = %lf\t\t dely = %lf\n\n",j, N_UDMI(v,4)-N_UDMI(v,2), N_UDMI(v,5)-N_UDMI(v,3)); */
			}
		}
	}

	#endif
}

/*---------- Move_Nodes_to_Stored_Positions ----------------------------
Purpose: Moves the nodes according to kinematic and dynamic motions
calculated in calcMeshMovement()

Input:	holdNodes	-	array of sorted nodes, from apex to margin
		i 			-	total number of nodes in this zone
		i1			-	index of flex point node in holdNodes
		a_sub		-	kinematic variable
		b_sub		-	kinematic variable
		
Output:	
---------------------------------------------------------------- */
void Move_Nodes_to_Stored_Positions(Thread *tf, int x_memslot, int y_memslot)
{
	#if !RP_HOST
	face_t f;
	Node *v;
	int n;
	int j = 0;
	
	
	/* -- Move nodes --  */
	Message0("\tMoving nodes to stored positions...\n");
	j = 0;
	begin_f_loop(f,tf) {
		if PRINCIPAL_FACE_P(f,tf) {
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);
				
				if (NODE_POS_NEED_UPDATE(v) && (NODE_Y(v) > 1e-7))
				{	/* this was commented out! may break the simulation, keep an eye out */
					NODE_X(v) = N_UDMI(v, x_memslot);
					NODE_Y(v) = N_UDMI(v, y_memslot);
					NODE_POS_UPDATED(v);
					
					
					/*Message("\nnode %i:  dist = %lf\n",j,N_UDMI(v,0)); */
					/* Message("node %i:  xold = %lf\t\t yold = %lf\n",v, N_UDMI(v,2), N_UDMI(v,3)); */
					/* Message("node %i:  xnew = %lf\t\t ynew = %lf\n",v, N_UDMI(v,4), N_UDMI(v,5)); */
					/* Message("node %i:  delx = %lf\t\t dely = %lf\n\n",j, N_UDMI(v,4)-N_UDMI(v,2), N_UDMI(v,5)-N_UDMI(v,3)); */
					j++;
				}
				/* else */
				/* { */
					/* Message("Already updated node %i:  xold = %lf\t\t yold = %lf\n",v, N_UDMI(v,2), N_UDMI(v,3)); */
					/* Message("Already updated node %i:  xnew = %lf\t\t ynew = %lf\n",v, N_UDMI(v,4), N_UDMI(v,5)); */
					
				/* } */
			}
		}
	}
	end_f_loop(f,tf);
	
	#endif
}


/*---------- Get_ArcLengths---------------------------------------
Purpose: Use to sorting approach to order the nodes after 	
finding the lengths to the tip. Array is sorted from margin to apex,
holdNodes[0] is the margin, and holdNodes[i-1] is the apex. Arc length
is stored from Smax at [0] and  Smin=0 at [i-1].

Input:	Thread *tf			-			
		Node *holdNodes		-	
		double xTip			-	
		double yTip			-
---------------------------------------------------------------- */
int Get_ArcLengths(Thread *tf, char meshZone, Node *holdNodes[], 
		double *SmaxUnflexed_ptr, int idArray[][2], double arclengthArray_unflex[])		
{
	Message0("\tGet_ArcLengths::first line\n");
	int i = 0;

	#if !RP_HOST
		double xApex;
		int j, n;
		double dist;
		Node *v, *w;
		
		/* Message("\tGet_ArcLengths::Rezero\n"); */
		Message0("\tGet_ArcLengths::Rezero\n");
		xApex = Rezero(tf);
		
		/* Loop over all mesh nodes on this computational node, put them in holdNodes array, 
			and store distance to tip if first iteration of sim  */
		Message0("\tGet_ArcLengths::Get_NodeDistance\n");
		i = Get_NodeDistances(tf, meshZone, holdNodes, xApex);
		
		/* sort mesh nodes and store in holdNodes[0:i] where i is local to this computational node  */
		Message0("\tGet_ArcLengths::quickSort\n");
		quickSort(holdNodes, 0, i-1, meshZone); 
		
		
		#if PARALLEL
			/* combines local node lists into ordered, global array from which the flexed and unflexed
				arc lengths are calculated. Other ordered arrays of node info are also avaible as results 
				(coordinates, radial distance, local arc lengths) but most of this info is stored at each node
				via N_UDMI and can be accessed using the idArray as a directory/lookup table  */
				/**
				 * ! This is one of the bugged parallel funcs!
				 * TODO: FIX THIS FUNCTION TO WORK WHEN A PROCESSOR HAS GROUPS OF NODES THAT AREN'T ALL NEXT TO EACH OTHER 
				 **/
			Message0("\tGet_ArcLengths::Get_ParArcLengths\n");
			Get_ParArcLengths(meshZone, holdNodes, i, SmaxUnflexed_ptr, idArray, arclengthArray_unflex);	
			Message0("\tGet_ArcLengths::Get_ParArcLengths -- complete\n");
		#endif

		
		#if !PARALLEL
			/* do it the old fashioned way...  */
		
			/* First iteration of arc length calc  */
			N_UDMI(holdNodes[0],1) = 0;
			idArray[0][0] = myid;
			idArray[0][1] = myid;
			/* Message("arclength = %lf, dist = %lf \n",N_UDMI(holdNodes[0],1), 0); */
			for (j = 1; j < i; j++) { 
				v = holdNodes[j-1];
				w = holdNodes[j]; 
				
				/* calculate distance from previous node  */
				if (NodeIsTip(w)) 
				{
					dist = sqrt(SQR(xTipOld - NODE_X(v)) + SQR(yTipOld - NODE_Y(v)));
					/* Message("\n**NODE IS TIP -- xTipOld = %lf, yTipOld = %lf***\n", xTipOld, yTipOld); */
					/* Message("arclength = %lf, dist = %lf \n",N_UDMI(w,1), dist); */
				}
				else 
				{ 
					dist = sqrt(SQR(NODE_X(w) - NODE_X(v)) + SQR(NODE_Y(w) - NODE_Y(v)));
					N_UDMI(w,1) = N_UDMI(v,1) + dist;  /* store arc length up to that node (cumulative)  */
					/* Message("arclength = %lf, dist = %lf \n",N_UDMI(w,1), dist); */
				}
				
				/* Message("j: %i  x = %lf, y = %lf, S = %lf, dist = %lf\n", j, NODE_X(w), NODE_Y(w), N_UDMI(w,1), dist); */
				
				idArray[j][0] = myid;
				idArray[j][1] = j;
			}
			
			/* Now calculate unflexed lengths  */
			*SmaxUnflexed_ptr = Get_Unflexed_ArcLengths(meshZone, holdNodes, i);
		#endif
	#endif

	return i;
}

/*---------- Get_Unflexed_ArcLengths---------------------------------------
Purpose: Use to sorting approach to order the nodes after 	
finding the lengths to the tip. Array is sorted from margin to apex,
holdNodes[0] is the margin, and holdNodes[i-1] is the apex. Arc length
is stored from Smax at [0] and  S=0 at [i-1].

Input:	Node *holdNodes		-	
		double xTip			-	
		double yTip			-
---------------------------------------------------------------- */
double Get_Unflexed_ArcLengths(char meshZone, Node *holdNodes[], int nNodes) 
{
	double SmaxUnflexed = 0.0;
	
	#if !RP_HOST
	int i, j, n;
	double dist, x;
	face_t f;
	Node *v, *w;
	
	
	/* First iteration of arc length calc  */
	N_UDMI(holdNodes[0],14) = 0;
	/* Message("arclength = %lf, dist = %lf \n",N_UDMI(holdNodes[0],1), 0); */
	for (j = 1; j < nNodes; j++) { 
		v = holdNodes[j-1];
		w = holdNodes[j]; 
		
		/* calculate distance from previous node  */
		if (NodeIsTip(w)) 
		{
			dist = sqrt(SQR(xTipUnflexedOld - N_UDMI(v,12)) + SQR(yTipUnflexedOld - N_UDMI(v,13)));
			N_UDMI(w,14) = N_UDMI(v,14) + dist;
			SmaxUnflexed = N_UDMI(w,14);
			Message0("\tNODE IS TIP -- S = SmaxUnflexed = %lf, xTipUnflexedOld = %lf, yTipUnflexedOld = %lf***\n",
					SmaxUnflexed, xTipUnflexedOld, yTipUnflexedOld);
			/* Message("arclength = %lf, dist = %lf \n",N_UDMI(w,1), dist); */
		}
		else 
		{ 
			dist = sqrt(SQR(N_UDMI(w,12) - N_UDMI(v,12)) + SQR(N_UDMI(w,13) - N_UDMI(v,13)));
			N_UDMI(w,14) = N_UDMI(v,14) + dist;  /* store arc length up to that node (cumulative)  */
			/* Message("arclength = %lf, dist = %lf \n",N_UDMI(w,1), dist); */
		}
		
		/* Message("unflex j: %i  x = %lf, y = %lf, S = %lf, dist = %lf\n", j, N_UDMI(w,12), N_UDMI(w,13), N_UDMI(w,14), dist); */
	}
	#endif
	
	return SmaxUnflexed;
}


/*---------- Rezero ---------------------------------------------
Purpose: 

Input:	Thread *tf - face thread for zone which touches axis of symmetry y=0
---------------------------------------------------------------- */
double Rezero (Thread *tf) {
	
	/* NEW  */
	double xapex = 0.0;
	
	#if !RP_HOST
	face_t f;
	int n, N, *k, *iworkN, i;
	Node *v;
	double *ylow, *xpos, *dworkN, ylowest;
	
	/* allocate memory for dynamic arrays  */
	ylow = (double *)malloc(sizeof(double)*compute_node_count);	
	if (ylow == NULL)
		Error("malloc failed!\n");
	
	xpos = (double *)malloc(sizeof(double)*compute_node_count);	
	if (xpos == NULL)
		Error("malloc failed!\n");
	
	dworkN = (double *)malloc(sizeof(double)*compute_node_count);	
	if (dworkN == NULL)
		Error("malloc failed!\n");
	
	k = (int *)malloc(sizeof(int)*compute_node_count);	
	if (k == NULL)
		Error("malloc failed!\n");
	
	iworkN = (int *)malloc(sizeof(int)*compute_node_count);	
	if (iworkN == NULL)
		Error("malloc failed!\n");
	
	
	/* initialize  */
	for (N = 0; N < compute_node_count; N++)	
	{
		ylow[N] = 0.0;
		xpos[N] = 0.0;
		dworkN[N] = 0.0;
		k[N] = 0;
		iworkN[N] = 0;
	}
	
	
	/* loop through mesh cells and compare node locations, storing the lowest y-val  */
	k[myid] = 0;
	begin_f_loop(f,tf) 
	{  
		if PRINCIPAL_FACE_P(f,tf)
		{
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);
				if (k[myid] == 0)
				{
					ylow[myid] = NODE_Y(v);
					xpos[myid] = NODE_X(v);
				}
				else
				{
					if (NODE_Y(v) < ylow[myid])
					{
						ylow[myid] = NODE_Y(v);
						xpos[myid] = NODE_X(v);
					}
				}
				k[myid]++;				
			}
		}
	}
	end_f_loop(f,tf)
	
	/* sync value to all comp nodes  */
	PRF_GRSUM(ylow, compute_node_count, dworkN);
	PRF_GRSUM(xpos, compute_node_count, dworkN);
	PRF_GISUM(k, compute_node_count, iworkN);

	i = 0;
	for (N = 0; N < compute_node_count; N++)	
	{
		/* Message("N:%i  ylow=%f, xpos=%f, k=%i\n", N, ylow[N], xpos[N], k[N]); */
		
		if (k[N] > 0)	/* this worker has nodes on the zone */
		{
			if (i==0)
			{
				ylowest = ylow[N];
				xapex = xpos[N];
				i++;
			}
			else
			{
				if (ylow[N] < ylowest)
				{
					ylowest = ylow[N];
					xapex = xpos[N];
				}
				i++;
			}
		}
	}
	/* Message("ylowest=%f, xapex=%f\n", ylowest, xapex); */
	
	free(ylow);
	free(xpos);
	free(dworkN);
	free(k);
	free(iworkN);
	#endif
	
	return xapex;
	
	
	/* OLD  */
	/* double xposition = 0.0; */
	
	/* #if !RP_HOST */
	/* face_t f; */
	/* int n; */
	/* Node *v; */
	
	
	/* begin_f_loop(f,tf)  */
	/* {   */
		/* if PRINCIPAL_FACE_P(f,tf) */
		/* { */
			/* f_node_loop(f,tf,n)  */
			/* { */
				/* v = F_NODE(f,tf,n); */
				/* if (NODE_Y(v) < 1E-10) */
					/* xposition = NODE_X(v); */
			/* } */
		/* } */
	/* } */
	/* end_f_loop(f,tf) */
	
	/* sync value to all comp nodes  */ 
	/* xposition = PRF_GRSUM1(xposition);  // = xposition + 0.0 + 0.0 ... + 0.0 */
	/* #endif */
	
	/* return xposition; */
	
}


/*---------- Get_NodeDistances ---------------------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
int Get_NodeDistances(Thread *tf, char meshZone, Node *holdNodes[], double xApex)
{
	int i = 0;
	
	#if !RP_HOST
		face_t f;
		Node *v;
		int n, j;
		int k = 0;
		double p[ND_ND]; 				/* midpoint from apex to tip  */
		double Xstar[ND_ND]; 			/* transformed node cartesian coordinates  */
		double Rstar[ND_ND];			/* transformed node polar coordinates  */
		
		/* Message("---------------\nGet_NodeDistances:\n"); */
		
		/* i. Draw line from apex to margin. Find midpoint.  */
		p[0] = 0.5*(xApex + xTipUnflexedOld);
		p[1] = 0.5*(0.0 + yTipUnflexedOld);		/* yApex is always 0.0 b/c axisymmetry  */
		
		Message0("\t\t p[0] = %f, p[1] = %f, xTipUnflexedOld = %f, yTipUnflexedOld = %f\n", p[0], p[1], xTipUnflexedOld, yTipUnflexedOld); 
		
		begin_f_loop(f,tf) {  
			if PRINCIPAL_FACE_P(f,tf) {
				f_node_loop(f,tf,n) {
					
					v = F_NODE(f,tf,n);
					
					if ((N_UDMI(v,0) == 0) || NodeIsTip(v)) {
						
						holdNodes[i] = v; 		/* store in an array of nodes  */	

						/* Get distances based on angle relative to the rough center of curvature  */
						/* ii. Define a new coordinate system (x*,y*) where the midpoint p is the origin. iii. Transform. */
						/* Xstar[0] = NODE_X(v) - p[0]; */
						/* Xstar[1] = NODE_Y(v) - p[1]; */
						Xstar[0] = N_UDMI(v,12) - p[0];
						Xstar[1] = N_UDMI(v,13) - p[1];
						
						
						/* 
							iv. Define a second transformation function g from cartesian to polar coordinates,
							i.e. such that g(x*,y*) -> R*(r,theta) where g(0,0) = (0,0).
						*/
						Rstar[0] = sqrt( SQR(Xstar[0]) + SQR(Xstar[1]) );		/*	radius, r */
						
						Rstar[1] = 0.0;											/* angle, theta */
						if ( (Xstar[1] >= 0.0) && (Rstar[0] != 0.0) )			/* if y>=0 and r!=0 */
							Rstar[1] = acos( Xstar[0]/Rstar[0] );
						
						else if ( Xstar[1] < 0.0 )								/* if y<0 */
							Rstar[1] = -acos( Xstar[0]/Rstar[0] );
						
						else if (Rstar[0] == 0.0)								/* if r=0 */
							Rstar[1] = 0.0;		/* undefined */

						else
							Error("\nError defining transformation for arc length calc \n");
						
						
						/* 
							v. Shift theta values to lie on [-2*pi,0] instead of [-pi, pi] (basically modulo).
								This gives the apex node the smallest value so that quickSort will correctly order
								holdNodes from apex to margin, i.e. smallest to largest theta "distance". 
						*/
						if (Rstar[1] < 0)
							Rstar[1] = -(Rstar[1] + 2*M_PI);
						else
							Rstar[1] = -Rstar[1];
						
						
						/* vi. Store "distance" (theta = radial distance) to use in sorting the nodes  */
						if (meshZone == 'e') 
						{
							N_UDMI(v,6) = Rstar[1]; /* Store initial distance to apex  */
							/* Message0("holdNodes[%i]: x = %lf, y = %lf, N_UDMI(v,6) = dist = %lf, NodeIsTip = %i \n",i, NODE_X(v), NODE_Y(v), N_UDMI(v,6), NodeIsTip(v));	 */
						}
						else if (meshZone == 's') 
						{
							N_UDMI(v,7) = Rstar[1]; /* Store initial distance to apex  */
							/* Message0("holdNodes[%i]:  x = %lf, y = %lf, N_UDMI(v,7) = dist = %lf, NodeIsTip = %i \n", i, NODE_X(v), NODE_Y(v), N_UDMI(v,7), NodeIsTip(v)); */
						}
						
						i++;
						N_UDMI(v,0) = 1;
					}
				}
			}
		}
		end_f_loop(f,tf)
	#endif
	Message0("\t\t i = %i\n", i);
	return i;
}


/**---------- Get_ParArcLengths -----------------------------
 * ! This is one of the bugged parallel funcs!
 * TODO: FIX THIS FUNCTION TO WORK WHEN A PROCESSOR HAS GROUPS OF NODES THAT AREN'T ALL NEXT TO EACH OTHER 
 *
---------------------------------------------------------------- */
void Get_ParArcLengths(char meshZone, Node *holdNodes[], int i, 
						double *SmaxUnflexed_ptr, int idArray[][2],
						double arclengthArray_unflex[] )
{
	#if RP_NODE
	
	/*--- 
	these arrays could be passed back as results but, if they aren't, they maybe should be
	allocated dynamically to avoid [nodes > MAX_NODES] out-of-bounds memory access errors 
	---	 */
	double coordArray_x_flex[MAX_NODES], coordArray_y_flex[MAX_NODES];
	double coordArray_x_unflex[MAX_NODES], coordArray_y_unflex[MAX_NODES];
	double arclengthArray_flex[MAX_NODES];
	/* double arclengthArray_unflex[MAX_NODES]; */

	int total_i, iHaveNodes = 0, totalHaveNodes;
	int *this_i, *iworkN, iworkj[MAX_NODES];
	double *dworkN, dworkj[MAX_NODES];
	int N, M, j, k = 0;
	int thisID;
	double thisArcDist_unflex, thisArcDist_flex;
	int rows, cols;

	int dist_mem_slot;
	int i_start, j_lcl;
	double masterArr_dist[MAX_NODES];
	double masterArr_x_un[MAX_NODES];
	double masterArr_y_un[MAX_NODES];
	double masterArr_x_fl[MAX_NODES];
	double masterArr_y_fl[MAX_NODES];
	double masterArr_myid[MAX_NODES];
	double masterArr_jlcl[MAX_NODES];
	double masterNodeIndex[7*MAX_NODES];

	if (meshZone == 'e')
	{
		dist_mem_slot = 6;
	}
	else
	{
		dist_mem_slot = 7;
	}

	/* calcaulte initial counters, misc variables  */
	/* Message("\tGet_ParArcLengths::calc initial counters, misc var\n"); */
	if (i > 0) 					/* if this processor's partition has any mesh nodes from this zone */
		iHaveNodes = 1;
	totalHaveNodes = PRF_GISUM1(iHaveNodes);	/* how many processors have mesh nodes from this zone */
	total_i = PRF_GISUM1(i);	/* total number of nodes documented for the mesh zone */

	/* dynamically allocate memory for arrays that need to be length of compute_node_count  */
	/* Message("\tGet_ParArcLengths::dynamically allocating memory\n"); */	
	this_i = (int *)malloc(sizeof(int)*compute_node_count);	/* allocate  */
	if (this_i == NULL)
		Error("malloc failed!\n");
	
	iworkN = (int *)malloc(sizeof(int)*compute_node_count);	/* allocate  */
	if (iworkN == NULL)
		Error("malloc failed!\n");
	
	dworkN = (double *)malloc(sizeof(double)*compute_node_count);	/* allocate  */
	if (dworkN == NULL)
		Error("malloc failed!\n");
	
	
	/* initialize arrays that contain one value per COMPUTATIONAL node */
	for (N = 0; N < compute_node_count; N++)
	{
		this_i[N] = 0;
		iworkN[N] = 0;
		dworkN[N] = 0.0;
	}
	
	/* sync an array that holds the # of mesh nodes contained on each processor */
	this_i[myid] = i;
	PRF_GISUM(this_i, compute_node_count, iworkN);	 
	Message0("\t\t I am Node%i and I have %i mesh nodes of mesh zone %c \n",myid, i, meshZone);
	
	/* initialize arrays that contain one value per MESH node */
	Message0("\t\t Initializing arrays\n");
	for (j = 0; j < MAX_NODES; j++)
	{
		idArray[j][0] = 0;
		idArray[j][1] = 0;
		coordArray_x_flex[j] = 0.0;
		coordArray_y_flex[j] = 0.0;
		coordArray_x_unflex[j] = 0.0;
		coordArray_y_unflex[j] = 0.0;
		arclengthArray_flex[j] = 0.0;
		arclengthArray_unflex[j] = 0.0;
		iworkj[j] = 0.0;
		dworkj[j] = 0.0;
		
		masterArr_dist[j] = 0.0;
		masterArr_x_un[j] = 0.0;
		masterArr_y_un[j] = 0.0;
		masterArr_x_fl[j] = 0.0;
		masterArr_y_fl[j] = 0.0;
		masterArr_myid[j] = 0;
		masterArr_jlcl[j] = 0;

		/* 2D initialization of all elements  */
		for (M = 0; M < 5; M++)
			masterNodeIndex[j*2 + M] = 0.0;
	}

	/**
	 * TODO: 1. populate masterNodeIndex w/ the following columns:
	 * TODO: dist, NODE_X, NODE_Y, myid, holdNodes index
	 * !CAN'T SYNC 2D ARRAYS or ARRAYS! Limitation of Fluent...
	 **/
	Message0("\t\t Populating masterArr vectors\n");
	i_start = 0;
	if (iHaveNodes)
	{
		for (N = 0; N < myid; N++)
		{
			i_start = i_start + this_i[N];
		}

		for (j = i_start; j < i_start + i; j++)
		{
			j_lcl = j - i_start;
			masterArr_dist[j] = N_UDMI(holdNodes[j_lcl], dist_mem_slot);
			masterArr_x_un[j] = N_UDMI(holdNodes[j_lcl], 12);
			masterArr_y_un[j] = N_UDMI(holdNodes[j_lcl], 13);
			masterArr_x_fl[j] = NODE_X(holdNodes[j_lcl]);
			masterArr_y_fl[j] = NODE_Y(holdNodes[j_lcl]);
			masterArr_myid[j] = (double) myid;
			masterArr_jlcl[j] = (double) j_lcl;

			if ( NodeIsTip(holdNodes[j_lcl]) )
			{
				/* Message0("\t\t Adding tip to masterArr! N_UDMI(v,6) = %lf, N_UDMI(v,7) = %lf\n",N_UDMI(holdNodes[j_lcl],6), N_UDMI(holdNodes[j_lcl],7)); */
				masterArr_x_un[j] = xTipUnflexedOld;
				masterArr_y_un[j] = yTipUnflexedOld;
				masterArr_x_fl[j] = xTipOld;
				masterArr_y_fl[j] = yTipOld;
			}
		}
	}
	Message0("\t\t Syncing masterArr vectors\n");
	PRF_GRSUM(masterArr_dist, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_x_un, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_y_un, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_x_fl, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_y_fl, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_myid, MAX_NODES, dworkj);
	PRF_GRSUM(masterArr_jlcl, MAX_NODES, dworkj);

	/* assemble vectors into one master array for sorting */
	Message0("\t\t Assembling masterArr \n");
	cols = 7;
	for (j = 0; j < total_i; j++)
	{	
		masterNodeIndex[j*cols + 0] = masterArr_dist[j];
		masterNodeIndex[j*cols + 1] = masterArr_x_un[j];
		masterNodeIndex[j*cols + 2] = masterArr_y_un[j];
		masterNodeIndex[j*cols + 3] = masterArr_x_fl[j];
		masterNodeIndex[j*cols + 4] = masterArr_y_fl[j];
		masterNodeIndex[j*cols + 5] = masterArr_myid[j];
		masterNodeIndex[j*cols + 6] = masterArr_jlcl[j];
	}

	/**
	 * TODO: 2. sort masterNodeIndex by dist using 2D bubble sort and reorder the vectors
	 **/
	/* sort the 2D array so the order goes from apex mesh nodes to tip mesh nodes 
			(i.e. ascending radial distance)  */
	Message0("\t\t Sorting masterArr \n");
	bubbleSort_double2D(masterNodeIndex, total_i, cols, 'a'); 

	Message0("\t\t Breaking masterArr back into vectors \n");
	for (j = 0; j < total_i; j++)
	{	
		masterArr_dist[j] = masterNodeIndex[j*cols + 0];
		masterArr_x_un[j] = masterNodeIndex[j*cols + 1];
		masterArr_y_un[j] = masterNodeIndex[j*cols + 2];
		masterArr_x_fl[j] = masterNodeIndex[j*cols + 3];
		masterArr_y_fl[j] = masterNodeIndex[j*cols + 4];
		masterArr_myid[j] = masterNodeIndex[j*cols + 5];
		masterArr_jlcl[j] = masterNodeIndex[j*cols + 6];

		idArray[j][0] = (int) masterArr_myid[j];
		idArray[j][1] = (int) masterArr_jlcl[j];
	}

	/**
	 * TODO: 3. Calc flexed and unflexed arclengths
	 * ? Are the old outputs actually used?
	 * ? *Smax_ptr and Smax -- NO -> deleted
	 * ? *SmaxUnflexed_ptr and SmaxUnflexed -- YES
	 * ? idArray[][2] -- YES
	 * ? coordArray_x_flex[] -- Load_b_sub? NOT USED external to this function --> del
	 * ? coordArray_y_flex[] -- Load_a_sub? NOT USED external to this function --> del
	 * ? arclengthArray_unflex[] -- YES (by Load_a/b_sub)
	 * ? arclengthArray_flex[]-->N_UDMI(v, 1) -- NO? But keeping it for now for simplicity
	 **/
	Message0("\t\t Calculating arc lengths \n");
	for (j = 1; j < total_i; j++)
	{	
		/* Unflexed  */
		thisArcDist_unflex = sqrt( SQR(masterArr_x_un[j] - masterArr_x_un[j-1]) + 
								SQR(masterArr_y_un[j] - masterArr_y_un[j-1]) );
		arclengthArray_unflex[j] = arclengthArray_unflex[j-1] + thisArcDist_unflex;

		/* Flexed  */
		thisArcDist_flex = sqrt( SQR(masterArr_x_fl[j] - masterArr_x_fl[j-1]) + 
								SQR(masterArr_y_fl[j] - masterArr_y_fl[j-1]) );
		arclengthArray_flex[j] = arclengthArray_flex[j-1] + thisArcDist_flex;
	}

	Message0("\t\t Printing out idArray:\n");
	for (j = 0; j < total_i; j++)
	{	
		thisID = (int) masterArr_myid[j];
		j_lcl = (int) masterArr_jlcl[j];

		/* Message0("\t\t j:%i  idArray = [%i][%i]  dist = %lf\n", j, thisID, j_lcl, masterArr_dist[j]); */
	}

	Message0("\t\t Saving locally \n");
	/* save locally */
	for (j = 0; j < total_i; j++)
	{	
		thisID = (int) masterArr_myid[j];
		j_lcl =  (int) masterArr_jlcl[j];
		if (myid == thisID)
		{
			N_UDMI(holdNodes[j_lcl], 1) = arclengthArray_flex[j];
			N_UDMI(holdNodes[j_lcl], 14) = arclengthArray_unflex[j];
		}
	}


	/* Debug */
	for (j = 0; j < total_i; j++)
	{
		/* Message0("\t\t arclengthArray_unflex[%i] = %10.10f (dist = %lf, x = %lf, y = %lf)\n", j, arclengthArray_unflex[j], masterArr_dist[j], masterArr_x_un[j], masterArr_y_un[j] );  */
	}
	
	/**
	 * TODO: 4. Calc SmaxUnflexed
	 **/

	/* Last index in global array is the tip node. Get SmaxUnflexed  */
	*SmaxUnflexed_ptr = arclengthArray_unflex[total_i - 1];
	Message0("\t\t Smaxunflexed = %lf\n", *SmaxUnflexed_ptr);

	/* deallocate dynamic memory  */
	free(this_i);
	free(iworkN);
	free(dworkN);

	
	#endif
}


/*---------- NodeIsTip -------------------------------------------
Purpose: 

Input:	FILE *stream	-	File stream
---------------------------------------------------------------- */
int NodeIsTip (Node *v) 
{
	#if !RP_HOST
		/* double x = NODE_X(v); */
		/* double y = NODE_Y(v); */
		double xold = N_UDMI(v,2);
		double yold = N_UDMI(v,3);
		double r1 = sqrt(SQR(xold - xTipOld) + SQR(yold - yTipOld));
		/* double r2 = sqrt(SQR(N_UDMI(v,12) - xTipUnflexedOld) + SQR(N_UDMI(v,13) - yTipUnflexedOld)); */
		
		/* if (N_TIME == 0) { */
			if (r1 < 1e-7 ) 
			{ 
				
				/*Message0("\tTIP FOUND @ xold = %f, yold = %f /* xnew = %f, ynew = %f\n", 
						xold, yold, NODE_X(v), NODE_Y(v));*/
				
				/* N_UDMI(v,15) = -1.0; */
				return 1;
			}
			else
				return 0;
		/* } */
		/* else { */
			/* if (N_UDMI(v,15) == -1.0) */
				/* return 1; */
			/* else */
				/* return 0; */
		/* } */
	
	#endif
	
	#if RP_HOST
		return 0;
	#endif
}



/*---------- Get_tparm -------------------------------------------
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
void Get_tparm(char meshZone, Node *holdNodes[],  double b, double b_sub[], int nNodes, 
				double SmaxUnflexed, int idArray[][2]) 
{
	#if !RP_HOST
	int j;
	int oldNode, lastOldNode, nextOldNode;
	int iL, iR;
	int iL_is_oldNode, iR_is_oldNode;
	Node *v;
	double Si, this_input = 0;
	double tparm_old_j, tparm_old_iL, tparm_old_iR;
	double SunflexedL, SunflexedR, Sunflexedj;
	
	int oldNode_Arr[MAX_NODES];
	double tparm_old_Arr[MAX_NODES], Sunflexed_Arr[MAX_NODES];
	int iworkj[MAX_NODES];
	double dworkj[MAX_NODES];
	int thisID, localj;
	int loop;
	double dif;
	
	
	/* initialize tparm arrays  */
	if (N_TIME == 0) {				/* also satisfies NewTimeStepForThisZone() b/c function call is Calc_Kinematics_and_Move in Calc_Mesh_Movement() */
		if (meshZone == 'e') {
			for (j = 0; j < MAX_NODES; j++) {
				tparm_ex[j] = 0.0;
				tparm_ex_old[j] = 0.0;
			}
		}
		else  /* (meshZone == 's') */
			for (j = 0; j < MAX_NODES; j++) {
				tparm_sub[j] = 0.0;
				tparm_sub_old[j] = 0.0;
			}
	}
	
	/* initialize other arrays  */
	for (j = 0; j < MAX_NODES; j++)
	{
		oldNode_Arr[j] = 0;
		tparm_old_Arr[j] = 0.0;
		Sunflexed_Arr[j] = 0.0;
	}
	
		
	/** Exumbrella * */
    if (meshZone == 'e') 
	{
		/* j = 0: set boundary condition at apex  */
		j = 0;
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
		
		if (myid == thisID)
		{
			v = holdNodes[localj];
			tparm_ex[j] = M_PI;
			N_UDMI(v,10) = M_PI;			/* tparm */
			N_UDMI(v,8) = N_UDMI(v,10);		/* tparm_old */
			N_UDMI(v,9) = 1;				/* old node flag (I guess assume that it is never remeshed?) */
			
			/** diagnostics * */
			Message0("\tj:%i  tp = %f, tpo = %f, x = %lf, y = %lf\n", 0, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v));
			/* Message("j:%i  x = %lf, x-b = %lf, b = %lf\n",0, NODE_X(holdNodes[0]), NODE_X(holdNodes[0])-b, b); */
		}
		
		
		/* 0 < j < nNodes: calculate values for all other nodes  */
		if (N_TIME==0) {	/* if first time step... */
			for (j = 1; j < nNodes; j++) 
			{	
				thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
				localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */
				
				if (myid == thisID)
				{
					v = holdNodes[localj];
					tparm_ex[j] = acos((NODE_X(v)-b)/b);	/* store in global array */
					N_UDMI(v,10) = tparm_ex[j];			/* store tparm (moving code away from global arrays like tparm_ex[]) */
					N_UDMI(v,8) = N_UDMI(v,10);		/* set tparm_old = tparm (execute_at_end macro doesn't work at N_TIME = 0) */
					
					/** diagnostics * */
					Message0("\t\tj:%i  tp = %f, tpo = %f, x = %lf, y = %lf\n", j, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v));
				}
			}
		}
		else {		/* if after first time step... */
			
			/* 1. Populate arrays of local node info  */
			for (j = 1; j < nNodes; j++) 
			{
				thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
				localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */
				
				if (myid == thisID)
				{
					v = holdNodes[localj];
					oldNode_Arr[j] = N_UDMI(v,9);		/* flag for old (i.e. not remeshed) node */
					tparm_old_Arr[j] = N_UDMI(v,8);
					Sunflexed_Arr[j] = N_UDMI(v,14);
				}
			}
			
			
			/* 2. Sync local info to create global arrays  */
			PRF_GISUM(oldNode_Arr, MAX_NODES, iworkj);
			PRF_GRSUM(tparm_old_Arr, MAX_NODES, dworkj);
			PRF_GRSUM(Sunflexed_Arr, MAX_NODES, dworkj);
			
			
			/* 3. Loop again to calc tparm  */
			for (j = 1; j < nNodes; j++) 
			{
				thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
				localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
				
				if (myid == thisID)
				{
					v = holdNodes[localj];
					Si = N_UDMI(v, 14);
					
					/* 3a. Handle new nodes added via remeshing/adaption  */
					if (!oldNode_Arr[j]) {	
						Message0("\tNew ex node detected\n");
						newNodeFlag = 1;
						
						/* new node added from remeshing, thus has no tparm_sub_old value! Must estimate it.  */
						lastOldNode = -1;
						nextOldNode = -1;
					
						/* scan left for last oldNode, get its node number  */
						iL = j;
						while ( (lastOldNode == -1) && (iL >= 0) )
						{
							iL--;
							iL_is_oldNode = oldNode_Arr[iL];	/* = N_UDMI(holdNodes[iL],9) for global/serial */
							if (iL_is_oldNode)
								lastOldNode = iL;
						}
						
						/* scan right for next oldNode, get its node number  */
						iR = j;
						while ( (nextOldNode == -1) && (iR < nNodes) )
						{
							iR++;
							iR_is_oldNode = oldNode_Arr[iR];	/* = N_UDMI(holdNodes[iR],9) for global/serial */
							if (iR_is_oldNode)
								nextOldNode = iR;
						}
						
						/* error check that lastOldNode and nextOldNode were found or j is a boundary node  */
						if ((lastOldNode == -1) || (nextOldNode == -1))
						{
							if (j == nNodes-1)  /* tip node */
								N_UDMI(v,8) = tparmTipOld;				/* set it to stored old tip value */
							else
								Error("Could not interpolate tparm_ex_old value of new node! Left oldNode/Right oldNode not found.");
						}
						
						/* interpolate tparmOld for the new node using lastOldNode and nextOldNode  */
						/* tparm_old_iL = tparm_old_Arr[iL];		// = N_UDMI(holdNodes[iL],8) */
						/* tparm_old_iR = tparm_old_Arr[iR];		// = N_UDMI(holdNodes[iR],8) */
						/* SunflexedL = Sunflexed_Arr[iL];		// = N_UDMI(holdNodes[iL],14) */
						/* SunflexedR = Sunflexed_Arr[iR];		// = N_UDMI(holdNodes[iR],14) */
						
						tparm_old_Arr[j] = tparm_old_Arr[iL] + (double)(j-iL)/(iR-iL) * (tparm_old_Arr[iR] - tparm_old_Arr[iL]);
						N_UDMI(v,8) = tparm_old_Arr[j];
						Message0("\tj: %i	iL: %i	iR: %i	tparm_old_iL: %lf	tparm_old_iR: %lf \n", j, iL, iR, tparm_old_Arr[iL], tparm_old_Arr[iR]);
						
					}
					
					/* calculate the new tparm value  */
					if (N_TIME == 1)
					{
						tparm_ex[j] = tparm_old_Arr[j] + Si/SmaxUnflexed * (tparmTip - tparmTipOld); 
						N_UDMI(v,15) = Si/SmaxUnflexed;
					}
					else
					{
						tparm_ex[j] = tparm_old_Arr[j] + N_UDMI(v,15) * (tparmTip - tparmTipOld); 
					}
					/* tparm_ex[j] = tparm_old_Arr[j] + (double)j/(nNodes-1)*(tparmTip - tparmTipOld);  */
					tparm_ex[j] = tparm_old_Arr[j] + Si/SmaxUnflexed * (tparmTip - tparmTipOld); 
					N_UDMI(v,10) = tparm_ex[j];
					
					/* old codes */
					/* tparm_ex[j] = tparm_old_j + Sunflexedj/SmaxUnflexed*(tparmTip - tparmTipOld);  */
					/* tparm_ex[j] = tparm_ex_old[j] + (double)j/(nNodes-1)*(tparmTip - tparmTipOld); */
					
					/* diagnostics  */
					/* Message("j:%i  x = %lf, x-b = %lf, b = %lf\n",j, NODE_X(v), NODE_X(v)-b, b); */
					/* Message("j:%i  Snorm = %lf, tp = %f, tpo = %lf, x = %lf, y = %lf\n",j, Sunflexedj/SmaxUnflexed, tparm_ex[j], tparm_old_j, NODE_X(v), NODE_Y(v)); */
					/**DEBUG* */
					/* if (N_TIME >= 63) */
						/* Message("j:%i  tp = %f, tpo = %lf, x = %lf, y = %lf\n",j, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v)); */
				}
			}
		}
	}
	

	/** Subumbrella * */
	else if (meshZone == 's') 
	{
		/* j = 0: set boundary condition at apex  */
		j = 0;
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
		
		if (myid == idArray[0][0])
		{
			v = holdNodes[localj];
			tparm_sub[j] = M_PI;
			N_UDMI(v, 10) = M_PI;			/* tparm */
			N_UDMI(v, 8) = N_UDMI(v,10);	/* tparm_old */
			N_UDMI(v, 9) = 1;				/* old node flag */
			
			/* diagnostics  */
			Message0("\tj:%i  tp = %f, tpo = %f, x = %lf, y = %lf\n", 0, N_UDMI(holdNodes[0],10), N_UDMI(holdNodes[0],8), NODE_X(holdNodes[0]), NODE_Y(holdNodes[0]));
			/* Message("j:%i  x = %lf, x-b = %lf, b = %lf\n",0, NODE_X(holdNodes[0]), NODE_X(holdNodes[0])-b, b); */
		}
		
		/* 0 < j < nNodes: calculate values for all other nodes   */
		if (N_TIME==0) {	/* if first time step... */
			for (loop = 0; loop < 2; loop++) {
				for (j = nNodes-1; j > 0; j--) 
				{
					thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
					localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */
					
					if (myid == thisID) 
					{
						v = holdNodes[localj];
						
						/* tparm calc LOOP 1 -- check b_sub  */
						if (loop == 0)
						{
							/* this_input = (NODE_X(v)-b)/b_sub[j]; */
							/* if (j ==1) Message("this_input = %lf, acos(this_input) = %lf, NODE_X = %lf, b = %lf, b_sub = %lf\n", this_input, acos(this_input), NODE_X(v), b, b_sub[j]); */
							/* tparm_sub[j] = acos(this_input); */
							tparm_sub[j] = tparmTip + (b_sub[j]-b)/(b-thickness - b) * (M_PI - tparmTip);
							
							/* Check b_sub  */
							if (this_input < -1.0)
							{
								/* new method  */
								/* interp b_sub method  */
								/* Message("Invalid b_sub value at j = %i. Interpolating b_sub from neighbors.\n",j); */
								/* b_sub[j] = (b - thickness) + (double)j/(j+1) * (b_sub[j+1] - (b - thickness)); */
								/* this_input = (NODE_X(v)-b)/b_sub[j]; */
								/* tparm_sub[j] = acos(this_input); */
								
								/* old method  */
								/* interp tparm method (NEEDS TO UPDATE B_SUB TOO!)  */
								Message0("Invalid b_sub value at j = %i. Interpolating tparm_sub from neighbors.\n",j);
								tparm_sub[j] = M_PI + (double)j/(j+1) * (tparm_sub[j+1]-M_PI);	/* assuming this can only happen near j=0 (i.e. apex node) */
								/* dif = (tparm_sub[j]-tparmTip)/(tparm_sub[0]-tparmTip) * thickness; */
								/* b_sub[j] = b - dif; */
								
							}
							
							/* store tparm  */
							N_UDMI(v,10) = tparm_sub[j];					/* store tparm (moving code away from global arrays like tparm_ex[]) */
							N_UDMI(v,8) = N_UDMI(v,10);		/* set tparm_old = tparm (execute_at_end macro doesn't work at N_TIME = 0) */
							
							/* diagnostics  */
							Message0("\t\tj:%i  tp = %f, tpo = %f, x = %lf, y = %lf\n", j, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v));
						}
						
						/* tparm calc LOOP 2 -- check tparm  */
						else if (loop == 1)
						{
							/* Check tparm  */
							/* if (tparm_sub[j] > tparm_sub[j-1]) */
							/* { */
								/* Message("Invalid tparm value = %f at j = %i. Interpolating tparm_sub from neighbors.\n",tparm_sub[j],j); */
								/* tparm_sub[j] = M_PI + (double)j/(j+1) * (tparm_sub[j+1]-M_PI);	// assuming this can only happen near j=0 (i.e. apex node) */
								/* N_UDMI(v,10) = tparm_sub[j];					// store tparm (moving code away from global arrays like tparm_ex[]) */
								/* N_UDMI(v,8) = N_UDMI(v,10);		// set tparm_old = tparm (execute_at_end macro doesn't work at N_TIME = 0) */
							
								/* diagnostics  */ 
								/* Message("\t\tj:%i  tp = %f, tpo = %f, x = %lf, y = %lf\n", */
									/* j, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v)); */
							/* } */
						}
					}
				}
			}
		}
		else {		/* if after first time step... */
			
			/* 1. Populate arrays of local node info  */
			for (j = 1; j < nNodes; j++) 
			{
				thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
				localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */
				
				if (myid == thisID)
				{
					v = holdNodes[localj];
					oldNode_Arr[j] = N_UDMI(v,9);						/* flag for old (i.e. not remeshed) node */
					tparm_old_Arr[j] = N_UDMI(v,8);
					Sunflexed_Arr[j] = N_UDMI(v,14);
				}
			}
			
			
			/* 2. Sync local info to create global arrays  */
			PRF_GISUM(oldNode_Arr, MAX_NODES, iworkj);
			PRF_GRSUM(tparm_old_Arr, MAX_NODES, dworkj);
			PRF_GRSUM(Sunflexed_Arr, MAX_NODES, dworkj);
			
			
			/* 3. Loop again to calc tparm  */
			for (j = 1; j < nNodes; j++) 
			{
				thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
				localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
				
				if (myid == thisID)
				{
					v = holdNodes[localj];
					Si = N_UDMI(v, 14);
					
					/* 3a. Handle new nodes added via remeshing/adaption  */
					if (!oldNode_Arr[j]) {
						Message0("New sub node detected\n");
						newNodeFlag = 1;
						
						/* new node added from remeshing, thus has no tparm_sub_old value! Must estimate it.  */
						lastOldNode = -1;
						nextOldNode = -1;
					
						/* scan left for last oldNode, get its node number  */
						iL = j;
						while ( (lastOldNode == -1) && (iL >= 0) )
						{
							iL--;
							iL_is_oldNode = oldNode_Arr[iL];	/* = N_UDMI(holdNodes[iL],9) for global/serial */
							if (iL_is_oldNode)
								lastOldNode = iL;
						}
						
						/* scan right for next oldNode, get its node number  */
						iR = j;
						while ( (nextOldNode == -1) && (iR < nNodes) )
						{
							iR++;
							iR_is_oldNode = oldNode_Arr[iR];	/* = N_UDMI(holdNodes[iR],9) for global/serial */
							if (iR_is_oldNode)
								nextOldNode = iR;
						}
						
						/* error check that lastOldNode and nextOldNode were found or j is a boundary node  */
						if ((lastOldNode == -1) || (nextOldNode == -1))
						{
							if (j == nNodes-1)
								N_UDMI(v,8) = tparmTipOld;				/* set it to stored old tip value */
							else
								Error("Could not interpolate tparm_sub_old value of new node! Left oldNode/Right oldNode not found.");
						}
						
						/* interpolate tparmOld for the new node using lastOldNode and nextOldNode  */
						tparm_old_iL = tparm_old_Arr[iL];		/* = N_UDMI(holdNodes[iL],8); */
						tparm_old_iR = tparm_old_Arr[iR];		/* = N_UDMI(holdNodes[iR],8); */
						SunflexedL = Sunflexed_Arr[iL];		/* = N_UDMI(holdNodes[iL],14); */
						SunflexedR = Sunflexed_Arr[iR];		/* = N_UDMI(holdNodes[iR],14); */
						
						tparm_old_Arr[j] = tparm_old_Arr[iL] + (double)(j-iL)/(iR-iL) * (tparm_old_Arr[iR] - tparm_old_Arr[iL]);
						N_UDMI(v,8) = tparm_old_Arr[j];
						Message0("\tj: %i	iL: %i	iR: %i	tparm_old_iL: %lf	tparm_old_iR: %lf \n", j, iL, iR, tparm_old_Arr[iL], tparm_old_Arr[iR]);
					
					}
					
					/* calculate the new tparm value  */
					if (N_TIME == 1)
					{
						tparm_sub[j] = tparm_old_Arr[j] + Si/SmaxUnflexed * (tparmTip - tparmTipOld); 
						N_UDMI(v,15) = Si/SmaxUnflexed;
					}
					else
					{
						tparm_sub[j] = tparm_old_Arr[j] + N_UDMI(v,15) * (tparmTip - tparmTipOld); 
					}
					/* tparm_sub[j] = tparm_old_Arr[j] + (double)j/(nNodes-1)*(tparmTip - tparmTipOld); */
					N_UDMI(v,10) = tparm_sub[j];
					
					/* old codes  */
					/* tparm_sub[j] = tparm_old_j + Sunflexedj/SmaxUnflexed*(tparmTip - tparmTipOld);  */
					/* tparm_sub[j] = tparm_sub_old[j] + (double)j/(nNodes-1)*(tparmTip - tparmTipOld); */
					
					/* diagnostics  */
					/* Message("j:%i  x = %lf, x-b = %lf, b = %lf\n",j, NODE_X(v), NODE_X(v)-b, b); */
					/* Message("j:%i  Snorm = %lf, tp = %f, tpo = %lf, x = %lf, y = %lf\n",j, Sunflexedj/SmaxUnflexed, tparm_sub[j], tparm_old_j, NODE_X(v), NODE_Y(v)); */
					/**DEBUG* */
					/* if (N_TIME >= 63) */
						/* Message("j:%i  tp = %f, tpo = %lf, x = %lf, y = %lf\n",j, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v));  */
				}
			}
		}
	}
	#endif
}


/*---------- get_a_sub ------------------------------------
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
void Get_a_sub(double *a_sub_ptr, double a, Node *holdNodes[], int nNodes, 
				double SmaxUnflexed, int idArray[][2]) 
{
	#if !RP_HOST
	int j;
	double Si, dif;
	Node *v;
	int thisID, localj;
	
	/* Loop over nodes, computing a_sub using change in tparm  */
	for (j = 0; j < nNodes; j++) 
	{
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
		
		if (myid == thisID)
		{
		/* THIS SECTION NOT USED ----------- */
		/*	v = holdNodes[localj];	*/			/* get node from ordered array */
			
			/* if (NodeIsTip(v))	*/		/* if tip node */ 
				/* Si = SmaxUnflexed; */					/* its arc length is the max arc length */ 
			/* else */
				/* Si = N_UDMI(v,14); */
		/* ---------------------------- */
			
			dif = (tparm_sub[j]-tparmTip)/(M_PI-tparmTip) * thickness;
			a_sub_ptr[j] = a - dif;
			
			/* if ((myid == 0) || (myid == 1)) */
			/* { */
				/* Message("j:%i  tparm_sub=%f, tparmTip=%f, tparm_sub[0]=%f, dif=%f\n", */
					/* j, tparm_sub[j], tparmTip, tparm_sub[0], dif); */
			/* } */
			/* Message("j:%i  j/(i-1) = %lf, pow(j/(i-1) = %lf, thickness = %lf\n",j, (double)j/(i-1), pow((double)j/(i-1), N),thickness ); */
			/* Message("j:%i  x = %lf, y = %lf, a_sub = %lf\n",j, NODE_X(v), NODE_Y(v), a_sub_ptr[j]); */
			/* Message("j:%i  x = %lf, y = %lf, a = %lf, a_sub = %lf, dif = %lf\n\n",j, NODE_X(v), NODE_Y(v), a, a_sub_ptr[j], dif); */
		}
	}
	#endif
}


/*---------- get_b_sub ------------------------------------
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
void Get_b_sub(double *b_sub_ptr, double b, Node *holdNodes[], int nNodes, 
				double SmaxUnflexed, int idArray[][2])  
{
	#if !RP_HOST
	int j;
	double Si, dif;
	Node *v;
	int thisID, localj;
	
	/* Loop over nodes, computing b_sub using change in tparm  */
	for (j = 0; j < nNodes; j++) 
	{
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];		/* local holdNodes index for node j on comp node thisID */
		
		if (myid == thisID)
		{
		/* THIS SECTION NOT USED -----------
			v = holdNodes[localj];
			if (NodeIsTip(v))
				Si = SmaxUnflexed;
			else
				Si = N_UDMI(v,14);
		---------------------------- */
			
			dif = (tparm_sub[j]-tparmTip)/(M_PI-tparmTip) * thickness;
			b_sub_ptr[j] = b - dif;
			
			/* if ((myid == 0) || (myid == 1)) */
			/* { */
				/* Message("j:%i  tparm_sub=%f, tparmTip=%f, tparm_sub[0]=%f, dif=%f\n", */
					/* j, tparm_sub[j], tparmTip, tparm_sub[0], dif); */
			/* } */
			/* Message("j:%i  j/(i-1) = %lf, pow(j/(i-1) = %lf, thickness = %lf\n",j, (double)j/(i-1), pow((double)j/(i-1), N),thickness ); */
			/* Message("j:%i  x = %lf, y = %lf, b = %lf, b_sub = %lf, dif = %lf\n",j, NODE_X(v), NODE_Y(v), b, b_sub_ptr[j], dif); */
		}
	}
	#endif
}


/*---------- Get_FlexPoint -------------------------------------
Purpose: Returns the tparm value for both the ex and sub surfaces

Input:	holdNodes			-	array of sorted nodes, from apex to margin
		i 					-	total number of nodes in this zone
		SmaxUnflexed		-	maximum arc length of this zone	BEFORE flexing	
		
Output:		i1		-	index of holdNodes where flex starts
---------------------------------------------------------------- */
void Get_FlexPoint(char meshZone, Node *holdNodes[], int nNodes, double SmaxUnflexed, 
					double sg, double g[], int idArray[][2]) 
{
	#if !RP_HOST
	int i1 = nNodes-1;
	int j;
	int lastOldNode, nextOldNode, iL, iR, iL_is_oldNode, iR_is_oldNode, oldNode;
	Node *v;
	double Si;
	double Sflex = pg*SmaxUnflexed;
	double flexpct_iL, flexpct_iR, flexpct_j, flexspacing_iL, flexspacing_iR, flexspacing_j;
	int thisID, localj;
	double dworkj[MAX_NODES];
	
	/* First loop: find first node that is flexing (i.e. node to the right of flex point).  */
	for (j = 0; j < nNodes; j++) 
	{ 		
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */

		if (myid == thisID)
		{
			v = holdNodes[localj];
			if (NodeIsTip(v)) {
				Si = SmaxUnflexed;
				N_UDMI(v,11) = 1.0;		/* percent flex */
			}
			else {
				Si = N_UDMI(v,14);
				N_UDMI(v,11) = 0.0;		/* percent flex (initialize) */
			}

			if (Si >= Sflex) {
				N_UDMI(v,11) = pow((Si - Sflex)/(SmaxUnflexed - Sflex), k);
				if (j < i1)
					i1 = j;
			}
			/* // Message("j:%i  Si = %lf, Sflex = %lf, SmaxUnflexed = %lf, i1 = %i\n", j, Si, Sflex, SmaxUnflexed, i1); */
		}
	}
	
	
	/* Second loop: Assign flex values.  */
	for (j = 0; j < nNodes; j++) 
	{
		thisID = idArray[j][0];		/* the comp node ID corresponding to node j in global holdNodes */
		localj = idArray[j][1];	/* local holdNodes index for node j on comp node thisID */

		if (myid == thisID)
		{
			v = holdNodes[localj];
			flexpct_j = N_UDMI(v,11);
			flexspacing_j = pow(flexpct_j, 1.0/k);
			
			g[j] = g1 + flexpct_j * (g2-g1) *sg;	
		}
	}
	
	/** DEBUG * */
	/* if ( (N_TIME > 63) || (N_TIME == 1) ) */
	/* { */
		/* for (j = 0; j < nNodes; j++)  Message("(before) j:%i  g[j] = %lf\n", j, g[j]); */
	/* } */
	
	
	
	/* -- sync local g arrays to make global array --  */
	#if PARALLEL
	/* Message("\n*Sync*\n"); */
	PRF_GRSUM(g, MAX_NODES, dworkj); 
	#endif
	
	
	
	
	/** DEBUG * */
	/* if ( (N_TIME > 63) || (N_TIME == 1) ) */
	/* { */
		/* for (j = 0; j < nNodes; j++)  Message("(after) j:%i  g[j] = %lf\n", j, g[j]); */
	/* } */

	/*  -- loop over global array and set g[j] = 0.0 to 1.0 --  */ 
	/* for (j = 0; j < nNodes; j++)  */
	/* { */
		/* thisID = idArray[j][0];		// the comp node ID corresponding to node j in global holdNodes */
		/* localj = idArray[j][1];	// local holdNodes index for node j on comp node thisID */

		/* if (myid == thisID) */
		/* { */
			/* if (g[j] == 0.0) */
				/* g[j] = 1.0; */
		/* } */
		
		/** DEBUG **/ 
		/* if ( N_TIME >= 63 ) */
			/* Message("j:%i  g[j] = %lf\n", j, g[j]);	 */
	/* } */
	

	#endif
}





/**---------- Store_OldKinematicVars ---------------------------
* TODO 1: Clean this up and make sure this isn't part of the problem w/ running on n>2 nodes
* TODO 2: Store the old nodal coordinates here in N_UDMI(v,2) and N_UDMI(v,3) 
---------------------------------------------------------------- */
void Store_OldKinematicVars(void)
{
	#if !RP_HOST
	Thread *tf_ex = Lookup_Thread(domain, EX_ZONE);
	Thread *tf_sub = Lookup_Thread(domain, SUB_ZONE);
	Thread *tf_array[NUM_ZONES];
	Thread *tf;
	face_t f;
	Node* v;
	int i, n, N, tipScan;
	int *tip_arr, *iworkN;
	
	#if PARALLEL
	
	/** Simple tip coord sync **/
	/* ----------------------  */
	/* store tip coords  */
	/* ----------------------  */
	/* (Parallel) sync xTip_temp, yTip_temp, etc.  */
	Message0("Syncing tip coordinates...\n");
	/* Message("(before) tip = (%f, %f) tipUnflexed = (%f, %f)\n", xTip, yTip, xTipUnflexed, yTipUnflexed); */
	/* Message("(before) tip_temp = (%f, %f) tipUnflexed_temp = (%f, %f)\n", xTip_temp, yTip_temp, xTipUnflexed_temp, yTipUnflexed_temp); */
	if (PRF_GRSUM1(iHaveTip) == 1)
	{
		xTip = PRF_GRSUM1(xTip_temp); 	/* = xTip_temp + 0.0 + 0.0 ... + 0.0 */
		yTip = PRF_GRSUM1(yTip_temp); 	/* = yTip_temp + 0.0 + 0.0 ... + 0.0 */
		xTipUnflexed = PRF_GRSUM1(xTipUnflexed_temp); 	/* = xTipUnflexed_temp + 0.0 + 0.0 ... + 0.0 */
		yTipUnflexed = PRF_GRSUM1(yTipUnflexed_temp); 	/* = yTipUnflexed_temp + 0.0 + 0.0 ... + 0.0 */
	}
	else
		Error("Tip either found on >1 node OR it was not found not at all! (iHaveTip = %i, PRF_GRSUM1(iHaveTip) = %i)", iHaveTip, PRF_GRSUM1(iHaveTip));
	
	iHaveTip = 0;
	
	#endif
	
	/* set local old = local current/new  */
	/** NOTE: Nomenclature is clearly redundant now. The "current/new" variables are never used
		other than to update the "old" variables, here, at the end of a time step. 
		Proposed variable changes:
		xTip_temp --> xTip_new	(will be = 0.0 unless tip node exists on comp node)
		xTipOld --> xTip	
	**/
	xTipOld = xTip;
	yTipOld = yTip;
	xTipUnflexedOld = xTipUnflexed;
	yTipUnflexedOld = yTipUnflexed;
	
	/* Message("(after) tip = (%f, %f) tipUnflexed = (%f, %f)\n", xTip, yTip, xTipUnflexed, yTipUnflexed); */
	/* Message("(after) tip_temp = (%f, %f) tipUnflexed_temp = (%f, %f)\n", xTip_temp, yTip_temp, xTipUnflexed_temp, yTipUnflexed_temp); */
	
	/* reset local temp  */
	xTip_temp = 0.0;
	yTip_temp = 0.0;
	xTipUnflexed_temp = 0.0;
	yTipUnflexed_temp = 0.0;
	
	
	/* ----------------------  */
	/* store tparm tip  */
	/* ----------------------  */
	tparmTipOld = tparmTip;
	
	
	/* ----------------------  */
	/* store tparm for all other nodes  */
	/* ----------------------  */
	tf_array[0] = tf_ex;
	tf_array[1] = tf_sub;
	
	Message0("Storing tparm_old values on surface ");
	for (i = 0; i < NUM_ZONES; i++)
	{
		Message0("%i... ", i);
		tf = tf_array[i];
		
		begin_f_loop(f, tf) {  
			f_node_loop(f, tf, n) 
			{
				v = F_NODE(f, tf, n);
				/* Message("end j:%i  tp = %f, tpo = %lf, x = %lf, y = %lf\n",v, N_UDMI(v,10), N_UDMI(v,8), NODE_X(v), NODE_Y(v)); */
				N_UDMI(v,8) = N_UDMI(v,10);			/* tparm_old = tparm				 */
			}
		}
		end_f_loop(f, tf)
	}
	
	Message0("Done.\n");
	#endif
}


