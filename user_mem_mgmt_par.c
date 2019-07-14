#include "user_mem_mgmt_par.h"

/** global variables **/


/** Function definitions **/
/*---------- Init_node_mem ------------------------------------
Purpose: 

Input:	Thread *tf - face thread where node memory is being accessed
----------------------------------------------------------------*/
int Init_node_mem(Thread *tf) 
{
	#if !RP_HOST
	face_t f;
	Node *v;
	int n;
	
	begin_f_loop(f,tf) 
	{  
		if PRINCIPAL_FACE_P(f,tf) 
		{
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);

				N_UDMI(v,0) = 0; 		/* Flag for sorting holdNodes */
				N_UDMI(v,1) = '\0';		/* arc length to node starting from apex */
				N_UDMI(v,2) = NODE_X(v); /* node's old x position at last time step/beginning of current time step */
				N_UDMI(v,3) = NODE_Y(v); /* node's old y position at last time step/beginning of current time step */
				N_UDMI(v,4) = '\0';		/* node's new x position */
				N_UDMI(v,5) = '\0';		/* node's new y position */
				N_UDMI(v,6) = '\0';		/* initial distance from node to bell apex, exumbrella surface */
				N_UDMI(v,7) = '\0';		/* initial distance from node to bell apex, subumbrella surface */
				N_UDMI(v,8) = '\0';		/* tparm_old */
				N_UDMI(v,9) = 1; 		/* old node flag */
				N_UDMI(v,10) = '\0';	/* tparm */
				N_UDMI(v,11) = '\0';	/* percent flex */
				N_UDMI(v,12) = NODE_X(v); /* node's UNFLEXED x position at last time step */	
				N_UDMI(v,13) = NODE_Y(v); /* node's UNFLEXED y position at last time step */
				N_UDMI(v,14) = '\0';	/* UNFLEXED arc length to node starting from apex */			
				N_UDMI(v,15) = '\0';		/* stores the % arc length for this node */
			}
		}
	}
	end_f_loop(f,tf)
	
	#endif
	return 1;
}


/*---------- Reinit_node_mem_int --------------------------------
Purpose: 

Input:	Thread *tf - face thread where node memory is being accessed
		int memslot - the N_UDMI location to be re-initialized
		int val  - node memory will be reinitialized to this integer value
					
----------------------------------------------------------------*/
int Reinit_node_mem_int(Thread *tf, int memslot, int val)
{	
	#if !RP_HOST
	face_t f;
	Node *v;
	int n;
	
	
	begin_f_loop(f,tf) 
	{  
		if PRINCIPAL_FACE_P(f,tf) 
		{
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);
				N_UDMI(v,memslot) = val; 
			}
		}
	}
	end_f_loop(f,tf)
	
	#endif
	return 1;
}

/*---------- Reinit_node_mem_double --------------------------------
Purpose: 

Input:	Thread *tf - face thread where node memory is being accessed
		int memslot - the N_UDMI location to be re-initialized
		int val  - node memory will be reinitialized to this double value
					
----------------------------------------------------------------*/
int Reinit_node_mem_double(Thread *tf, int memslot, double val)
{
	#if !RP_HOST
	face_t f;
	Node *v;
	int n;
	
	
	begin_f_loop(f,tf) 
	{  
		if PRINCIPAL_FACE_P(f,tf)
		{
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);
				N_UDMI(v,memslot) = val; 
			}
		}
	}
	end_f_loop(f,tf)
	
	#endif
	return 1;
}

/*---------- Reinit_node_mem_char --------------------------------
Purpose: 

Input:	Thread *tf - face thread where node memory is being accessed
		int memslot - the N_UDMI location to be re-initialized
		int val  - node memory will be reinitialized to this char value
					
----------------------------------------------------------------*/
int Reinit_node_mem_char(Thread *tf, int memslot, char val)
{
	#if !RP_HOST
	face_t f;
	Node *v;
	int n;
	
	
	begin_f_loop(f,tf) 
	{  
		if PRINCIPAL_FACE_P(f,tf)
		{
			f_node_loop(f,tf,n) 
			{
				v = F_NODE(f,tf,n);
				N_UDMI(v,memslot) = val; 
			}
		}
	}
	end_f_loop(f,tf)
	
	#endif
	return 1;
}


/*---------- Fix_tip_node_coords --------------------------------
Purpose: Utility to set the tip points to given coordinates before
the simulation starts.

Background: I forgot to write out the initial bell profile coordinates
with >8 digits of precision. This wouldn't be a problem if the same
point was referenced for the tip on both exum and subum during mesh
creation, but unfortunately it was not. The result is the exum and
subum tip points in the mesh differ by O(1e-8). This lack of 
precision can break the simulation depending on the values of 
kinematic variables a, b, d, g.

Input:	Thread *tf - face thread where node memory is being accessed
					
----------------------------------------------------------------*/
int Fix_tip_node_coords(Thread *tf, double x_target, double y_target)
{
	
	#if !RP_HOST
		face_t f;
		Node *this_node, *tip;
		int n;
		int k = 0;
		double this_x;
		double this_y;
		double this_euc_dist;
		double dist_max = 0.0;
		double dist_tip;
		double x_tip;
		double y_tip;
		
		/** Exumbrella **/
		begin_f_loop(f,tf) 
		{
			if PRINCIPAL_FACE_P(f,tf)
			{
				f_node_loop(f,tf,n) 
				{
					this_node = F_NODE(f,tf,n);
					this_x = NODE_X(this_node);
					this_y = NODE_Y(this_node);
					this_euc_dist = sqrt(pow(this_x,2) + pow(this_y,2));
					
					if (this_euc_dist >= dist_max)			
					{
						tip = this_node;		
						dist_max = this_euc_dist;	
						k++;
					}
				}
			} 
		}
		end_f_loop(f,tf); 
		
		k = PRF_GISUM1(k);
		if (k == 0)
		{
			Error("(Fix_tip_node_coords) Tip not found on surface! \n");
		}
		
		/* ----- x coords ------ */
		#if !RP_NODE
			dist_tip = dist_max;
		#endif
		
		/* sync if parallel */
		#if RP_NODE
			dist_tip = PRF_GRHIGH1(dist_max);
		#endif

		/* if you were the processor w/ the tip node, move the tip
			to the target coordinates */
		if ( fabs(dist_tip - dist_max) < 1e-10 )
		{
			Message("Tip node coords BEFORE move: x = %lf, y = %lf\n", NODE_X(tip), NODE_Y(tip));
			Message("Tip node coords AFTER move: x = %lf, y = %lf\n", x_target, y_target);
			Message("Tip node coords DELTA move: x = %lf, y = %lf\n", x_target - NODE_X(tip), y_target - NODE_Y(tip));
			NODE_X(tip) = x_target;
			NODE_Y(tip) = y_target;
			
		}


	#endif
	return 1;
}
