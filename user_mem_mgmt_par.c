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
				N_UDMI(v,15) = 0.0;		
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
