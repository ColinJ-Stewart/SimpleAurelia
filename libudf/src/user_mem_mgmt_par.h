#ifndef USER_MEM_MGMT_H
#define USER_MEM_MGMT_H

#include <udf.h>
/* #include <unsteady.h> */

#include "global_var_par.h"

/* prototypes */
int Init_node_mem(Thread *tf);
int Reinit_node_mem_int(Thread *tf, int memslot, int val);
int Reinit_node_mem_double(Thread *tf, int memslot, double val);
int Reinit_node_mem_char(Thread *tf, int memslot, char val);

#endif
