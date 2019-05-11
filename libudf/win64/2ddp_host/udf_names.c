/* This file generated automatically. */ 
/* Do not modify. */ 
#include "udf.h" 
#include "prop.h" 
#include "dpm.h" 
extern DEFINE_EXECUTE_ON_LOADING(get_domain, libname);
extern DEFINE_INIT(init_node_mem, domain);
extern DEFINE_GRID_MOTION(Jelly_motion,domain,dt,time,dtime);
extern DEFINE_GEOM(axis_top, domain, dt, position);
extern DEFINE_GEOM(axis_bot, domain, dt, position);
extern DEFINE_CG_MOTION(CG_motion_intfc, dt, vel, omega, time, dtime);
extern DEFINE_EXECUTE_AT_END(forces_at_end);
__declspec(dllexport) UDF_Data udf_data[] = { 
{"get_domain", (void (*)())get_domain, UDF_TYPE_EXECUTE_ON_LOADING},
{"init_node_mem", (void (*)())init_node_mem, UDF_TYPE_INIT},
{"Jelly_motion", (void (*)())Jelly_motion, UDF_TYPE_GRID_MOTION},
{"axis_top", (void (*)())axis_top, UDF_TYPE_GEOM},
{"axis_bot", (void (*)())axis_bot, UDF_TYPE_GEOM},
{"CG_motion_intfc", (void (*)())CG_motion_intfc, UDF_TYPE_CG_MOTION},
{"forces_at_end", (void (*)())forces_at_end, UDF_TYPE_EXECUTE_AT_END},
}; 
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); 
#include "version.h" 
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) 
{ 
*major = RampantReleaseMajor; 
*minor = RampantReleaseMinor; 
*revision = RampantReleaseRevision; 
} 
