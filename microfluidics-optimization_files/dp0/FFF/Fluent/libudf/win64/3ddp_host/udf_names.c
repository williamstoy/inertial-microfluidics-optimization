/* This file generated automatically. */ 
/* Do not modify. */ 
#include "udf.h" 
#include "prop.h" 
#include "dpm.h" 
extern DEFINE_DPM_BODY_FORCE(inertial_lift, p, i);
extern DEFINE_PROFILE(inlet_x_velocity, thread, position);
extern DEFINE_INIT(precalculate_inertial_lift_coefficients_and_umax, domain);
extern DEFINE_ON_DEMAND(port_input_vars_to_scheme);
extern DEFINE_DPM_INJECTION_INIT(particle_init,I);
__declspec(dllexport) UDF_Data udf_data[] = { 
{"inertial_lift", (void (*)(void))inertial_lift, UDF_TYPE_DPM_BODY_FORCE},
{"inlet_x_velocity", (void (*)(void))inlet_x_velocity, UDF_TYPE_PROFILE},
{"precalculate_inertial_lift_coefficients_and_umax", (void (*)(void))precalculate_inertial_lift_coefficients_and_umax, UDF_TYPE_INIT},
{"port_input_vars_to_scheme", (void (*)(void))port_input_vars_to_scheme, UDF_TYPE_ON_DEMAND},
{"particle_init", (void (*)(void))particle_init, UDF_TYPE_DPM_INJECTION_INIT},
}; 
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); 
#include "version.h" 
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) 
{ 
*major = RampantReleaseMajor; 
*minor = RampantReleaseMinor; 
*revision = RampantReleaseRevision; 
} 
