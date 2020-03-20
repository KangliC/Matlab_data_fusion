/******************************************************
*include files for dcm6 with complimentary filter
*
*data type and defines
******************************************************/
#ifndef DCM6_H_
#define DCM6_H_
///////////////////////////////////////////////////////
#include "math_port_config.h"
///////////////////////////////////////////////////////
#define CHX 0
#define CHY 1
#define CHZ 2

extern void dcm6_init(void);

extern int32_t DCM6_updateIMU9(float ux, float uy, float uz, 
					float ax, float ay, float az, 
					float mx, float my, float mz, uint32_t mag_update,
					float dt);
					
extern void DCM6_getGlobalAccel(float *dat);
extern void DCM6_getAttitude(float *dat);				
#endif
