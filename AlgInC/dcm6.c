#include "dcm6.h"
#include "math_port.h"
#include "assert.h"
///////////////////////////////////////////////////////
#define  deltat		0.001f
#define  deltat2	0.000001f
#define	 PIx2		(2.0f*PI)

static float r3x3[9];
static float p6x6[36];
static float rot3x3[9] = {
	1.0F,0,0,
	0,1.0F,0,
	0,0,1.0F
};

struct dcm6{
	//complimentary filter param
	float Ki;//0.001
	float Kp;//0.01
	float err_b;
	//dcm param
	float g0;
	float q_dcm2;//0.00003
	float q_gyro_bias2;//0.001^2
	float r_acc2;//0.5^2
	float r_a2;//10^2
	float q_dcm2_init;//1^2
	float q_gyro_bias2_init;//0.1^2
	float state[6];//(C31 C32 C33 wb1 wb2 wb3)
	float accel_global[3];//x,y,z, in g
	float yaw;
	float pitch;
	float roll;
	float first_row[3];//first row in DCM matrix
	
	float Q[6];//simplified here. 6x6 for matrix
	float H;//simplified here. 3x6 for matrix
	float R;//simplified here. 3x3 for matrix
	matrix_instance_t P;//6x6
	matrix_instance_t Rot;//3x3
};

static struct dcm6 obj;
/////////////////////////////////////////////////////////////////////
void dcm6_init(void)
{
	//complimentary filter
	obj.Ki = 0.001f;
	obj.Kp = 0.01f;
	//dcm
	obj.g0 = 9.8f;
	obj.q_dcm2 = 0.00003f;
	obj.q_gyro_bias2 = 1e-6f;
	obj.r_acc2 = 0.25f;
	obj.r_a2 = 100.0f;
	obj.q_dcm2_init = 1.0f;
	obj.q_gyro_bias2_init = 0.01f;
	//arm_mat_init_f32(&obj.R,3,3,r3x3);
	arm_mat_init_f32(&obj.P,6,6,p6x6);
	arm_mat_init_f32(&obj.Rot,3,3,rot3x3);
	obj.P.pData[0] = obj.q_dcm2_init;
	obj.P.pData[7] = obj.q_dcm2_init;
	obj.P.pData[14] = obj.q_dcm2_init;
	obj.P.pData[21] = obj.q_gyro_bias2_init;
	obj.P.pData[28] = obj.q_gyro_bias2_init;
	obj.P.pData[35] = obj.q_gyro_bias2_init;
	obj.H = obj.g0;
	obj.first_row[0] = 1.0f;
	obj.first_row[1] = 0;
	obj.first_row[2] = 0;
	obj.state[2] = 1.0f;//all others are 0
	obj.Q[0]=obj.Q[1]=obj.Q[2]=obj.q_dcm2*deltat2;
	obj.Q[3]=obj.Q[4]=obj.Q[5]=obj.q_gyro_bias2*deltat2;
	return;
}

int32_t DCM6_updateIMU9(float ux, float uy, float uz, 
					float ax, float ay, float az, 
					float mx, float my, float mz, uint32_t mag_update,
					float dt)
{
	uint8_t i,j;
	float tmp3x3[9],tmp6x6[36],tmp3x3_[9],ktmp[36];
//	float tmp_s;//scale tmp
	float P_predict[36];
	float *x,x_predict[6];
	x = obj.state;//point to obj.state

	float ux_bx = (ux - x[3])*dt;
	float uy_by = (uy - x[4])*dt;
	float uz_bz = (uz - x[5])*dt;
	float x0dt = x[0]*dt;
	float x1dt = x[1]*dt;
	float x2dt = x[2]*dt;
	//predict
	x_predict[0] = x[0] - x[2]*uy_by + x[1]*uz_bz;
	x_predict[1] = x[1] + x[2]*ux_bx - x[0]*uz_bz;
	x_predict[2] = x[2] - x[1]*ux_bx + x[0]*uy_by;
	x_predict[3] = x[3];
	x_predict[4] = x[4];
	x_predict[5] = x[5];
	//F, jacobian matrix
	float f[36] = {
		1.0f,	uz_bz,	-uy_by,	0,	x2dt,	-x1dt,
		-uz_bz,	1.0f,	ux_bx,	-x2dt,	0,	x0dt,
		uy_by,	-ux_bx,	1.0f,	x1dt,	-x0dt,	0,
		0,0,0,			1.0f,0,0,
		0,0,0,			0,1.0f,0,
		0,0,0,			0,0,1.0f
	};

	matrix_instance_t F,tmp_mx,P_pre,S,Sinv,K,IKH;
	
	arm_mat_init_f32(&F,6,6,f);
	arm_mat_init_f32(&tmp_mx,6,6,tmp6x6);
	arm_mat_init_f32(&P_pre,6,6,P_predict);
	//F*P*F'
	arm_mat_mult_f32(&F,&obj.P,&tmp_mx);
	f[1] = -f[1];
	f[2] = -f[2];
	f[4] = f[5] = 0;
	f[6] = -f[6];
	f[8] = -f[8];
	f[9] = f[11] = 0;
	f[12] = -f[12];
	f[13] = -f[13];
	f[15] = f[16] = 0;
	f[19] = -x2dt;
	f[20] = x1dt;
	f[24] = x2dt;
	f[26] = -x0dt;
	f[30] = -x1dt;
	f[31] = x0dt;
//	ft = {
//		1.0f,	-uz_bz,	uy_by,	0,	0,	0,
//		uz_bz,	1.0f,	-ux_bx,	0,	0,	0,
//		-uy_by,	ux_bx,	1.0f,	0,	0,	0,
//		0,	-x2dt,	x1dt,			1.0f,0,0,
//		x2dt,	0,	-x0dt,			0,1.0f,0,
//		-x1dt,	x0dt,	0,			0,0,1.0f
//	};
	arm_mat_mult_f32(&tmp_mx,&F,&P_pre);
	//P_prediction = F*P*F' + Q
	matrix_diag_plus(&P_pre, obj.Q);
	float dev[3];
	dev[0] = (ax-x_predict[0])*obj.g0;
	dev[1] = (ay-x_predict[1])*obj.g0;
	dev[2] = (az-x_predict[2])*obj.g0;
	//size of linear accel(m/s2)
	float val_size;
	arm_sqrt_f32(dev[0]*dev[0]+dev[1]*dev[1]+dev[2]*dev[2], &val_size);
	obj.R = val_size*obj.r_a2 + obj.r_acc2;
	//Kalman innovation y(m/s2), same to dev above
	//y = a - H*x_predict;
	arm_mat_init_f32(&tmp_mx,6,3,tmp6x6);
	//tmp_mx = P_pre*obj.H';simplified below for higher performance.
	for(i=0;i<6;++i){
		//tmp6x6 here use as 6x3,P is 6x6
		tmp6x6[i*3] = P_predict[i*6]*obj.H;
		tmp6x6[i*3+1] = P_predict[i*6+1]*obj.H;
		tmp6x6[i*3+2] = P_predict[i*6+2]*obj.H;
	}
	arm_mat_init_f32(&S,3,3,tmp3x3);
	arm_mat_init_f32(&Sinv,3,3,tmp3x3_);
	//S = obj.H * P_pre * obj.H' + R -> S = obj.H * tmp_mx + R
	for(i=0;i<3;++i){
		tmp3x3[i*3] = obj.H*tmp6x6[i*3];
		tmp3x3[i*3+1] = obj.H*tmp6x6[i*3+1];
		tmp3x3[i*3+2] = obj.H*tmp6x6[i*3+2];
	}
	tmp3x3[0] += obj.R;
	tmp3x3[4] += obj.R;
	tmp3x3[8] += obj.R;
	//Kalman gain = tmp/S -> tmp*Sinv
	if(arm_mat_inverse_f32(&S,&Sinv))//!=ARM_MATH_SUCCESS IS ERR
		assert(0);//mismatch or non invertible
		
	arm_mat_init_f32(&K,6,3,ktmp);
	arm_mat_mult_f32(&tmp_mx,&Sinv,&K);
	//update a posteriori x = x_pre + K*y; y same to dev here 
	for(i=0;i<6;++i){
		//ktmp array here is content for K
		x[i] = x_predict[i] + ktmp[i*3]*dev[0]+ktmp[i*3+1]*dev[1]+ktmp[i*3+2]*dev[2];
	}
	//tmp6x6,tmp3x3_,tmp3x3,f available
	arm_mat_init_f32(&IKH,6,6,f);
	//IKH = eye(6) - K*obj.H; content in f is very similar to IKH,
	for(i = 0;i<6;++i){
		f[i*6] = -ktmp[i*3] * obj.H;
		f[i*6+1] = -ktmp[i*3+1] * obj.H;
		f[i*6+2] = -ktmp[i*3+2] * obj.H;
	}
	f[0] += 1.0f; 
	f[7] += 1.0f; 
	f[14] += 1.0f; 
	//obj.P = IKH * P_pre * IKH' + K*R*K';
	arm_mat_init_f32(&tmp_mx,6,6,tmp6x6);
	arm_mat_mult_f32(&IKH,&P_pre,&tmp_mx);
	for(i = 0;i<6;++i){
		for(j = 0;j<6;++j){
			P_predict[i*6+j] = f[j*6+i];//IKH'(i,j) = IKH(j,i)
		}
	}
	arm_mat_mult_f32(&tmp_mx,&P_pre,&IKH);//IKH content IKH*P*IKH'
	//K*R*K'
	for(i = 0;i<3;++i){
		for(j = 0;j<6;++j){
			tmp6x6[i*6+j] = obj.R*ktmp[j*3+i];//K'(i,j) = K(j,i)
		}
	}
	arm_mat_init_f32(&tmp_mx,3,6,tmp6x6);//set tmp_mx as K'
	arm_mat_mult_f32(&K,&tmp_mx,&P_pre);//P_pre content K*R*K'
	//update a posteriori covariance
	arm_mat_add_f32(&IKH,&P_pre,&obj.P);
	//normalize x & P (divide by DCM vector length)
	float x0x0 = x[0]*x[0];
	float x1x1 = x[1]*x[1];
	float x2x2 = x[2]*x[2];	
	float x0x1 = x[0]*x[1];
	float x0x2 = x[0]*x[2];
	float x1x2 = x[1]*x[2];
	val_size = rsqrt(x0x0+x1x1+x2x2);//dcm vector length.
	float val_size3 = val_size*val_size*val_size;
	//J
	tmp6x6[0] = (x1x1 + x2x2)*val_size3;
	tmp6x6[1] = -x0x1*val_size3;
	tmp6x6[2] = -x0x2*val_size3;
	tmp6x6[3] = tmp6x6[4] = tmp6x6[5] = 0;
	tmp6x6[6] = -x0x1*val_size3;
	tmp6x6[7] = (x0x0 + x2x2)*val_size3;
	tmp6x6[8] = -x1x2*val_size3;
	tmp6x6[9] = tmp6x6[10] = tmp6x6[11] = 0;
	tmp6x6[12] = -x0x2*val_size3;
	tmp6x6[13] = -x1x2*val_size3;
	tmp6x6[14] = (x0x0 + x1x1)*val_size3;
	tmp6x6[15] = tmp6x6[16] = tmp6x6[17] = 0;
	tmp6x6[18] = tmp6x6[19] = tmp6x6[20] = 0;
	tmp6x6[21] = 1.0f;
	tmp6x6[22] = tmp6x6[23] = 0;
	tmp6x6[24] = tmp6x6[25] = tmp6x6[26] = tmp6x6[27] = tmp6x6[29] = 0;
	tmp6x6[28] = 1.0f;
	tmp6x6[30] = tmp6x6[31] = tmp6x6[32] = tmp6x6[33] = tmp6x6[34] = 0;
	tmp6x6[35] = 1.0f;
	//obj.P = J*P*J';J = J'
	arm_mat_mult_f32(&tmp_mx,&obj.P,&P_pre);//P_pre = J*P
	arm_mat_mult_f32(&P_pre,&tmp_mx,&obj.P);//obj.P = P_pre*J'
	x[0] *= val_size;
	x[1] *= val_size;
	x[2] *= val_size;
	//update UX
	ux_bx = (ux - x[3])*dt;
	uy_by = (uy - x[4])*dt;
	uz_bz = (uz - x[5])*dt;
	//x1 with gyro data
	tmp3x3[0] = obj.first_row[0] + uz_bz*obj.first_row[1] - uy_by*obj.first_row[2];//x1(1)
	tmp3x3[1] = obj.first_row[1] - uz_bz*obj.first_row[0] + ux_bx*obj.first_row[2];//x1(2)
	tmp3x3[2] = obj.first_row[2] + uy_by*obj.first_row[0] - ux_bx*obj.first_row[1];//x1(3)
	//normalize x1
	float rnorm = rsqrt(tmp3x3[0]*tmp3x3[0]+tmp3x3[1]*tmp3x3[1]+tmp3x3[2]*tmp3x3[2]);
	tmp3x3[0] *= rnorm;
	tmp3x3[1] *= rnorm;
	tmp3x3[2] *= rnorm;
	//x2
	tmp3x3[3] = -x[2]*tmp3x3[1] + x[1]*tmp3x3[2];//x2(1)
	//tmp3x3[4] = x[2]*tmp3x3[0] - x[0]*tmp3x3[2];//x2(2)
	//tmp3x3[5] = -x[1]*tmp3x3[0] + x[0]*tmp3x3[1];//x2(3)
	//euler
	tmp3x3_[0] = atan2f(tmp3x3[3],tmp3x3[0]);//yaw_g // fatan2_deg
	//obj.pitch = asinf(-x[0]);//fasin_deg, not need in this turn
	obj.roll = atan2f(x[1],x[2]);//fatan2_deg
	tmp3x3_[1] = arm_cos_f32(obj.roll);//cos_roll 
	if(x[2] == 0.0f){
		if(x[1]>0)
			tmp3x3_[2] = 1.0f;
		else 
			tmp3x3_[2] = -1.0f;
		return -1;//error
	}else
		tmp3x3_[2] = tmp3x3_[1]*x[1]/x[2];//sin_roll
	
	tmp3x3_[3] = x[1]/tmp3x3_[2];//cos_pitch
	if(mag_update){
		//complimentary filter
		tmp3x3[0] = x[2]*tmp3x3_[1]+x[1]*tmp3x3_[2];//x1(1)
		tmp3x3[1] = -x[0]*tmp3x3_[2];//x1(2)
		tmp3x3[2] = -x[0]*tmp3x3_[1];//x1(3)
		tmp3x3[3] = tmp3x3[0]*mx+tmp3x3[1]*my+tmp3x3[2]*mz;//magx
		tmp3x3[4] = tmp3x3_[1]*my-tmp3x3_[2]*mz;//magy
		tmp3x3_[4] = atan2f(-tmp3x3[4],tmp3x3[3]);//yaw_m
		tmp3x3_[5] = tmp3x3_[4] - obj.yaw;//yaw_m-obj.yaw
		//solve the discontinue between angle (-pi) and pi.
		if(tmp3x3_[5] > PI)
			tmp3x3_[5] -= PIx2;
		else if(tmp3x3_[5] < -PI)
			tmp3x3_[5] += PIx2;
		
		obj.err_b = obj.err_b + obj.Ki*dt*(tmp3x3_[5]);
		obj.yaw = tmp3x3_[0] + obj.Kp*tmp3x3_[5] + obj.err_b;
	}else{
		obj.yaw = tmp3x3_[0];
	}
	//limit yaw between -pi ~ pi
	if(obj.yaw>PI)
		obj.yaw = PI;
	else if(obj.yaw<-PI)
		obj.yaw = -PI;
	
	float cos_yaw = arm_cos_f32(obj.yaw);
	float sin_yaw = arm_sin_f32(obj.yaw);
	//update first row
	obj.first_row[0] = tmp3x3_[3]*cos_yaw;
	obj.first_row[1] = -tmp3x3_[1]*sin_yaw-tmp3x3_[2]*x[0]*cos_yaw;
	obj.first_row[2] = tmp3x3_[2]*sin_yaw-tmp3x3_[1]*x[0]*cos_yaw;
	//update x2
	tmp3x3[3] = -x[2]*obj.first_row[1] + x[1]*obj.first_row[2];//x2(1)
	tmp3x3[4] = x[2]*obj.first_row[0] - x[0]*obj.first_row[2];//x2(2)
	tmp3x3[5] = -x[1]*obj.first_row[0] + x[0]*obj.first_row[1];//x2(3)
	obj.accel_global[0] = (obj.first_row[0]*ax+obj.first_row[1]*ay+obj.first_row[2]*az);
	obj.accel_global[1] = (tmp3x3[3]*ax+tmp3x3[4]*ay+tmp3x3[5]*az);
	obj.accel_global[2] = (x[0]*ax+x[1]*ay+x[2]*az - 1.0f);
	
	return 0;
}

// return dat in (m/s^2)
void DCM6_getGlobalAccel(float *dat)
{
	dat[0] = obj.accel_global[0]*obj.g0;
	dat[1] = obj.accel_global[1]*obj.g0;
	dat[2] = obj.accel_global[2]*obj.g0;
}
void DCM6_getAttitude(float *dat)
{
	obj.pitch = asinf(-obj.state[0]);//fasin_deg
	dat[0] = obj.pitch*rad2deg;//(-90 ~ 90)
	dat[1] = obj.roll*rad2deg;//(-180 ~ 180)
	dat[2] = obj.yaw*rad2deg;//(-180 ~ 180)
}

