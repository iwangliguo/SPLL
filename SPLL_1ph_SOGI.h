#ifndef SPLL_1ph_SOGI_H_
#define SPLL_1ph_SOGI_H_
#define SPLL_SOGI_Q _IQ23
#define SPLL_SOGI_Qmpy _IQ23mpy
#define SPLL_SOGI_SINE _IQ23sin
#define SPLL_SOGI_COS _IQ23cos
//*********** Structure Definition ********//
typedef struct{
int32 osg_k;
int32 osg_x;
int32 osg_y;
int32 osg_b0;
int32 osg_b2;
int32 osg_a1;
int32 osg_a2;
int32 osg_qb0;
int32 osg_qb1;
int32 osg_qb2;
}SPLL_SOGI_OSG_COEFF;
typedef struct{
int32 B1_lf;
int32 B0_lf;
int32 A1_lf;
}SPLL_SOGI_LPF_COEFF;
typedef struct{
int32 u[3]; // Ac Input
int32 osg_u[3];
int32 osg_qu[3];
int32 u_Q[2];
int32 u_D[2];
int32 ylf[2];
int32 fo; // output frequency of PLL
int32 fn; //nominal frequency
int32 theta[2];
int32 cos;
int32 sin;
int32 delta_T;
SPLL_SOGI_OSG_COEFF osg_coeff;
SPLL_SOGI_LPF_COEFF lpf_coeff;
}SPLL_1ph_SOGI;
//*********** Function Declarations *******//
inline void SPLL_1ph_SOGI_init(int Grid_freq, long DELTA_T, volatile SPLL_1ph_SOGI *spll,
volatile SPLL_SOGI_LPF_COEFF lpf_coeff);
inline void SPLL_1ph_SOGI_coeff_update(float delta_T, float wn, volatile SPLL_1ph_SOGI *spll);
inline void SPLL_1ph_SOGI_run_FUNC(SPLL_1ph_SOGI *spll1);
//*********** Structure Init Function ****//
inline void SPLL_1ph_SOGI_init(int Grid_freq, long DELTA_T, volatile SPLL_1ph_SOGI *spll_obj,
volatile SPLL_SOGI_LPF_COEFF lpf_coeff)
{
spll_obj->u[0]=SPLL_SOGI_Q(0.0);
spll_obj->u[1]=SPLL_SOGI_Q(0.0);
spll_obj->u[2]=SPLL_SOGI_Q(0.0);
spll_obj->osg_u[0]=SPLL_SOGI_Q(0.0);
spll_obj->osg_u[1]=SPLL_SOGI_Q(0.0);
spll_obj->osg_u[2]=SPLL_SOGI_Q(0.0);
spll_obj->osg_qu[0]=SPLL_SOGI_Q(0.0);
spll_obj->osg_qu[1]=SPLL_SOGI_Q(0.0);
spll_obj->osg_qu[2]=SPLL_SOGI_Q(0.0);
spll_obj->u_Q[0]=SPLL_SOGI_Q(0.0);
spll_obj->u_Q[1]=SPLL_SOGI_Q(0.0);
spll_obj->u_D[0]=SPLL_SOGI_Q(0.0);
spll_obj->u_D[1]=SPLL_SOGI_Q(0.0);
spll_obj->ylf[0]=SPLL_SOGI_Q(0.0);
spll_obj->ylf[1]=SPLL_SOGI_Q(0.0);
spll_obj->fo=SPLL_SOGI_Q(0.0);
spll_obj->fn=SPLL_SOGI_Q(Grid_freq);
spll_obj->theta[0]=SPLL_SOGI_Q(0.0);
spll_obj->theta[1]=SPLL_SOGI_Q(0.0);
spll_obj->sin=SPLL_SOGI_Q(0.0);
spll_obj->cos=SPLL_SOGI_Q(0.0);
//coefficients for the loop filter
spll_obj->lpf_coeff.B1_lf=lpf_coeff.B1_lf;
spll_obj->lpf_coeff.B0_lf=lpf_coeff.B0_lf;
spll_obj->lpf_coeff.A1_lf=lpf_coeff.A1_lf;
spll_obj->delta_T=DELTA_T;
}
//*********** Structure Coeff Update *****//
inline void SPLL_1ph_SOGI_coeff_update(float delta_T, float wn, volatile SPLL_1ph_SOGI *spll)
{
float osgx,osgy,temp;
spll->osg_coeff.osg_k=SPLL_SOGI_Q(0.5);
osgx=(float)(2.0*0.5*wn*delta_T);
spll->osg_coeff.osg_x=SPLL_SOGI_Q(osgx);
osgy=(float)(wn*delta_T*wn*delta_T);
spll->osg_coeff.osg_y=SPLL_SOGI_Q(osgy);
temp=(float)1.0/(osgx+osgy+4.0);
spll->osg_coeff.osg_b0=SPLL_SOGI_Q((float)osgx*temp);
spll->osg_coeff.osg_b2=SPLL_SOGI_Qmpy(SPLL_SOGI_Q(-1.0),spll->osg_coeff.osg_b0);
spll->osg_coeff.osg_a1=SPLL_SOGI_Q((float)(2.0*(4.0-osgy))*temp);
spll->osg_coeff.osg_a2=SPLL_SOGI_Q((float)(osgx-osgy-4)*temp);
spll->osg_coeff.osg_qb0=SPLL_SOGI_Q((float)(0.5*osgy)*temp);
spll->osg_coeff.osg_qb1=SPLL_SOGI_Qmpy(spll-
>osg_coeff.osg_qb0,SPLL_SOGI_Q(2.0));
spll->osg_coeff.osg_qb2=spll->osg_coeff.osg_qb0;
}
//*********** Function Definition ********//
inline void SPLL_1ph_SOGI_run_FUNC(SPLL_1ph_SOGI *spll_obj)
{
// Update the spll_obj->u[0] with the grid value before calling this routine
//-------------------------------//
// Orthogonal Signal Generator //
//-------------------------------//
spll_obj->osg_u[0]=SPLL_SOGI_Qmpy(spll_obj->osg_coeff.osg_b0,(spll_obj->u[0]-
spll_obj->u[2]))+SPLL_SOGI_Qmpy(spll_obj->osg_coeff.osg_a1,spll_obj-
>osg_u[1])+SPLL_SOGI_Qmpy(spll_obj->osg_coeff.osg_a2,spll_obj->osg_u[2]);
spll_obj->osg_u[2]=spll_obj->osg_u[1];
spll_obj->osg_u[1]=spll_obj->osg_u[0];
spll_obj->osg_qu[0]=SPLL_SOGI_Qmpy(spll_obj->osg_coeff.osg_qb0,spll_obj-
>u[0])+SPLL_SOGI_Qmpy(spll_obj->osg_coeff.osg_qb1,spll_obj->u[1])+SPLL_SOGI_Qmpy(spll_obj-
>osg_coeff.osg_qb2,spll_obj->u[2])+SPLL_SOGI_Qmpy(spll_obj-
>osg_coeff.osg_a1,spll_obj->osg_qu[1])+SPLL_SOGI_Qmpy(spll_obj-
>osg_coeff.osg_a2,spll_obj->osg_qu[2]);
spll_obj->osg_qu[2]=spll_obj->osg_qu[1];
spll_obj->osg_qu[1]=spll_obj->osg_qu[0];
spll_obj->u[2]=spll_obj->u[1];
spll_obj->u[1]=spll_obj->u[0];
//-------------------------------------------------------//
// Park Transform from alpha beta to d-q axis //
//-------------------------------------------------------//
spll_obj->u_Q[0]=SPLL_SOGI_Qmpy(spll_obj->cos,spll_obj-
>osg_u[0])+SPLL_SOGI_Qmpy(spll_obj->sin,spll_obj->osg_qu[0]);
spll_obj->u_D[0]=SPLL_SOGI_Qmpy(spll_obj->cos,spll_obj->osg_qu[0])-
SPLL_SOGI_Qmpy(spll_obj->sin,spll_obj->osg_u[0]);
//---------------------------------//
// Loop Filter //
/---------------------------------//
spll_obj->ylf[0]=spll_obj->ylf[1]+SPLL_SOGI_Qmpy(spll_obj-
>lpf_coeff.B0_lf,spll_obj->u_Q[0])+SPLL_SOGI_Qmpy(spll_obj->lpf_coeff.B1_lf,spll_obj->u_Q[1]);
//spll_obj->ylf[0]=(spll_obj->ylf[0]>SPLL_SOGI_Q(20.0))?SPLL_SOGI_Q(20.0):spll_obj-
>ylf[0];
//spll_obj->ylf[0]=(spll_obj->ylf[0] < SPLL_SOGI_Q(-20.0))?SPLL_SOGI_Q(-
20.0):spll_obj->ylf[0];
spll_obj->ylf[1]=spll_obj->ylf[0];
spll_obj->u_Q[1]=spll_obj->u_Q[0];
//spll_obj->u_D[1]=spll_obj->u_D[0];
//---------------------------------//
// VCO //
//---------------------------------//
spll_obj->fo=spll_obj->fn+spll_obj->ylf[0];
spll_obj->theta[0]=spll_obj->theta[1]+SPLL_SOGI_Qmpy(SPLL_SOGI_Qmpy(spll_obj-
>fo,spll_obj->delta_T),SPLL_SOGI_Q(2*3.1415926));
if(spll_obj->theta[0]>SPLL_SOGI_Q(2*3.1415926))
spll_obj->theta[0]=spll_obj->theta[0]-SPLL_SOGI_Q(2*3.1415926);
spll_obj->theta[1]=spll_obj->theta[0];
spll_obj->sin=SPLL_SOGI_SINE(spll_obj->theta[0]);
spll_obj->cos=SPLL_SOGI_COS(spll_obj->theta[0]);
}
//*********** Macro Definition ***********//
#define SPLL_1ph_SOGI_run_MACRO(v) \
v.osg_u[0]=SPLL_SOGI_Qmpy(v.osg_coeff.osg_b0,(v.u[0]-
v.u[2]))+SPLL_SOGI_Qmpy(v.osg_coeff.osg_a1,v.osg_u[1])+SPLL_SOGI_Qmpy(v.osg_coeff.osg_a2,v.osg_u[2
]);\
v.osg_u[2]=v.osg_u[1];\
v.osg_u[1]=v.osg_u[10];\
v.osg_qu[0]=SPLL_SOGI_Qmpy(v.osg_coeff.osg_qb0,v.u[0])+SPLL_SOGI_Qmpyu(v.osg_coeff.osg_qb1,v.u[1])
+SPLL_SOGI_Qmpy(v.osg_coeff.osg_qb2,v.u[2])+SPLL_SOGI_Qmpy(v.osg_coeff.osg_a1,v.osg_qu[1])+
SPLL_SOGI_Qmpy(v.osg_coeff.osg_a2.v.osg_qu[2]);
\
v.osg_qu[2]=v.osg_qu[1]; \
v.osg_qu[1]=v.osg_qu[0]; \
v.u[2]=v.u[1]; \
v.u[1]=v.u[0]; \
v.u_Q[0]=SPLL_SOGI_Qmpy(v.cos,v.osg_u[0])+SPLL_SOGI_Qmpy(v.sin,v.osg_qu[0]); \
v.u_D[0]=SPLL_SOGI_Qmpy(v.cos,v.osg_qu[0])-SPLL_SOGI_Qmpy(v.sin,v.osg_u[0]); \
v.ylf[0]=v.ylf[1]+SPLL_SOGI_Qmpy(v.lpf_coeff.B0_lf,v.u_Q[0])+SPLL_SOGI_Qmpy(v.lpf_coeff.B1_lf,v.u_
Q[1]); \
v.ylf[1]=v.ylf[0]; \
v.u_Q[1]=v.u_Q[0]; \
v.fo=v.fn+v.ylf[0]; \
v.theta[0]=v.theta[1]+SPLL_SOGI_Qmpy(SPLL_SOGI_Qmpy(v.fov.v.delta_T),SPLL_Q(2*3.1415926)); \
if(v.theta[0]>SPLL_SOGI_Q(2*3.1415926)) \
v.theta[0]=SPLL_SOGI_Q(0.0); \
v.theta[1]=v.theta[0]; \
v.sin=SPLL_SOGI_SINE(v.theta[0]); \
v.cos=SPLL_SOGI_COS(v.theta[0]);
#endif/* SPLL_1ph_SOGI_H_ */