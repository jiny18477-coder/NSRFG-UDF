#include "udf.h"
#include "math.h"
#include "string.h"

#define INLET_FACE_ZONE_NAME "inlet"
#define uref 10.0
#define zref 0.05
#define alfa 0.15
#define uxing 0.787
#define omu 72.6e-6
#define fai 22.25
#define I10 0.14
#define z10 10.0
#define PI 3.141592653589793
#define z0 0.05
#define SCALE  (1.0/200.0)
#define NBANDS  1000
#define DELTF   0.25
#define TAU0    1
#define C1 (10.0*2.25)   /* 22.5 */
#define C2 (12.0*2.10)   /* 25.2 */
#define C3 (12.0*2.10)   /* 25.2 */
static long s1=13579, s2=24680;
static inline real urand01_1(void){ s1 = 1664525L*s1 + 1013904223L; return (real)(s1 & 0x00FFFFFF)/16777216.0; }
static inline real urand01_2(void){ s2 = 1103515245L*s2 + 12345L;   return (real)((s2>>8) & 0x00FFFFFF)/16777216.0; }
static int   bands_ready = 0;
static real  fband[NBANDS];     
static real  th1  [NBANDS];     
static real  phi2 [NBANDS];     
static inline void ensure_bands(void){
    if (bands_ready) return;
    for (int n=0; n<NBANDS; ++n){
        fband[n] = (2.0*(n+1)-1.0) * 0.5 * DELTF; /* (2n-1)/2 Δf 
        th1  [n] = 2.0*PI * urand01_1();
        phi2 [n] = 2.0*PI * urand01_2();
    }
    bands_ready = 1;
}
static inline void nsrfg_at_point(const real X[ND_ND], real tnow,
                                  real *uout, real *vout, real *wout)
{
    const real x = (ND_ND>0? X[0]:0.0);
    const real y = (ND_ND>1? X[1]:0.0);
    real       z = (ND_ND>2? X[2]:1e-6);
    if (z < 1e-6) z = 1e-6;
    const real Uav = uref * pow(z/zref, alfa);
    const real z10s = z10 * SCALE;
    const real fcor = 2.0*omu*sin(fai*PI/180.0);
    const real h    = uxing/(6.0*fcor);
    const real hs   = h * SCALE;
    const real carg   = PI*z/(2.0*hs);
    const real c4     = pow(cos(carg), 4.0);
    const real ratio2 = 1.0 - 0.22*c4;
    const real ratio3 = 1.0 - 0.45*c4;
    const real Iu = I10 * pow(z/z10s, -alfa);
    const real Iv = Iu * ratio2;
    const real Iw = Iu * ratio3;
    const real z0s = z0;
    const real Lu0 = 300.0*SCALE * pow(z/(300.0*SCALE), 0.46 + 0.074*log(z0s));
    const real Lu  = Lu0;
    const real Lv  = 0.5 * pow(ratio2,3.0) * Lu;
    const real Lw  = 0.5 * pow(ratio3,3.0) * Lu;
    real su=0.0, sv=0.0, sw=0.0;
    for (int n=0; n<NBANDS; ++n){
        const real f = fband[n];
        const real Su = (4.0*pow(Iu*Uav,2.0)*(Lu/Uav))
                        / pow(1.0 + 70.8*pow((Lu*f)/Uav,2.0), 5.0/6.0);
        const real Sv = (4.0*pow(Iv*Uav,2.0)*(Lv/Uav))
                        * (1.0 + 188.4*pow(2.0*f*Lv/Uav,2.0))
                        / pow(1.0 + 70.8*pow(2.0*f*Lv/Uav,2.0), 11.0/6.0);
        const real Sw = (4.0*pow(Iw*Uav,2.0)*(Lw/Uav))
                        * (1.0 + 188.4*pow(2.0*f*Lw/Uav,2.0))
                        / pow(1.0 + 70.8*pow(2.0*f*Lw/Uav,2.0), 11.0/6.0);
        const real P1 = sqrt(fmax(2.0*Su*DELTF,0.0));
        const real P2 = sqrt(fmax(2.0*Sv*DELTF,0.0));
        const real P3 = sqrt(fmax(2.0*Sw*DELTF,0.0));
        const real L1 = Uav/(C1*f);
        const real L2 = Uav/(C2*f);
        const real L3 = Uav/(C3*f);
        const real Q1=P1/L1, Q2=P2/L2, Q3=P3/L3;
        const real Q1_2=Q1*Q1, Q2_2=Q2*Q2, Q3_2=Q3*Q3;
        const real A = sqrt((Q2_2+Q3_2)*(Q2_2+Q3_2) + Q1_2*Q2_2 + Q1_2*Q3_2) + 1e-16;
        const real B = sqrt(Q2_2+Q3_2) + 1e-16;
        const real th = th1[n];
        const real K1 = -(Q2_2+Q3_2)*sin(th)/A;
        const real K2 =  (Q1*Q2*sin(th))/A + Q3*cos(th)/B;
        const real K3 =  (Q1*Q3*sin(th))/A - Q2*cos(th)/B;
        const real phase = K1*(x/L1) + K2*(y/L2) + K3*(z/L3)
                         + 2.0*PI*f*(tnow/TAU0) + phi2[n];
        const real s = sin(phase);
        su += P1*s;  sv += P2*s;  sw += P3*s;
    }

    *uout = Uav + su;  
#if RP_2D
    *vout = 0.0; *wout = 0.0;
#else
    *vout = sv;  *wout = sw;
#endif
}
DEFINE_PROFILE(inlet_u_nsrfg, t, i)
{
    if (!THREAD_NAME(t) || strcmp(THREAD_NAME(t), INLET_FACE_ZONE_NAME)!=0) return;
    ensure_bands();
    face_t f;
    begin_f_loop(f, t){
        real X[ND_ND]; F_CENTROID(X,f,t);
        real u,v,w; nsrfg_at_point(X, CURRENT_TIME, &u,&v,&w);
        F_PROFILE(f,t,i) = u;
    } end_f_loop(f, t)
}
DEFINE_PROFILE(inlet_v_nsrfg, t, i)
{
    if (!THREAD_NAME(t) || strcmp(THREAD_NAME(t), INLET_FACE_ZONE_NAME)!=0) return;
    ensure_bands();
    face_t f;
    begin_f_loop(f, t){
        real X[ND_ND]; F_CENTROID(X,f,t);
        real u,v,w; nsrfg_at_point(X, CURRENT_TIME, &u,&v,&w);
        F_PROFILE(f,t,i) = v;
    } end_f_loop(f, t)
}
DEFINE_PROFILE(inlet_w_nsrfg, t, i)
{
    if (!THREAD_NAME(t) || strcmp(THREAD_NAME(t), INLET_FACE_ZONE_NAME)!=0) return;
    ensure_bands();
    face_t f;
    begin_f_loop(f, t){
        real X[ND_ND]; F_CENTROID(X,f,t);
        real u,v,w; nsrfg_at_point(X, CURRENT_TIME, &u,&v,&w);
        F_PROFILE(f,t,i) = w;
    } end_f_loop(f, t)
}
