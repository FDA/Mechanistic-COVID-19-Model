#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <time.h>
static double parms[32];
#define K14 parms[0]
#define K16 parms[1]
#define K13 parms[2]
#define K11 parms[3]
#define K18 parms[4]
#define K48 parms[5]
#define K29 parms[6]
#define K15 parms[7]
#define K42 parms[8]
#define K34 parms[9]
#define K40 parms[10]
#define K31 parms[11]
#define K27 parms[12]
#define K24 parms[13]
#define K47 parms[14]
#define K26 parms[15]
#define K25 parms[16]
#define K28 parms[17]
#define K288 parms[18]
#define K45 parms[19]
#define K37 parms[20]
#define K32 parms[21]
#define K322 parms[22]
#define K3222 parms[23]
#define unnamed parms[24]
#define K49 parms[25]
#define K53 parms[26]
#define K57 parms[27]
#define K59 parms[28]
#define K245 parms[29]
#define timeout parms[30]
#define starttime parms[31]

void lagvalue(double *T, int *nr, int N, double *yout) {
static void(*fun)(double*, int*, int, double*) = NULL;
if(fun==NULL)
fun =  (void(*)(double*, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
return fun(T, nr, N, yout);
}
void lagderiv(double *T, int *nr, int N, double *yout) {
static void(*fun)(double*, int*, int, double*) = NULL;
if (fun == NULL)
fun =  (void(*)(double*, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
return fun(T, nr, N, yout);
}

void initmod(void (* odeparms)(int *, double *)){
int N=32;
odeparms(&N, parms);}

void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
if (ip[0] < 3 ) error("nout not enough!");
time_t s = time(NULL);
if((int) s - (int) starttime > timeout) error("timeout!");
double ytau[3] = {0,0,0};
int Nout1 = 1;
int nr1[1]= {8};
double ytau1[1] = {0};
double T1 = *t - 12;
if (*t >  12 ) {
lagvalue(&T1, nr1, Nout1, ytau1);
}
ytau[0] = ytau1[0];
int Nout2 = 1;
int nr2[1]= {4};
double ytau2[1] = {0};
double T2 = *t - 24;
if (*t >  24 ) {
lagvalue(&T2, nr2, Nout2, ytau2);
}
ytau[1] = ytau2[0];
int Nout3 = 1;
int nr3[1]= {1};
double ytau3[1] = {0};
double T3 = *t - 48;
if (*t >  48 ) {
lagvalue(&T3, nr3, Nout3, ytau3);
}
ytau[2] = ytau3[0];
ydot[0] = 1/unnamed*(-K14*y[0]+K13*y[2]-K34*y[0]*0.15*ytau[0]-K28*y[7]*y[0]-K288*y[5]*y[0]);
ydot[1] = 1/unnamed*(K16*y[0]-K18*y[1]);
ydot[2] = 1/unnamed*(-K13*y[2]+K11*y[1]*y[3]-K34*y[2]*0.15*ytau[0]-K28*y[7]*y[2]-K288*y[5]*y[2]);
ydot[3] = 1/unnamed*(-K11*y[1]*y[3]+K27*y[6]-K26*y[7]*y[3]-K288*y[5]*y[3]);
ydot[4] = 1/unnamed*(K48*y[1]*y[10]-K53*y[4]-0);
ydot[5] = 1/unnamed*(-K29*y[5]+K15*y[4]+K49*y[0]+K245*ytau[2]);
ydot[6] = 1/unnamed*(-K27*y[6]+K26*y[7]*y[3]-K288*y[5]*y[6]);
ydot[7] = 1/unnamed*(K24*0.15*ytau[0]+K47*y[0]+K57*y[4]-K25*y[7]);
ydot[8] = 1/unnamed*(-K45*y[8]+K32*y[9]*y[14]-K3222*y[1]*y[8]);
ydot[9] = 1/unnamed*(K37*ytau[1]-K42*y[9]);
ydot[10] = K59*(1E5-y[10])-K48*y[1]*y[10];
ydot[11] = K28*y[7]*y[0]+K28*y[7]*y[2];
ydot[12] = K288*y[5]*y[0]+K288*y[5]*y[2]+K288*y[5]*y[3]+K288*y[5]*y[6];
ydot[13] = K34*y[0]*0.15*ytau[0]+K34*y[2]*0.15*ytau[0];
ydot[14] = K322*(1E5-y[14])-K3222*y[1]*y[14];
yout[0] = ytau[0];
yout[1] = ytau[1];
yout[2] = ytau[2];
}
