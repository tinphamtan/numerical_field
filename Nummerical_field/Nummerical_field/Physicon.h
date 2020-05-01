#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define C_0       2.99792458E+8
#define Cvak      C_0
#define Z_0       376.991
#define PI        3.14159265358979
#define K_PHI     0.0174533
#define MY_0      4*PI*1.E-7
#define EPS_0     8.8542E-12

#define EVecX  RVector(1.,0.,0.)
#define EVecY  RVector(0.,1.,0.)
#define EVecZ  RVector(0.,0.,1.)

#define QUAT(x)  ((x)*(x))
#define CUBE(x)  ((x)*(x)*(x))

//Leitfähigkeit verschiedener Metalle
// in Siemens/Meter
#define SIGMA_CU  5.88E+7
#define SIGMA_AG  6.21E+7
#define SIGMA_AU  4.55E+7
