#ifndef MYMATLIB_H
#define MYMATLIB_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include "cmplx.h"

//Biliothek für mathematische Funktionen
#define SWAP(a,b) tempr =(a);(a) = (b); (b) = tempr
#define QUAT(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
void rDFT_Algo(Cmplx  pOut[], double pIn[], int nDim); //Diskrete Fourieranalyse zum Test
void cFFT_Algo(double data[], unsigned long nn, int isign);//FFT-Algorithmus(complex Funktion)
void rFFT_Algo(double data[], unsigned long nn, int isign);//FFT-Algorithmus(reelle Funktion)
bool isEven(int n);
double intXuToXo( double f(double, double, double ),
                double a, double b, double xU, double xO);
double intCoshCos(double a, double b, double x);
double intCoshCos(double a, double x);
double intCoshSin(double a, double b, double x);
double intCoshSin(double a, double x);
double intSinhCos(double a, double b, double x);
double intSinhCos(double a, double x);
double intSinhSin(double a, double b, double x);
double intSinhSin(double a, double x);
double intCosCos(double a, double b, double x);
double intCosCos(double a, double x);
double intSinSin(double a, double b, double x);
double intSinSin(double a, double x);
double intSinCos(double a, double b, double x);
double intSinCos(double a, double x);
double intCosSin(double a, double x);
double intCosSin(double a, double b, double x);
Cmplx  ccotbess(Cmplx x, Cmplx y);
double cotbess(double x, double y);
Cmplx  hank_1(int n, double x);
Cmplx  hank_2(int n, double x);
double bessJ0(double x);
double bessJ1(double x);
double bessJn(int n,double x);
Cmplx  cBessJ0(Cmplx z);
Cmplx  cBessJ1(Cmplx z);

double bessY0(double x);
double bessY1(double x);
double bessYn(int n, double x);
Cmplx  cBessY0(Cmplx z);
Cmplx  cBessY1(Cmplx z);

double bessI0(double x);
double bessI1(double x);
double bessIn(int n, double x);

double bessK0(double x);
double bessK1(double x);
double bessKn(int n, double x);

double  dBToMag(double x);//dB in mag-Umwandeln
double  coth(double x);
double  acot(double x);//Umkehrfunktion von Cot;
double  getPi();
double  elip_K(double x);
double  asym_Jn(int n, double x );
double  fakultaet(int n);
double  ipotnz( double, int );
double  potnz( double, double );
double  nRoot(double, int); //Berechnet n-te Wurzel
double  invSinh(double x );
double  invCosh(double x );
double  invTanh(double x );//Umkehrfunktion vom Tanh
double  invCoth(double x );//Umkehrfunktion vom Coth
//double  exp10( double a );
#endif