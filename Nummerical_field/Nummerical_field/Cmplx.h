/***************************************************************************/
/*                                                                         */
/*  Header-Datei zu Cmplx.cpp                                              */
/*  G. Zimmer                                                              */
/*                                                                         */
/***************************************************************************/

#ifndef CMPLX_H
#define CMPLX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <iostream>

//#ifndef KONSOL_APP
//#include <afxwin.h>
//#endif

// Makrodefinitionen:
#define Z_MIN    1.E-36
#define Z_MAX    1.E+36
#define REAL     double
#define J        Cmplx(0.,1)

class Cmplx 
{
//DECLARE_SERIAL(Cmplx)
protected:
  REAL a, b;

public:
	 Cmplx(REAL x=0, REAL y=0);
     Cmplx(int  aInt);
     Cmplx(const Cmplx &aZ);
     Cmplx(const Cmplx *pZ);
	 Cmplx invers(void)const;
	 Cmplx operator=(const Cmplx& aZ);
	 Cmplx operator+=(const Cmplx& z);
	 Cmplx operator+ (const Cmplx& z);
	 Cmplx operator-=(const Cmplx& z);
	 Cmplx operator- (const Cmplx& z);
	 Cmplx operator*=(const Cmplx& z);
	 Cmplx operator* (const Cmplx& z);
	 Cmplx operator/=(const Cmplx& z);
	 Cmplx operator/ (const Cmplx& z);
     int   operator!=(const Cmplx& z);
     int   operator==(const Cmplx& z);
     Cmplx operator-()const;
     Cmplx operator*=(const REAL aR);
	 Cmplx operator* (const REAL aR);
     Cmplx operator*=(const int  aI);
	 Cmplx operator* (const int  aI);
	 Cmplx operator/=(const REAL aR);
	 Cmplx operator/ (const REAL aR);
     Cmplx operator/=(const int  aI);
	 Cmplx operator/ (const int  aI);
	 Cmplx operator= (const REAL z);
	 Cmplx operator= (const int aI);
	 Cmplx operator+=(const REAL aR);
	 Cmplx operator+ (const REAL aR);
     Cmplx operator+=(const int  aI);
	 Cmplx operator+ (const int  aI);
	 Cmplx operator-=(const REAL aR);
	 Cmplx operator- (const REAL aR);
     Cmplx operator-=(const int  aI);
	 Cmplx operator- (const int  aI);
	 Cmplx conjg(void) const;
	 REAL  abs()const;
	 REAL  winkel(void)const;
	 REAL  MagTodB();// Umwandlung in dB
	 REAL  real()const
		{return(a);};
	 REAL imag()const
		{return(b);};
     bool IsReal(); //Nur Reell
     bool IsImag(); //Nur Imaginär
	 void  show(void);
     void  SetRealPart(REAL aVal){a = aVal;};
     void  SetImagPart(REAL aVal){b = aVal;};
#ifndef MY_KONSOL_APP
	 //void  Serialize(CArchive &ar);
#endif
	 friend Cmplx polarToRect(REAL betrag, REAL winkel);
	 friend Cmplx cmplx(REAL a, REAL b);
	 friend Cmplx csqrt(const Cmplx& z);
	 friend REAL abs(Cmplx z)
		{return(sqrt(z.a*z.a + z.b*z.b));};
	 friend Cmplx exp (const Cmplx& z);
	 friend Cmplx sinh(const Cmplx& z);
	 friend Cmplx cosh(const Cmplx& z);
	 friend Cmplx cos (const Cmplx& z);
	 friend Cmplx sin (const Cmplx& z);
	 friend Cmplx tanh(const Cmplx& z)
		{return(sinh(z)/cosh(z));};
	 friend Cmplx tan (const Cmplx& z)
		{return(sin(z)/cos(z));};
	 friend int gauss(Cmplx *z, int n);
	 friend int gauss(Cmplx *ma, int dim, int nSol);
	 friend int gauss(Cmplx *out_vek, Cmplx *in_ma, Cmplx *in_vek,
			  int dim, int maxDim);
	 friend int gauss(Cmplx *loe_ma, Cmplx *in_ma, Cmplx *er_ma,
			  int dim, int nSol,int maxDim);
};

// Ausgabe zu einem ostream
std::ostream& operator<<(std::ostream& os, const Cmplx& a);

// Eingabe von einem istream
std::istream& operator>>(std::istream& is,  Cmplx& a);

#endif
