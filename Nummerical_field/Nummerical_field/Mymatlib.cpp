#include "stdafx.h"
#include <math.h>
#include<stdlib.h>

#ifndef MYMATLIB_H
#include "mymatlib.h"
#endif

#include "physicon.h"
//--------------------------------------------------
double dBToMag(double x)
{
	x /= 20.;

	return potnz(10.,x);
}
//------------------------------------------------------------
void rDFT_Algo(Cmplx* pOut, double* pIn, int nDim)
// 
//  Realisiert diskrete Fourieranalyse zu Test Zwecken
//  
//  pOut[] - Komplexes Output Array der Dimension nDim
//  pIn[]  - reelles Input Array der Dimension    nDim
//  nDim   - Anzahl der Abtastwerte (Linksseitige Abtastung)
//           Gleichung realisiert nach EG-Vorlesung DF-3
{
	Cmplx  sum;
	double pi2  = 2.*PI;
	double quot = 1./((double) nDim);
	for(int m = 0; m < nDim; m++)
	{
		sum         = 0;
		Cmplx alpha = -J*pi2*m*quot;
		for(int n = 0; n < nDim; n++)
		{
			//sum += exp(alpha*(n+1))*pIn[n]*quot;
			sum += exp(alpha*n)*pIn[n]*quot;
		}
		pOut[m] = sum;
	}
}
//------------------------------------------------------------
void rFFT_Algo(double data[], unsigned long n, int isign)
//
//  FFT-Algorithmus aus Numerical Recipies in C, Cambridge
//  University Press Seite 513
//  Ersetzt data[1 ...n] durch diskrete Fouriertransformierte
//                       komplexe Werte der positiven Frequenzen
//                       (n/2)  
//  data  - reelle Abtastwerte der Dimension n
//  n     -  muss 2er Potenz sein wird nicht überprüft 
//  isign - Vorzeichen der Exponentialfunktion
{
	unsigned long i, i1,i2,i3,i4,np3;
	double   c1 = 0.5,c2,h1r,h1i,h2r,h2i;
	double   wr,wi,wpr,wpi,wtemp,theta;

	theta = 3.141592653589793/(double)(n>>1);
	if(isign == 1)
	{
		c2 = -0.5;
		cFFT_Algo(data,n>>1,1);
	}
	else
	{
		c2    = 0.5;
		theta = -theta;
	}

	wtemp = sin(0.5*theta);
	wpr   = -2.0*wtemp*wtemp;
	wpi   = sin(theta);
	wr    = 1.0 + wpr;
	wi    = wpi;
	np3   = n + 3;
	for(i = 2; i <= (n>>2); i++)
	{
		i4  =  1 + (i3=np3-(i2=1+(i1=i+i-1)));
		h1r =  c1*(data[i1] + data[i3]);
		h1i =  c1*(data[i2] - data[i4]);
		h2r = -c2*(data[i2] + data[i4]);
		h2i =  c2*(data[i1] - data[i3]);

		data[i1] =  h1r + wr*h2r - wi*h2i;
		data[i2] =  h1i + wr*h2i + wi*h2r;
		data[i3] =  h1r - wr*h2r + wi*h2i;
		data[i4] = -h1i + wr*h2i + wi*h2r;

		wr = (wtemp=wr)*wpr - wi*wpi + wr;
		wi =  wi*wpr + wtemp*wpi + wi;
	}
	if(isign == 1)
	{
		data[1] = (h1r=data[1]) + data[2];
		data[2] = h1r - data[2];
	}
	else
	{
		data[1] = c1*((h1r=data[1]) + data[2]);
		data[2] = c1*(h1r-data[2]);
		cFFT_Algo(data,n>>1,-1);
	}
}
//------------------------------------------------------------
void cFFT_Algo(double data[], unsigned long nn, int isign)
//
//  FFT-Algorithmus aus Numerical Recipies in C, Cambridge
//  University Press Seite 507
//  Ersetzt data[1 ...2*nn] durch diskrete Fouriertransformierte
//  data  - complexes Array der Dimension nn
//  nn    -  muss 2er Potenz sein wird nicht überprüft 
//  isign - Vorzeichen der Exponentialfunktion
{
	unsigned long n, mmax, m, j, istep, i;
	double   wtemp, wr, wpr, wpi, wi, theta;
	double   tempr, tempi;

	//Bit-reversal Section
	n = nn << 1;
	j = 1;
	for(i = 1; i < n; i += 2)
	{
		if(j > i)
		{
			SWAP(data[j]    ,data[i]);
			SWAP(data[j + 1],data[i + 1]);
		}
		m = n >> 1;
		while( m >= 2 && j > m)
		{
			j -= m;
			m >>=1;
		}
		j += m;
	}
	//Ab hier Danielson-Lancos 
	mmax = 2;
	while( n > mmax)
	{
		istep = mmax << 1;
		theta = -isign*(PI*2./mmax);
		wtemp = sin(0.5*theta);
		wpr   = -2.0*wtemp*wtemp;
		wpi   = sin(theta);
		wr    = 1.0;
		wi    = 0.0;
		for( m = 1; m < mmax; m += 2)
		{
			for( i = m; i <= n; i+= istep)
			{
				j          = i + mmax;
				tempr      = wr*data[j]   - wi*data[j + 1];
				tempi      = wr*data[j+1] + wi*data[j];
				data[j]    = data[i]   - tempr;
				data[j+1]  = data[i+1] - tempi;
				data[i]   += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp=wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
}
//------------------------------------------------------------
bool isEven(int n)
{
  if(n == 0) return true;
  if(n%2)
       return true;
  else
       return false;
}
//------------------------------------------------------------
//  Bestimmtes Integral in den Grenzen xU bis xO
//------------------------------------------------------------
double intXuToXo(double f(double,double,double),double a, double b,double xU, double xO)
{
  return f(a,b,xO) - f(a,b,xU);
}
//------------------------------------------------------------
// Unbestimmtes Int(Cosh(ax)*sin(bx) dx
//------------------------------------------------------------
double intCoshSin(double a, double b, double x)
{
   if( a == b) return intCoshCos(a,x);

   double val1  = exp(a*x)/(a*a + b*b);
          val1 *= (-b*cos(b*x) + a*sin(b*x));

   double val2  = exp(-a*x)/(a*a + b*b);
          val2 *= (-b*cos(b*x) - a*sin(b*x));

   return 0.5*(val1 + val2);
}
//------------------------------------------------------------
double intCoshSin(double a, double x)
{
    if( fabs(a) < Z_MIN ) return 0;

    double res1  = 0.25*exp(a*x)/a;
           res1 *= (-cos(a*x) + sin(a*x));

    double res2  = 0.25*exp(-a*x)/a;
           res2 *= (-sin(a*x) - cos(a*x));

    return res1 + res2;
}
//------------------------------------------------------------
// Unbestimmtes Int(Cosh(ax)*cos(bx) dx
//------------------------------------------------------------
double intCoshCos(double a, double b, double x)
{
   if( a == b) return intCoshCos(a,x);

   double val1  = exp(a*x)/(a*a + b*b);
          val1 *= (a*cos(b*x) + b*sin(b*x));

   double val2  = exp(-a*x)/(a*a + b*b);
          val2 *= (-a*cos(b*x) + b*sin(b*x));

   return 0.5*(val1 + val2);
}
//------------------------------------------------------------
double intCoshCos(double a, double x)
{
    if( fabs(a) < Z_MIN ) return x;

    double res1  = 0.25*exp(a*x)/a;
           res1 *= (cos(a*x) + sin(a*x));

    double res2  = 0.25*exp(-a*x)/a;
           res2 *= (sin(a*x) - cos(a*x));

    return res1 + res2;
}
//------------------------------------------------------------
// Unbestimmtes Int(sinh(ax)*cos(bx) dx
//------------------------------------------------------------
double intSinhCos(double a, double b, double x)
{
   if( a == b) return intSinhCos(a,x);

   double val1  = exp(a*x)/(a*a + b*b);
          val1 *= (a*cos(b*x) + b*sin(b*x));

   double val2  = exp(-a*x)/(a*a + b*b);
          val2 *= (-a*cos(b*x) + b*sin(b*x));

   return 0.5*(val1 - val2);
}
//------------------------------------------------------------
double intSinhCos(double a, double x)
{
    if( fabs(a) < Z_MIN ) return 0;

    double res1  = 0.25*exp(a*x)/a;
           res1 *= (cos(a*x) + sin(a*x));

    double res2  = 0.25*exp(-a*x)/a;
           res2 *= (cos(a*x) - sin(a*x));

    return res1 + res2;
}
//------------------------------------------------------------
// Unbestimmtes Int(sinh(ax)*sin(bx) dx
//------------------------------------------------------------
double intSinhSin(double a, double b, double x)
{
   if( a == b) return intSinhSin(a,x);

   double val1  = exp(a*x)/(a*a + b*b);
          val1 *= (a*sin(b*x) - b*cos(b*x));

   double val2  = exp(-a*x)/(a*a + b*b);
          val2 *= (-a*sin(b*x) - b*cos(b*x));

   return 0.5*(val1 - val2);
}
//------------------------------------------------------------
double intSinhSin(double a, double x)
{
    if( fabs(a) < Z_MIN ) return 0;

    double res1  = 0.25*exp(a*x)/a;
           res1 *= (-cos(a*x) + sin(a*x));

    double res2  = 0.25*exp(-a*x)/a;
           res2 *= (cos(a*x) + sin(a*x));

    return res1 + res2;
}
//------------------------------------------------------------
//  Unbestimmtes Int(cos(ax)*cos(bx) dx
//------------------------------------------------------------
double intCosCos(double a, double b, double x)
//  a  - Konstante des Sinus
//  b  - Konstante des Sinus
{
   if( fabs(a) == fabs(b) )
      {
        if( a == b)   return intCosCos(a,x);

        if( a < 0)    return intCosCos(-a,x);

        if( b < 0)    return intCosCos(a,x);
      }
   double res1  = 0.5*sin((b-a)*x);
          res1 /=(b-a);

   double res2  = 0.5*sin((a + b)*x);
          res2 /=(a + b);

   return res1 + res2;
}
//------------------------------------------------------------
double intCosCos(double a, double x)
{
   if(fabs(a) < Z_MIN) return x;

   double res  = a*x;
          res += sin(a*x)*cos(a*x);

   return 0.5*res/a;
}
//------------------------------------------------------------
//  Unbestimmtes Int(sin(ax)*sin(bx) dx
//------------------------------------------------------------
double intSinSin(double a, double b, double x)
//  a  - Konstante des Sinus
//  b  - Konstante des Sinus
{
   if( fabs(a) == fabs(b) )
      {
        if( a == b)   return intSinSin(a,x);

        if( a < 0)    return intSinSin(-a,x);

        if( b < 0)    return intSinCos(a,x);
      }
   double res1  = 0.5*sin((b-a)*x);
          res1 /=(b-a);

   double res2  = -0.5*sin((a + b)*x);
          res2 /=(a + b);

   return res1 + res2;
}
//------------------------------------------------------------
double intSinSin(double a, double x)
{
   if(fabs(a) < Z_MIN) return 0;

   double res  = a*x;
          res -= sin(a*x)*cos(a*x);

   return 0.5*res/a;
}
//------------------------------------------------------------
//  Unbestimmtes Int(sin(ax)*cos(bx) dx
//------------------------------------------------------------
double intSinCos(double a, double b, double x)
//  a  - Konstante des Sinus
//  b  - Konstante des Cosinus
{
   if( fabs(a) == fabs(b) )
      {
        if( a == b)   return intSinCos(a,x);

        if( a < 0)    return -intSinCos(a,x);

        if( b < 0)    return  intSinCos(a,x);
      }
   double res1  = -0.5*cos((a+b)*x);
          res1 /=(a+b);

   double res2  = 0.5*cos((b - a)*x);
          res2 /=(b - a);

   return res1 + res2;
}
//-------------------------------------------------------------
double intCosSin(double a, double x)
{
    return intSinCos(a,x);
}
//-------------------------------------------------------------
double intCosSin(double a, double b, double x)
{
    return intSinCos(b,a,x);
}
//------------------------------------------------------------
double intSinCos(double a, double x)
{
   if(fabs(a) < Z_MIN) return 0;

   return 0.5*QUAT( sin(a*x))/a;
}
//------------------------------------------------------------
// Funktion coth(x)
//------------------------------------------------------------
double coth(double x)
{
   double nenr = tanh(x);

   if(fabs(nenr) < Z_MIN) return Z_MAX;
   else
      return 1./nenr;
}
//------------------------------------------------------------
//   Large Radial Cotangent Function
//   siehe: Transmission Line Design Handbook p. 303
//------------------------------------------------------------
Cmplx ccotbess(Cmplx x, Cmplx y)
{
    Cmplx zaelr = cBessY0(x)*cBessJ1(y) - cBessJ0(x)*cBessY1(y);
    Cmplx nenr  = cBessJ1(x)*cBessY1(y) - cBessY1(x)*cBessJ1(y);
    if( nenr.abs() < Z_MIN) nenr = Z_MIN;

    return zaelr/nenr;
}
//------------------------------------------------------------
double cotbess(double x, double y)
{
    double zaelr = bessY0(x)*bessJ1(y) - bessJ0(x)*bessY1(y);
    double nenr  = bessJ1(x)*bessY1(y) - bessY1(x)*bessJ1(y);
    if(fabs(nenr) < Z_MIN) nenr = Z_MIN;

    return zaelr/nenr;
}
//------------------------------------------------------------
//   Hankelfunktion Hn_2(x) = Jn(x) - j*Yn(x)
//------------------------------------------------------------
Cmplx hank_2(int n, double x)
{
    Cmplx res  = bessJn(n,x);
          res -= J*bessYn(n,x);

    return res;
}
//------------------------------------------------------------
//   Hankelfunktion Hn_1(x) = Jn(x) + j*Yn(x)
//------------------------------------------------------------
Cmplx hank_1(int n, double x)
{
    Cmplx res  = bessJn(n,x);
          res += J*bessYn(n,x);

    return res;
}
//------------------------------------------------------------
// Literaturstelle zu Besselfunktionen:
//
//       William, H. Press, Brian P. Flannery
//       Saul A. Teukolsky, William T. Vetterling
//
//       Numerical Recipes in Pascal
//       Cambridge University Press
//------------------------------------------------------------
//          Besselfunktion Kn(x)
//------------------------------------------------------------
double bessKn(int n, double x)
{
  if(n == 0)
      return bessK0(x);
  else if(n == 1)
      return bessK1(x);
  else
   {
     double tox,bkp,bkm,bk;
     int j;

     tox = 2.0/x;
     bkm = bessK0(x);
     bk  = bessK1(x);

     for( j = 1; j <= n-1; j++)
	{
	   bkp = bkm + j*tox*bk;
	   bkm = bk;
	   bk  = bkp;
	}
     return bk;
   }
}
//------------------------------------------------------------
//          Besselfunktion K1(x)
//------------------------------------------------------------
double bessK1(double x)
{
   double y, ans;

   if( x <= 2.0 )
     {
	y    = x*x/4.0;

	ans  = (log(x/2.0)*bessI1(x))+ (1.0/x)*(1.0+y*(0.154431144+y*(-0.67278579
		+y*(-0.18156897 + y*(-0.1919402e-1+y*(-0.110404e-2
		+y*(-0.4866e-4)))))));
     }
   else
     {
	y  = 2.0/x;

	ans = (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1
	      +y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
     }
  return ans;
}
//------------------------------------------------------------
//          Besselfunktion K0(x)
//------------------------------------------------------------
double bessK0(double x)
{
   double y, ans;

   if( x <= 2.0 )
     {
	y    = x*x/4.0;

	ans  = (-log(x/2.0)*bessI0(x)) + (-0.57721566 + y*(0.42278420
		+y*(0.23069756 + y*(0.3488590e-1+y*(0.262698e-2
		+y*(0.10750e-3+y*0.74e-5))))));
     }
   else
     {
	y  = 2.0/x;

	ans = (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1
	      +y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
     }
  return ans;
}
//------------------------------------------------------------
//          Besselfunktion In(x)
//------------------------------------------------------------
double bessIn(int n, double x)
{
  const int     iacc  = 40;
  const double  bigno = 1.0e+10;
  const double  bigni = 1.0e-10;

  double bi, bim, bip, tox, ans;
  int    j, m;

  if(n == 0)
     return bessI0(x);
  else if(n == 1)
     return bessI1(x);
  else
   {
      if(x == 0.0) return 0.0;

      ans  = 0.0;
      tox  = 2.0/fabs(x);
      bip  = 0.0;
      bi   = 1.0;

      div_t h = div(n+int(sqrt(1.0*(iacc*n))),2);
      m       = 2*h.quot;

      for( j = m; j >= 1; j--)
	{
	   bim = bip + j*tox*bi;
	   bip = bi;
	   bi  = bim;
	   if(fabs(bi) > bigno)
	     {
		ans = ans*bigni;
		bi  = bi*bigni;
		bip = bip*bigni;
	     }
	   if( j == n) ans = bip;
	}
    if ( (x < 0) && (fmod((double)n,2.0) != 0) ) ans = -ans;
    return ans*bessI0(x)/bi;
   }
}
//------------------------------------------------------------
//          Besselfunktion I1(x)
//------------------------------------------------------------
double bessI1(double x)
{
  double ax;
  double y, ans;
  if (fabs(x) < 3.75)
     {  //Approximation ber Polynom
	y    = pow((x/3.75),2);
	ans = x*(0.5 + y*(0.87890594+y*(0.51498869 +y*(0.15084934
	       +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));

       return ans;
     }
  else
     {    //Approximation ber Polynom
	  ax = fabs(x);
	  y  = 3.75/ax;

	  ans = 0.2282967e-1 +y*(-0.2895312e-1 +y*(0.1787654e-1 -y*0.420059e-2));
	  ans = 0.39894228 + y*(-0.3988024e-1 +y*(-0.362018e-2
	      + y*(0.163801e-2+y*(-0.1031555e-1 +y*ans))));

	  ans = (exp(ax)/sqrt(ax))*ans;
	  if( x < 0.0) ans = -ans;
      return ans;
  }
}
//------------------------------------------------------------
//          Besselfunktion I0(x)
//------------------------------------------------------------
double bessI0(double x)
{
  double ax;
  double y, ans;
  if (fabs(x) < 3.75)
     {  //Approximation ber Polynom
	y    = pow((x/3.75),2);
	ans = 1.0+y*(3.5156229+y*(3.0899424 +y*(1.2067492
	      +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));

       return ans;
     }
  else
     {    //Approximation ber Polynom
	  ax = fabs(x);
	  y  = 3.75/ax;

	  ans=(exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 +y*(0.225319e-2
	      + y*(-0.157565e-2+y*(0.916281e-2 +y*(-0.2057706e-1 + y*(0.2635537e-1
	      + y*(-0.1647633e-1 + y*0.392377e-2))))))));

      return ans;
  }
}
//------------------------------------------------------------
//          Besselfunktion Yn(x)
//------------------------------------------------------------
double bessYn(int n, double x)
{
  double by, bym, byp, tox;
  int j;

  if(n == 0)
     return bessY0(x);
  else if(n == 1)
     return bessY1(x);
  else
   {
      tox = 2.0/x;
      by  = bessY1(x);
      bym = bessY0(x);
      for( j = 1; j <= n-1; j++)
	{
	  byp = j*tox*by - bym;
	  bym = by;
	  by  = byp;
	}
      return by;
   }
}
//------------------------------------------------------------
//          Besselfunktion Y1(x)
//------------------------------------------------------------
Cmplx cBessY1(Cmplx z) //Näherung für komplexe Besselfunktion Y1
{                      //     Imagiärteil << Realteil
   REAL a = z.real();
   REAL b = z.imag();

   Cmplx res = bessY1(a);

   res +=  J*b*(bessY0(a) - bessY1(a)/a );

   return res;
}
//------------------------------------------------------------
double bessY1(double x)
{
  double xx, z;
  double y, ans, ans1, ans2;
  if ( x < 8.0)
     {  //Approximation ber Polynom
	y    = pow(x,2);
	ans1 = x*(-0.4900604943e13+y*(0.1275274390e13+y*(-0.5153438139e11
	      +y*(0.7349264551e9+y*(-0.4237922726e7+y*0.8511937935e4)))));

	ans2 = 0.249958057e14+y*(0.424419664e12+y*(0.3733650367e10
	      +y*(0.2245904002e8+y*(0.1020426050e6
	      +y*(0.3549632885e3 + y*1.0)))));

	ans  = ans1/ans2 + 0.636619772*(bessJ1(x)*log(x) -1./x);
       return ans;
     }
  else
     {    //Approximation ber cos sin Summe
	  z  = 8.0/x;
	  y  = pow(z,2);
	  xx = x-2.356194491;

	 ans1 = 1.0+y*(0.183105e-2+y*(-0.351639646e-4
		+y*(0.2457520174e-5+y*(-0.240337019e-6))));

	 ans2 = 0.04687499995+y*(-0.2002690873e-3
	     +y*(0.8449199096e-5+y*(-0.8822987e-6
	     -y*0.105787412e-6)));

	 ans = sqrt(0.636619772/x)*(sin(xx)*ans1 + z*cos(xx)*ans2);

      return ans;
  }
}
//------------------------------------------------------------
//          Besselfunktion Y0(x)
//------------------------------------------------------------
Cmplx cBessY0(Cmplx z) //Näherung für komplexe Besselfunktion Y0
{                      //     Imagiärteil << Realteil
   REAL a = z.real();
   REAL b = z.imag();

   Cmplx res = bessY0(a);

   res -=  J*b*bessY1(a);

   return res;
}
//-------------------------------------------------------------
double bessY0(double x)
{
  double xx, z;
  double y, ans, ans1, ans2;
  if (fabs(x) < 8.0)
     {  //Approximation ber Polynom
	y    = pow(x,2);
	ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
	      +y*(10879881.29+y*(-86327.92757+y*228.4622733))));

	ans2 = 40076544269.0+y*(745249964.8+y*(7189466.438
	      +y*(47447.26470+y*(226.1030244+y*1.0))));

	ans  = ans1/ans2 + 0.636619772*bessJ0(x)*log(x);
       return ans;
     }
  else
     {    //Approximation ber cos sin Summe
	  z  = 8.0/x;
	  y  = pow(z,2);
	  xx = x-0.785398164;

	 ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
		+y*(-0.2073370639e-5+y*0.2093887211e-6)));

	 ans2 = -0.1562499995e-1+y*(0.1430488765e-3
	     +y*(-0.6911147651e-5+y*(0.7621095161e-6
	     -y*0.934945152e-7)));

	 ans = sin(xx)*ans1 + z*cos(xx)*ans2;

	 ans = sqrt(0.636619772/x)*ans;

      return ans;
  }
}
//------------------------------------------------------------
//          Besselfunktion J0(x)
//------------------------------------------------------------
Cmplx cBessJ0(Cmplx z) //Näherung für komplexe Besselfunktion J0
{                      //     Imagiärteil << Realteil
   REAL a = z.real();
   REAL b = z.imag();

   Cmplx res = bessJ0(a);

   res -=  J*b*bessJ1(a);

   return res;
}
//--------------------------------------------------------------
double bessJ0(double x)
{
  double ax, xx,z;
  double y,ans,ans1,ans2;

  if (fabs(x) < 8.0)
     {  //Approximation ber Polynom
	y    = x*x;
	ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7
	      +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));

	ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718
	      +y*(59272.64853+y*(267.8532712+y*1.0))));
       return ans1/ans2;
     }
  else
     {    //Approximation ber cos sin Summe
	  ax = fabs(x);
	  z  = 8.0/ax;
	  y  = z*z;
	  xx = ax-0.785398164;

	 ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
		+y*(-0.2073370639e-5+y*0.2093887211e-6)));

	 ans2 = -0.1562499995e-1+y*(0.1430488765e-3
	     +y*(-0.6911147651e-5+y*(0.7621095161e-6
	     -y*0.934945152e-7)));

	 ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);

      return ans;
  }
}
//----------------------------------------------------------------
//            Besselfunktion J1(x)
//----------------------------------------------------------------
Cmplx cBessJ1(Cmplx z) //Näherung für komplexe Besselfunktion J1
{                      //     Imagiärteil << Realteil
   REAL a = z.real();
   REAL b = z.imag();

   Cmplx res = bessJ1(a);

   res +=  J*b*(bessJ0(a) - bessJ1(a)/a );

   return res;
}
//------------------------------------------------------------------
double bessJ1(double x)
{
   double ax,xx,z;
   double y,ans,ans1,ans2;

  if (fabs(x) < 8.0)
   {               //Approximation durch Polynom
      y = pow(x,2);

      ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
	     +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));

      ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74
	     +y*(99447.43394+y*(376.9991397+y*1.0))));

      return ans1/ans2;
   }
   else
   {
      ax = fabs(x);
      z  = 8.0/ax;
      y  = pow(z,2);
      xx = ax-2.356194491;

      ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
	     +y*(0.2457520174e-5+y*(-0.240337019e-6))));

      ans2 = 0.04687499995+y*(-0.2002690873e-3
	     +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));

      ans = sqrt(0.636619772/ax)*(cos(xx)*ans1
	    -z*sin(xx)*ans2);
      if(x < 0.0) ans = -ans;
      return ans;
    }
}
//--------------------------------------------------------------
//            Besselfunktion Jn(x)
//--------------------------------------------------------------
double bessJn(int n, double x)
{

const int    iacc  = 40;
const double bigno = 1.0e+10;
const double bigni = 1.0e-10;

double bj,bjm,bjp,sum,tox,ans;
int    j,jsum,m;

if(n == 0)
     return bessJ0(x);
else if( n == 1)
     return  bessJ1(x);
else
  {
    if( x == 0.0) ans = 0.0;
    else if (fabs(x) > n)
      {
	tox = 2./fabs(x);
	bjm = bessJ0(fabs(x));
	bj  = bessJ1(fabs(x));

	//Rekursion
	for (j=1;j <= n-1;j++)
	{
	   bjp = j*tox*bj-bjm;
	   bjm = bj;
	   bj = bjp;
	}
	ans = bj;
      }
      else
      //Rekursion mit Normalisierung
      {
	tox     = 2.0/fabs(x);
	div_t h = div(n+int(sqrt(1.0*(iacc*n))),2);
	m       = 2*h.quot;

	ans  = 0.0;
	jsum = 0;
	sum  = 0.0;
	bjp  = 0.0;
	bj   = 1.0;

	for (j = m; j>=1; j--)
	{
	  bjm = j*tox*bj-bjp;
	  bjp = bj;
	  bj  = bjm;

	  if (fabs(bj) > bigno)
	  {
	    bj  = bj*bigni;
	    bjp = bjp*bigni;
	    ans = ans*bigni;
	    sum = sum*bigni;
	  }

	  if (jsum != 0) sum = sum + bj;

	  jsum = 1-jsum;

	  if (j == n) ans = bjp;
	}

	sum = 2.0*sum-bj;
	ans = ans/sum;
      }

    if ((x < 0) && (fmod((double)n,2.0) != 0)) ans = -ans;
    return ans;
  }
}

//--------------------------------------------------------------------------
double invSinh(double x)
{
  return log(x + sqrt(x*x + 1) );
}
//--------------------------------------------------------------------------
double invCosh(double x)
{
  return log( x + sqrt(x*x - 1) );
}
//--------------------------------------------------------------------------
double invTanh(double x)
{
	double aVal = (1 + x )/(1 - x );

	return 0.5*log(aVal);
}
//--------------------------------------------------------------------------
double invCoth(double x)
{
	double aVal = (x + 1)/(x - 1);

	return 0.5*log(aVal);
}
//--------------------------------------------------------------------------
double acot(double x)
{
	return 0.5*PI - atan(x);
}
//--------------------------------------------------------------------------
double nRoot(double x, int n) //Berechnet n-teWurzel
{
  if(x < 0 ) x = -1*x;
  if(n < 0 ) n = -1*n;
  if(n == 0) return x;
  return exp(log(x)/n);
}
//--------------------------------------------------------------------------
double getPi()
{
  return 4*atan(1.);
}
 //-------------------------------------------------------------------------
 // Name     : elipfktK
 // Bemerkung: Berechnet vollständiges elliptisches Integral 1. Art
 // x        : Argument |x| < 1
 //-------------------------------------------------------------------------
 double elip_K(double k)
 {  //Berechnung nach Reihendarstellung Bronstein Seite 124
    double epsilon = 1.E-10;
    if(fabs(k) >= 1.) return 1.E+30;

    //Abbruch Index bestimmen
    double nTimes = log(epsilon*(1. -k*k))/log(k*k);
    if(nTimes > 2147483647) return 1.E+30;
    long int i   = 1;
    double aI    = 0.5*k;
    double sum   = aI*aI;
    while(i < nTimes)
      {
	 i++;
	 aI  *= (2*i - 1)*k/(2*i);
	 sum += aI*aI;
      }
    sum += 1;
    sum *= 0.5*PI;
    return sum;
 }
//--------------------------------------------------------------------------
//Asymptotische Darstellung der Besselfunktion Jn
//---------------------------------------------------------------------------
double asym_Jn(int n, double x)
{
   double arg = x - 0.25*PI - 0.5*PI*n;
   double fak = sqrt(2./(PI*x));
   return fak*cos(arg);
}
//--------------------------------------------------------------------------
double fakultaet(int n)
{
   if(n == 0) return 1.;
   if(n == 1) return 1.;
   double aValue = 1.;
   for(int i = 2; i <= n; i++)
      {
	  aValue *= i;
      }
   return aValue;
}
/******************************-Funktion-***********************************/
/*  Name     : ipotnz                                                      */
/*  Parameter: siehe unten                                                 */
/*  Bemerkung: Berechnet ganzzahlige Potenzen                              */
/***************************************************************************/
/*  x   Von diesem Wert wird die Potenz berechnet  */
/*  n   Ganzzahliger Exponent                      */
/***************************************************/
double ipotnz(double x, int n)
{
	 double z;
	 int sgn=0, i;
	 if( n==0 )  return(1.);
	 if( n<0 )
		{
		  n =  -n;
		  sgn = 1;
		}
	 z = x;
	 for( i=1 ; i<n ; i++ )
		{
		  z *= x;
		}
	 if(sgn)
		return(1./z);
	 else
		return(z);
}

/******************************-Funktion-***********************************/
/*  Name     : potnz                                                       */
/*  Parameter: siehe unten                                                 */
/*  Bemerkung: Berechnet rationale Potenzen einer Zahl x                   */
/***************************************************************************/
/*  x  Von diesem Wert wird die Potenz berechnet  */
/*  z  Exponent                                   */
/**************************************************/
double potnz( double x, double z )
{
	return( exp( z*log(x) ) );
}

