/***************************************************************************/
/*                                                                         */
/*  CMPLXBIB.CPP Komplexe Zahlen               */
/*  G. Zimmer                                                              */
/*                                                                         */
/***************************************************************************/
//#ifndef KONSOL_APP
// #include <afx.h>
//#endif
#include "stdafx.h"

#ifndef  CMPLX_H
#include "cmplx.h"
#endif

#include<float.h>
#include "physicon.h"

// Definition komplexe Eingabe re, im:
std::istream& operator >> ( std::istream& is, Cmplx& z )
{
	char c;
	REAL a, b;

	is >> c >> a >> c >> b>>c;
	z = cmplx(a,b);
	return is;
}
/*-----------------------------------------------*/
// Definition komplexe Ausgabe
std::ostream& operator << ( std::ostream& os, const Cmplx& a )
{
	os << '(' << a.real() << ',' << a.imag() << ')';
	return os;
}
#ifndef MY_KONSOL_APP
//-------------------------------------------------
//IMPLEMENT_SERIAL( Cmplx, CObject, 1 )
//-------------------------------------------------
void Cmplx::Serialize(CArchive &ar)
{
	if(ar.IsStoring())
	{
		ar << a;
		ar << b;
	}
	else
	{
		ar >> a;
		ar >> b;
	}
}
#endif
//------------------------------------------------//
REAL  Cmplx::abs()const
{
   if(fabs(a) == 0 && fabs(b) == 0) return 0.;
   if(fabs(a) >= fabs(b))
       return fabs(a)*sqrt(1. + QUAT(b/a));
   else
       return fabs(b)*sqrt(1. + QUAT(a/b));
   //return sqrt( a*a + b*b);
}
/*-----------------------------------------------*/
Cmplx polarToRect(REAL betrag, REAL winkel)
{
  const REAL kw = 0.0174533;
       Cmplx z;
       z.a  = betrag*cos(kw*winkel);
       z.b  = betrag*sin(kw*winkel);
       return z;
}
/*-----------------------------------------------*/
Cmplx cmplx(REAL a, REAL b)
{
	Cmplx z;

	z.a = a;
	z.b = b;
	return z;
}
/*-----------------------------------------------*/
Cmplx csqrt(const Cmplx& z)
{
   Cmplx wurzel;
   REAL betrag = sqrt(z.abs());
   REAL wi     = z.winkel()/2.;
   wurzel      = polarToRect(betrag,wi);
   return wurzel;
}
/*-----------------------------------------------*/
Cmplx sin(const Cmplx& z)
{
	 Cmplx s,c;

	 c = J*z;
	 s = J*(exp(c) - exp( -c ))*(-0.5);
	 return s;
}
/*-----------------------------------------------*/
Cmplx cos(const Cmplx& z)
{
	 Cmplx s,c;

	 c = J*z;
	 s = (exp(c) + exp( -c ))*0.5;
	 return s;
}
/*-----------------------------------------------*/
Cmplx cosh(const Cmplx& z)
{
	 Cmplx s;

	 s = ( exp(z) + exp(-z) )*0.5;
	 return s;
}
/*-----------------------------------------------*/
Cmplx sinh(const Cmplx& z)
{
	 Cmplx s;

	 s = ( exp(z) - exp(-z) )*0.5;
	 return s;
}
/*-----------------------------------------------*/
Cmplx exp(const Cmplx& z)
{
	Cmplx s;
	REAL re,im;

	re = z.real();
	im = z.imag();
	s  = ( J*sin(im) + cos(im) )*exp(re);
	return(s);
}
//------------------------------------------------
//Konstruktor
Cmplx::Cmplx(REAL x, REAL y)
{
	a = x;
	b = y;
}
//-------------------------------------------------
Cmplx::Cmplx(int aInt)
{
   a = aInt;
   b = 0.;
}
//-------------------------------------------------
Cmplx::Cmplx(const Cmplx& aZ)
{
   a = aZ.a;
   b = aZ.b;
}
//-------------------------------------------------
Cmplx::Cmplx(const Cmplx* pZ)
{
   a = pZ->a;
   b = pZ->b;
}
//-------------------------------------------------
bool Cmplx::IsReal()
{
  if(fabs(b) <= DBL_MIN)
                 return true;
  else
                 return false;
}
//-------------------------------------------------
bool Cmplx::IsImag()
{
  if(fabs(a) <= DBL_MIN)
                 return true;
  else
                 return false;
}
//-------------------------------------------------
int Cmplx::operator==(const Cmplx &z)
{
   if(a ==z.a && b == z.b)
                    return 1;
   else
                    return 0;

}
//-------------------------------------------------
int Cmplx::operator!=(const Cmplx& z)
{
   if(a != z.a) return 1;
   if(b != z.b) return 1;

   return 0;
}
//-------------------------------------------------
Cmplx Cmplx::operator-()const
{
  Cmplx res;
  res.a = -a;
  res.b = -b;

  return res;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator*=( const Cmplx& aZ)
{
	REAL c = a*aZ.a - b*aZ.b;
    REAL d = a*aZ.b + b*aZ.a;
    a      = c;
    b      = d;
	return *this;
}
/*-------------------------------------------------*/
Cmplx Cmplx::operator*(const Cmplx& aZ) // Komplexe Multiplikation
{
	Cmplx result = *this;

    result *= aZ;
	return result;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator/=( const Cmplx& aZ)
{
   Cmplx inv = aZ.invers();

   *this *= inv;

   return *this;
}
/*-------------------------------------------------*/
Cmplx Cmplx::operator/(const Cmplx& aZ) //Komplexe Division
{
  Cmplx result = *this;

  result /= aZ;

  return result;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator-=( const Cmplx& aZ)
{
	a -= aZ.a;
    b -= aZ.b;
	return *this;
}
/*---------------------------------------------------*/
Cmplx Cmplx::operator-(const Cmplx& aZ) //Subtraktion
{
	Cmplx result = *this;

	result -= aZ;

	return result;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator+=(const Cmplx& aZ)
{
	a += aZ.a;
    b += aZ.b;
	return *this;
}
/*----------------------------------------------------*/
Cmplx Cmplx::operator+(const Cmplx& aZ)  //Addition
{
	Cmplx result = *this;

	result += aZ;

	return result;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator*=(const REAL aR) // Reelle Multiplikation
{

	a *= aR;
	b *= aR;
	return *this;
}
/*-----------------------------------------------*/
Cmplx Cmplx::operator*(const REAL aR) // Reelle Multiplikation
{
	Cmplx s = *this;

     s *= aR;

	return s;
}
//---------------------------------------------------
Cmplx Cmplx::operator*=(const int aInt)
{
   return operator*=((REAL)aInt);
}
//---------------------------------------------------
Cmplx Cmplx::operator*(const int aInt)
{
   return operator*((REAL)aInt);
}
/*--------------------------------------------------*/
Cmplx Cmplx::operator/(const REAL aR) //Division durch reelle Zahl
{
   Cmplx res = *this;

   res /= aR;

   return res;
}
//----------------------------------------------------------
Cmplx Cmplx::operator/=(const REAL aR) //Division durch reelle Zahl
{
	REAL r, sgn;

	r = aR; sgn = 1;
	if(aR < 0)
	  {
	      r = -aR; sgn = -1;
	  }
	if( r <= DBL_MIN)
	  {
	     a = 0.7*FLT_MAX*sgn;
	     b = 0.7*FLT_MAX*sgn;
	     return *this;
	  }
	else
	  {
	     a /= aR;
	     b /= aR;
	     return *this;
	  }
}
//---------------------------------------------------
Cmplx Cmplx::operator/=(const int aInt)
{
   return operator/=((REAL)aInt);
}
//---------------------------------------------------
Cmplx Cmplx::operator/(const int aInt)
{
   return operator/((REAL)aInt);
}
//---------------------------------------------------
Cmplx Cmplx::operator=(const int aInt)
{
   return operator=((REAL)aInt);
}
/*-------------------------------------------------*/
Cmplx Cmplx::operator=(const REAL aR) //Gleichheitszeichen
{
	 a =  aR;
	 b =  0.;

	 return *this;
}
/*--------------------------------------------------*/
Cmplx Cmplx::operator-=(const REAL aR) //Subtraktion
{
	a -= aR;
	return *this;
}
//----------------------------------------------------
Cmplx Cmplx::operator-(const REAL aR) //Subtraktion
{
	Cmplx res = *this;

    res   -= aR;

	return  res;
}
//---------------------------------------------------
Cmplx Cmplx::operator-=(const int aInt)
{
   return operator-=((REAL)aInt);
}
//---------------------------------------------------
Cmplx Cmplx::operator-(const int aInt)
{
   return operator-((REAL)aInt);
}
/*--------------------------------------------------*/
Cmplx Cmplx::operator+=(const REAL aR) //Addition
{
	a += aR;
	return *this;
}
//----------------------------------------------------
Cmplx Cmplx::operator+(const REAL aR) //Addition
{
	Cmplx res = *this;

    res   += aR;

	return  res;
}
//---------------------------------------------------
Cmplx Cmplx::operator+=(const int aInt)
{
   return operator+=((REAL)aInt);
}
//---------------------------------------------------
Cmplx Cmplx::operator+(const int aInt)
{
   return operator+((REAL)aInt);
}
/*-----------------------------------------------------*/
Cmplx Cmplx::operator=(const Cmplx& aZ)  //Gleichsetzen
{
	 a = aZ.a;
	 b = aZ.b;
	 return *this;
}
/*-----------------------------------------------*/
Cmplx Cmplx::conjg(void) const
{
	 Cmplx s;

	 s.a = a;
	 s.b = -b;
	 return s;
}
/*-----------------------------------------------------*/
Cmplx Cmplx::invers(void)const
{
	 Cmplx s;
	 REAL r;

	 if( (r = abs()) <= DBL_MIN)
		{
		   s = FLT_MAX;
		   return s;
		}
	 else
		{
		   r   = r*r;
		   s.a = a/r;
		   s.b = -b/r;
		   return s;
		}
}
/*----------------------------------------------------*/
REAL Cmplx::winkel()const //Berechnet den Winkel der komplexen Zahl
{
	Cmplx z;
	int   sgn_r = 1, sgn_i = 1;
	REAL  wi = 0., re, im, kw=57.296;

	z = Cmplx(a,b);
	if(z.abs() == 0.)
	  return(0.);

	re = z.real();
	if( re< 0. )  sgn_r = -1;
	if( re==0. )  sgn_r =  0;

	im = z.imag();
	if( im< 0  )  sgn_i  = -1;
	if( im==0. )  sgn_i  =  0;

	switch(sgn_r )
	  {
	  case 1:                /* 1.  und  4. Quadrant */
		  wi = kw*atan(im/re);
		  break;
	  case -1:               /* 2. und 3. Quandrant  */
		  switch(sgn_i)
			 {
			 case 1:           /* 2. Quandrant         */
				 wi = 180. + kw*atan(im/re);
				  break;
			 case -1:          /* 3. Quandrant         */
				 wi = -180. + kw*atan(im/re);
				 break;
			 case 0:
				 wi = 180.;
				 break;
			 }
		  break;
	  case  0:
		  switch( sgn_i )
			 {
			 case 1:             /* 2. Quandrant   */
				 wi =  90.;
				 break;
			 case -1:            /* 3. Quandrant   */
				 wi = -90.;
				 break;
			 case 0:
				 wi =  0.;
				 break;
			 }
		  break;
	  }
 return(wi);
}
//-----------------------------------------------------
REAL Cmplx::MagTodB()
{
   REAL mag = this->abs();
   if(mag > 1.E-10) mag = 20*log10(mag);
   else
     mag               = -200.;
   return mag;
}
/*-----------------------------------------------------*/
void Cmplx::show(void)
{
	std::cout << '(' << a << ',' << b << ')';
}
/*------------------------------------------------------------*/
int gauss( Cmplx *loe_ma, Cmplx *in_ma, Cmplx *er_ma,
	   int dim,int nSol, int maxDim)
//-------------------------------------------------------------
/*    Funktion realisiert den Gauss-Algorithmus
		loe_ma   - Loesungsmatrix der Dimension n x nSol
		in_ma    - Matrix der Dimension n x n
		er_ma    - Erregungsmatrix der Dimension n x nSol
		dim      - Dimension der Matrix
		nSol     - Anzahl der L”sungen
		maxDim   - Dimension wie im rufenden Programm definiert

		Rueckgabewert 0 - Determinante verschwindet
			      1 - Erfolgreiche Loesung
		Uebergebene Matrix wird nicht veraendert  */
{
   int iZeil, nSpalt,jSol;
   Cmplx *ma = NULL;
   ma = new Cmplx[dim*(dim + nSol)];
   if(ma != NULL)
   { // Uebernehmen der neuen Matrix
     for( iZeil = 0; iZeil < dim; iZeil++)
       {
	for( nSpalt = 0; nSpalt < dim; nSpalt++)
	  {
	    ma[iZeil*(dim + nSol) + nSpalt] = in_ma[iZeil*maxDim + nSpalt];
	  }
       }
     for( jSol  = 0; jSol < nSol; jSol++)
       {
	for( iZeil = 0; iZeil < dim; iZeil++)
	  { //Uebernahme der Erregung
	   ma[iZeil*(dim + nSol) + dim + jSol]
			 = er_ma[ iZeil*nSol + jSol];
	  }
	}
    if(gauss(ma,dim,nSol))
      {
	for( jSol  = 0; jSol < nSol; jSol++)
	{
	for( iZeil = 0; iZeil < dim; iZeil++)
	 { //Uebernahme des Ergebnisses
	    loe_ma[iZeil*nSol + jSol]
			 = ma[iZeil*(dim + nSol) + dim +jSol];
	 }
	}
	delete ma;
	return 1;
      }
    delete ma;
   }
 return 0;
}
/*------------------------------------------------------------*/
int gauss( Cmplx *out_vek, Cmplx *in_ma, Cmplx *in_vek,
	   int dim, int maxDim)
//-------------------------------------------------------------
/*    Funktion realisiert den Gauss-Algorithmus
		out_vek  - Ergebnissvektor der Dimension n
		in_ma    - Matrix der Dimension n x n
		in_vek   - Erregungsvektor
		dim      - Dimension der Matrix
		maxDim   - Dimension wie im rufenden Programm definiert

		Rueckgabewert 0 - Determinante verschwindet
			      1 - Erfolgreiche Loesung
		Uebergebene Matrix wird nicht veraendert  */
{
   int iZeil, nSpalt;
   Cmplx *ma = NULL;
   ma = new Cmplx[dim*(dim + 1)];
   //Cmplx ma[6*9];
   if(ma != NULL)
   { // Uebernehmen der neuen Matrix
     for( iZeil = 0; iZeil < dim; iZeil++)
     {
	     for( nSpalt = 0; nSpalt < dim; nSpalt++)
		 {
	        ma[iZeil*(dim + 1) + nSpalt] = in_ma[iZeil*maxDim + nSpalt];
		 }
     }
     for( iZeil = 0; iZeil < dim; iZeil++)
     { //Uebernahme der Erregung
	     ma[iZeil*(dim + 1) + dim] = in_vek[iZeil];
     }
     if(gauss(ma,dim))
     {
	    for( iZeil = 0; iZeil < dim; iZeil++)
		{ //Uebernahme des Ergebnisses
	       out_vek[iZeil] = ma[iZeil*(dim + 1) + dim];
		}
      delete ma;
	  ma = NULL;
	  return 1;
	 }
    if(ma)
	{
	  delete ma;
	}
  }
 return 0;
}
/*------------------------------------------------------------*/
int gauss( Cmplx* ma, int n , int nSol)
/*------------------------------------------------------------*/
/*    Funktion realisiert den Gauss-Algorithmus
		ma   - Matrix der Dimension n x n+nSol
		n    - Dimension der Grundmatrix nxn
		nSol - Anzahl der gleichzeitigen Erregungen und
		       L”sungen

		Rueckgabewert 0 - Determinante verschwindet
			      1 - Erfolgreiche Loesung
		Uebergebene Matrix wird veraendert !!!          */
{
	int   i, k;
	int   kmax, j;
	int   m = n + nSol;
	REAL  smax;
	Cmplx pil, zz, det(1.,0.);

	for( i=0 ; i<n ; i++ ) /* Groesstes Element der Spalte bestimmen */
	  {
	  kmax = 0; smax = 0.;
	  for( k=i ; k<n ; k++ )
		 {
		 if( abs(ma[k*m+i]) > smax)
			{
			kmax = k; smax = abs(ma[k*m+i]);
			}
		 }
	  if( smax == 0)
		 {
		 //cout << "\n Funktion GAUSS: Determinante Null \n";
		 return(0);
		 }
	  else
		 {
		 if(kmax != i)
			{            /*   Zeilentausch  */
			for( j=i ; j<m ; j++ )
			  {
			  zz           = ma[i*m+j];
			  ma[i*m+j]    = ma[kmax*m+j];
			  ma[kmax*m+j] = zz;
			  }
			det = det*(-1);
			}
	  }

	/* Dividieren der Pilotzeile */
	pil = ma[i*m+i];
	det = pil*det;
	for( j=i+1 ; j<m ; j++ )
		ma[i*m+j] = ma[i*m+j]/pil;

	/* Eliminationsschritt       */
	for( k=i+1 ; k<n ; k++ )
		{
		zz = ma[k*m+i];
		for( j=i+1 ; j<n+1 ; j++ )
		  {
		  ma[k*m+j] = ma[k*m+j] - zz*ma[i*m+j];
		  }
		}
	}
	/*  Rueckwaertseinsetzen     */
	for( i=n-1 ; i>=0 ; i-- )
	  {
	  for( k=i-1 ; k>=0 ; k-- )
	    {
	     for(j = 0; j < nSol; j++)
	       {
		 ma[k*m+n+j] = ma[k*m+n+j] - ma[k*m+i]*ma[i*m+n+j];
	       }
	    }
	  }
	return(1);
}
/*------------------------------------------------------------*/
int gauss( Cmplx* ma, int n )
/*----------------------------------------------------*/
/*    Funktion realisiert den Gauss-Algorithmus
		ma - Matrix der Dimension n x n+1
		n  - Dimension der Matrix

		Rueckgabewert 0 - Determinante verschwindet
			      1 - Erfolgreiche Loesung
		Uebergebene Matrix wird veraendert !!!          */
{
	int   i, k;
	int   kmax, j, m = n + 1;
	REAL  smax;
	Cmplx pil, zz, det(1.,0);

	for( i=0 ; i<n ; i++ ) /* Groesstes Element der Spalte bestimmen */
	  {
	  kmax = 0; smax = 0.;
	  for( k=i ; k<n ; k++ )
		 {
		 if( abs(ma[k*m+i]) > smax)
			{
			kmax = k; smax = abs(ma[k*m+i]);
			}
		 }
	  if( smax == 0)
		 {
		 //cout << "\n Funktion GAUSS: Determinante Null \n";
		 return(0);
		 }
	  else
		 {
		 if(kmax != i)
			{            /*   Zeilentausch  */
			for( j=i ; j<m ; j++ )
			  {
			  zz           = ma[i*m+j];
			  ma[i*m+j]    = ma[kmax*m+j];
			  ma[kmax*m+j] = zz;
			  }
			det = det*(-1);
			}
	  }

	/* Dividieren der Pilotzeile */
	pil = ma[i*m+i];
	det = pil*det;
	for( j=i+1 ; j<m ; j++ )
		ma[i*m+j] = ma[i*m+j]/pil;

	/* Eliminationsschritt       */
	for( k=i+1 ; k<n ; k++ )
		{
		zz = ma[k*m+i];
		for( j=i+1 ; j<n+1 ; j++ )
		  {
		  ma[k*m+j] = ma[k*m+j] - zz*ma[i*m+j];
		  }
		}
	}
	/*  Rueckwaertseinsetzen     */
	for( i=n-1 ; i>=0 ; i-- )
	  {
	  for( k=i-1 ; k>=0 ; k-- )
		 {
		 ma[k*m+n] = ma[k*m+n] - ma[k*m+i]*ma[i*m+n];
		 }
	  }
	return(1);
}
