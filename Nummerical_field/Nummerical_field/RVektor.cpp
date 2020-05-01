#include "stdafx.h"
#include <iostream>


#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

#ifndef MYMATLIB_H
#include "mymatlib.h"
#endif
// ErrHdlr ist eine statische Elementfunktion
// zu einer void Funktion mit char* Argument
void (*Vektor<double>::ErrHdlr)(char*) = Vektor<double>::DefaultErrHdlr;
char Vektor<double>::ErrFlag = 0;
//---------------------------------------------------------
//Implementierung der Klasse RVektor
//---------------------------------------------------------
RVektor::RVektor(void) : Vektor<REAL>(1)
{
}
//---------------------------------------------------------
RVektor::RVektor(int n) : Vektor<REAL>(n)
{
}
//---------------------------------------------------------
void RVektor::SetZero()
{
	for(int i= 0; i < Size; i++)
		                    a[i] = 0;
}
//---------------------------------------------------------
REAL RVektor::Norm()
{
	REAL sum = 0;
	for(int i = 0; i < Size; i++)
		                   sum += a[i]*a[i];
	return sqrt(sum);
}
//---------------------------------------------------------
RVektor& RVektor::operator /=(const double aDouble)
{
	//ASSERT(fabs(aDouble) > 0);
	
	for(int i = 0; i < Size; i++)
	{
		a[i] /= aDouble;
	}

	return *this;
}
//---------------------------------------------------------
RVektor& RVektor::operator -=(const RVektor& aVec)
{
	//ASSERT(Size == aVec.Size);

	for(int i = 0; i < Size; i++)
	{
		a[i] -= aVec.a[i];
	}

	return *this;
}
//---------------------------------------------------------
RVektor RVektor::operator -(const RVektor& aVec)
{
	//ASSERT(Size == aVec.Size);

	RVektor res = *this;

	res -= aVec;

	return res;
}
//-----------------------------------------------------
RVektor& RVektor::operator*=(const double& r)
{
  for ( int i=0; i<Size; i++ )
    a[i] *= r;
  return *this;
}
//---------------------------------------------------------------
RVektor RVektor::operator*(const double& t) const
{
  RVektor Result = *this;
  Result *= t;
  return Result;
}
//---------------------------------------------------------------
RVektor operator*(const double& t, const RVektor& v)
{ 
	return v * t; 
}
//---------------------------------------------------------
RVektor& RVektor::operator+=(const RVektor& v)
{
  if ( v.Size != Size )
    Error( "Vektoren verschiedener Größen können nicht addiert werden" );
  else
    for ( int i=0; i<Size; i++ )
      a[i] += v.a[i];
  return *this;
}
//---------------------------------------------------------------
RVektor RVektor::operator+(const RVektor& v) const
{
  RVektor Result = *this;
  Result += v;
  return Result;
}
//-------------------------------------------------------------
RVektor RVektor::operator-(const RVektor& v) const
{
  RVektor Result = *this;
  Result -= v;
  return Result;
}
//----------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const RVektor& v)
{
  //os << '(';						comment out for compatibility with .m format
  int n = v.nElems() - 1;
  int i;
  for (i=0; i<n; i++ )
          os << v.Elem(i) << ' ';
  os << v.Elem(i) ;					//reduce one ')' for compatibility with .m format
  return os;
}
//------------------------------------------------------------
std::istream& operator>>(std::istream& is, RVektor& v)
{
  char c;
  int n = v.nElems();
  is >> c;  // lese '('
  for (int i=0; i<n; i++ )
    is >> v[i];
  is >> c;  // lese ')'
  return is;
}

#ifndef MY_KONSOL_APP
void RVektor::Serialize(CArchive& ar)
{
	if(ar.IsStoring())
	{
		ar << Size;
		for(int i = 0; i < Size; i++)
		{
			ar << a[i];
		}
	}
	else
	{
		ar >> Size;
		if(a)
		{
			delete a;
			a = NULL;
		}
		a = new REAL[Size];
		if(a)
		{
			for(int i = 0; i < Size; i++)
			{
				ar >> a[i];
			}
		}
		else
		   Size = 0;
	}
}
#endif
